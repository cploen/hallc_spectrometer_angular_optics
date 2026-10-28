#!/usr/bin/env python3
"""Apply an HMS matrix to saved TFit events and run the legacy residual macros."""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

import numpy as np
import uproot

from compare_matrices import read_matrix
from core_sample import PROJECT, digest, read_tsv
from sieve_afterburner import reconstruct
from spectrometer_config import from_campaign, hms_xptar_truth, run_metadata

REQUIRED = set('entry xfp xpfp yfp ypfp reactycalc ztar xtar xptar yptar ytar '
               'xs ys xsT ysT xptarT yptarT ytarT ztarT delta'.split())
UPDATED = ('xtar', 'xptar', 'ytar', 'yptar', 'xs', 'ys')


def prepare(campaign, fit_dir, metadata=None, fileid=-1):
    if from_campaign(campaign).name != 'HMS':
        raise ValueError('The legacy afterburner adapter currently supports HMS only')
    tables = list((campaign / 'config').glob('rungroups_*_inputs.tsv'))
    if len(tables) != 1:
        raise ValueError('Expected one campaign rungroup table')
    candidates = [fit_dir.parent / 'optics.dat', campaign / 'config/optics.dat',
                  PROJECT / 'DATfiles/list_of_optics_run.dat']
    metadata = Path(metadata) if metadata else next(p for p in candidates if p.is_file())
    empty = set()
    build = fit_dir.parent / 'build.json'
    if build.exists():
        empty = {r['rungroup'] for r in json.loads(build.read_text()).get('groups', [])
                 if r.get('status') == 'empty' and r.get('entries') == 0}
    groups = []
    for row in read_tsv(tables[0]):
        name = row['rungroup']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name):
            raise ValueError('Invalid rungroup name')
        if name in empty:
            continue
        path = fit_dir / f'Optics_{int(row["optics_id"])}_{fileid}_fit_tree_gmm.root'
        meta = run_metadata(int(row['optics_id']), metadata)
        if not np.isclose(meta['angle_deg'], float(row['hms_angle_deg']), rtol=0, atol=1e-6):
            raise ValueError(f'{name}: campaign/metadata angle mismatch')
        with uproot.open(path) as root:
            tree = root['TFit']
            missing = REQUIRED - set(tree.keys())
            if missing:
                raise ValueError(f'{name}: missing TFit branches {sorted(missing)}; rebuild fit ntuples')
            entries = tree.num_entries
        if entries:
            groups.append(dict(rungroup=name, optics_id=int(row['optics_id']),
                               path=path, angle=meta['angle_deg'], entries=entries))
    if not groups:
        raise ValueError('No nonempty fit trees')
    if len({r['rungroup'] for r in groups}) != len(groups):
        raise ValueError('Duplicate rungroups')
    return groups, metadata, tables[0]


def export_tree(source, target, rows, angle, offsets, xmis=None):
    """Keep membership/truth, using reactycalc and ztar saved by the TFit builder."""
    total = unfinished = 0
    with uproot.open(source) as root, uproot.recreate(target) as output:
        tree = root['TFit']
        for a in tree.iterate(step_size=50000, library='np'):
            if any(v.ndim != 1 or v.dtype.kind not in 'fiub' for v in a.values()):
                raise ValueError('Expected scalar numeric TFit branches')
            if any(not np.isfinite(a[k]).all() for k in REQUIRED):
                raise ValueError('Nonfinite TFit inputs; no events silently discarded')
            if a['entry'].dtype.kind not in 'iu':
                raise ValueError('TFit entry IDs must be integers')
            expected = hms_xptar_truth(a['xsT'], a['ztarT'], angle, a['reactycalc'])
            if not np.allclose(a['xptarT'], expected, rtol=1e-12, atol=1e-13):
                raise ValueError('Stale or inconsistent HMS xptarT; rebuild fit ntuples with corrected geometry')
            if any('replay_' + k in a for k in UPDATED):
                raise ValueError('Input already contains afterburner replay_* branches; use the original TFit')
            inputs = {k: a[k] for k in ('xfp', 'xpfp', 'yfp', 'ypfp')}
            inputs.update(ry=a['reactycalc'], rz=a['ztar'])
            xy, iterations, active, values, xtar = reconstruct(
                inputs, rows, angle, offsets, xmis=xmis, return_target=True)
            for k in UPDATED:
                a['replay_' + k] = a[k].copy()
            a.update(xtar=xtar, xptar=values[:, 0], ytar=100 * values[:, 1],
                     yptar=values[:, 2], xs=xy[:, 0], ys=xy[:, 1],
                     afterburner_iterations=iterations.astype(np.int32))
            if total == 0:
                output.mktree('TFit', {k: v.dtype for k, v in a.items()})
            output['TFit'].extend(a)
            total += len(xtar)
            unfinished += int(active.sum())
        output['afterburner_convention'] = (
            'Positions cm; slopes rad; complete matrix including 00000; external offsets '
            + json.dumps(list(offsets)) + ' [xptar mrad, ytar cm, yptar mrad]. '
            'Saved reactycalc/ztar vertex, original delta/truth/event selection retained.')
    return dict(entries=total, iteration_limit=unfinished)


def plot_commands(groups, out, metadata, executable):
    for group in groups:
        tag = group['rungroup']
        for variable in ('xptar', 'yptar', 'ytar'):
            macro = PROJECT / f'diagnostics/fit/plot_{variable}_residuals.C'
            args = [str(group['optics_id']), json.dumps(str(out / group['output'])),
                    json.dumps(str(metadata)), json.dumps(str(out / 'residuals' / tag)),
                    json.dumps(tag)]
            yield group, variable, [executable, '-l', '-b', '-q', str(macro) + '(' + ','.join(args) + ')']


def plot(groups, out, metadata, executable, save):
    (out / 'logs').mkdir(exist_ok=True)
    for group, variable, cmd in plot_commands(groups, out, metadata, executable):
        tag = group['rungroup']
        print(shlex.join(cmd), flush=True)
        log = out / 'logs' / f'{tag}_{variable}.log'
        with log.open('w') as stream:
            result = subprocess.run(cmd, cwd=PROJECT, stdout=stream, stderr=subprocess.STDOUT)
        folder = out / 'residuals' / tag
        hist = folder / f'{tag}_{variable}_residuals.root'
        expected = [folder / f'{tag}_{variable}_residuals.pdf',
                    folder / f'{tag}_{variable}_residuals_by_sieve.pdf']
        if variable == 'ytar':
            expected.append(folder / f'{tag}_ytar_residuals_delta.pdf')
        if result.returncode or not hist.is_file() or any(not p.is_file() for p in expected):
            raise RuntimeError(f'Legacy {variable} diagnostics failed; inspect {log}')
        with uproot.open(hist) as root:
            if not root.keys():
                raise RuntimeError(f'Empty histogram output; inspect {log}')
        save()


def run(campaign, name, matrix, fit_dir, offsets=(0., 0., 0.), metadata=None,
        fileid=-1, xmis=None, export_only=False, dry_run=False):
    groups, metadata, table = prepare(campaign, fit_dir, metadata, fileid)
    rows = read_matrix(matrix)
    out = campaign / '07_diagnostics' / name
    if out.exists():
        raise FileExistsError(f'Output exists; choose a new name: {out}')
    executable = shutil.which(os.environ.get('ROOT', 'root'))
    if not export_only and not dry_run and not executable:
        raise ValueError('Load the ROOT environment or set ROOT=/path/to/root')
    for group in groups:
        group['output'] = f'root/{group["rungroup"]}.root'
    if dry_run:
        print(f'Matrix: {matrix}; external offsets: {list(offsets)}; metadata: {metadata}')
        for g in groups:
            print(f'Export {g["path"]} -> {out / g["output"]} ({g["entries"]} events)')
        if not export_only:
            for _, _, cmd in plot_commands(groups, out, out / 'optics.dat', executable or 'root'):
                print(shlex.join(cmd))
        return
    for sub in ('root', 'logs'):
        (out / sub).mkdir(parents=True, exist_ok=True)
    manifest = dict(schema='legacy_afterburner_v1', status='running', matrix=str(matrix),
                    matrix_sha256=digest(matrix), external_offsets=list(offsets), xmis_cm=xmis,
                    offset_units=['xptar mrad', 'ytar cm', 'yptar mrad'],
                    fitted_constants='All matrix terms including 00000 retained',
                    reconstruction='Fixed saved reaction vertex; iterative xtar; original delta and truth retained',
                    sample=str(fit_dir), metadata_sha256=digest(metadata), table_sha256=digest(table), groups=[])
    def save():
        (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    save()
    try:
        shutil.copy2(matrix, out / 'matrix.dat')
        shutil.copy2(metadata, out / 'optics.dat')
        shutil.copy2(table, out / table.name)
        for group in groups:
            print(f'Apply matrix: {group["rungroup"]} ({group["entries"]} events)', flush=True)
            report = export_tree(group['path'], out / group['output'], rows,
                                 group['angle'], offsets, xmis)
            if report['entries'] != group['entries']:
                raise RuntimeError('Export membership count changed')
            manifest['groups'].append(dict(rungroup=group['rungroup'], source=str(group['path']),
                source_sha256=digest(group['path']), output=group['output'], **report))
            if report['iteration_limit']:
                print(f'WARNING: {report["iteration_limit"]} events reached the xtar iteration limit', flush=True)
            save()
        if not export_only:
            plot(groups, out, out / 'optics.dat', executable, save)
        manifest['status'] = 'exported' if export_only else 'complete'
        manifest['outputs'] = {str(p.relative_to(out)): digest(p) for p in out.rglob('*')
                               if p.is_file() and p.name != 'manifest.json'}
        save()
        print(f'Complete: {out}', flush=True)
    except BaseException as exc:
        manifest.update(status='incomplete', error=str(exc))
        save()
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign', type=Path)
    p.add_argument('name', nargs='?', default='legacy_afterburner')
    p.add_argument('--matrix', required=True, type=Path)
    p.add_argument('--fit-dir', required=True, type=Path, help='Original TFit root/ directory')
    p.add_argument('--metadata', type=Path)
    p.add_argument('--offsets', nargs=3, type=float, default=[0., 0., 0.], metavar=('XP_MRAD', 'Y_CM', 'YP_MRAD'))
    p.add_argument('--xmis', type=float, help='Override HMS vertical mispointing in cm')
    p.add_argument('--fileid', type=int, default=-1)
    p.add_argument('--export-only', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    try:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', a.name):
            raise ValueError('Invalid output name')
        if not np.isfinite(a.offsets).all() or (a.xmis is not None and not np.isfinite(a.xmis)):
            raise ValueError('Offsets and xmis must be finite')
        run(a.campaign.resolve(), a.name, a.matrix.resolve(), a.fit_dir.resolve(),
            a.offsets, a.metadata, a.fileid, a.xmis, a.export_only, a.dry_run)
    except (OSError, ValueError, KeyError, RuntimeError, StopIteration) as exc:
        p.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
