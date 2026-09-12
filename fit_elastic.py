#!/usr/bin/env python3
"""Fit balanced HMS cores; evaluate protected samples only with --evaluate."""
import argparse
import importlib.metadata
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path
import numpy as np
import uproot
from core_sample import PROJECT, digest, read_tsv, write_tsv
from spectrometer_config import from_campaign, run_metadata
from preallocated_svd import verify_build
from elastic_net import DEFAULTS, TARGETS, SCALES, validate, folds_for, tune


def matrix_rows(path):
    rows, inside = [], False
    for line in path.read_text().splitlines():
        if line.strip().startswith('----'):
            if inside: break
            inside = True
            continue
        if not inside or not line.strip() or line.lstrip().startswith('!'): continue
        fields = line.split()
        if len(fields) != 5 or not re.fullmatch(r'\d{5}', fields[4]):
            raise ValueError(f'Malformed transport row: {line}')
        coeff = np.array([float(v.replace('D', 'E').replace('d', 'e')) for v in fields[:4]])
        if not np.isfinite(coeff).all(): raise ValueError('Nonfinite matrix coefficient')
        rows.append((coeff, tuple(map(int, fields[4]))))
    if not rows or len({t for c, t in rows}) != len(rows):
        raise ValueError('Missing or duplicate matrix terms')
    return rows


def basis(rows):
    return [(0, 0, 0, 0, 0)] + [t for c, t in rows if t[4] == 0 and any(t)]


def design(a, terms):
    variables = [a['xfp']/100., a['xpfp'], a['yfp']/100., a['ypfp'], a['xtar']/100.]
    x = np.ones((len(a['entry']), len(terms)))
    for j, term in enumerate(terms):
        for v, power in zip(variables, term):
            if power: x[:, j] *= v**power
    if not np.isfinite(x).all(): raise ValueError('Nonfinite design matrix')
    return x


def problem(a, rows):
    fixed = [(c, t) for c, t in rows if t[4]]
    offset = design(a, [t for c, t in fixed]) @ np.array([c[:3] for c, t in fixed]).reshape(-1, 3)
    truth = np.column_stack((a['xptarT'], a['ytarT']/100., a['yptarT']))
    return design(a, basis(rows)), truth-offset, offset


def write_matrix(path, rows, terms, coeff):
    old = {t: c for c, t in rows}
    with path.open('w') as out:
        out.write('! Candidate angular fit; delta and xtar-dependent terms retained from frozen seed.\n --------------------------------------------\n')
        for t, c in zip(terms, coeff):
            values = list(c) + [old.get(t, np.zeros(4))[3]]
            out.write(' '+' '.join(f'{v:.16e}' for v in values)+' '+''.join(map(str, t))+'\n')
        for c, t in rows:
            if t[4]: out.write(' '+' '.join(f'{v:.16e}' for v in c)+' '+''.join(map(str, t))+'\n')
        out.write(' --------------------------------------------\n')


def checked(path, sha, hashes):
    actual = digest(path)
    if actual != sha: raise ValueError(f'Frozen input changed: {path}')
    hashes[str(path)] = actual


def load_sample(campaign, tag, sample):
    """Exact entry join against the allocation. No GMM fallback or resampling."""
    source = campaign/'05c_core_sample'/tag
    build = campaign/'06c_core_ntuple'/tag/sample
    hashes = {}
    manifest = json.loads((source/'manifest.json').read_text())
    if (manifest.get('schema'), manifest.get('mode'), manifest.get('status')) != ('core_balance_v1', 'event_allocation', 'quotas_ready'):
        raise ValueError('Requires a completed event-level balanced allocation')
    bm = json.loads((build/'build.json').read_text())
    checked(source/'manifest.json', bm['sample_manifest'], hashes)
    hashes[str(build/'build.json')] = digest(build/'build.json')
    if bm['sample'] != sample: raise ValueError('Wrong sample build')
    if sample == 'fit':
        _, wanted = verify_build(build)
    for rel, h in bm.get('outputs', {}).items(): checked(build/rel, h, hashes)
    tables = list(source.glob('rungroups_*_inputs.tsv'))
    if len(tables) != 1: raise ValueError('Expected one frozen campaign table')
    for p in (tables[0], source/'metadata/optics.dat', source/'metadata/sieve_mask.json'):
        checked(p, manifest['outputs'][str(p.relative_to(source))], hashes)
    spec = from_campaign(campaign)
    if spec.name != 'HMS': raise ValueError('Elastic TFit adapter currently supports HMS only')
    blocked = {tuple(v) for v in json.loads((source/'metadata/sieve_mask.json').read_text())['blocked']}
    code = {'fit': 1, 'holdout': 2}[sample]
    pieces = []
    settings = read_tsv(tables[0])
    if len({r['rungroup'] for r in settings}) != len(settings):
        raise ValueError('Duplicate setting in frozen table')
    edges = None
    for row in sorted(settings, key=lambda r: r['rungroup']):
        name = row['rungroup']
        maskpath = source/f'root/CoreSample_{name}.root'
        checked(maskpath, manifest['outputs'][str(maskpath.relative_to(source))], hashes)
        with uproot.open(maskpath) as f: mask = f['CoreSample'].arrays(library='np')
        idx = np.flatnonzero(mask['sample'] == code)
        idx = idx[np.argsort(mask['entry'][idx])]
        if not len(idx): continue
        p = build/f'root/Optics_{row["optics_id"]}_-1_fit_tree_gmm.root'
        rel = str(p.relative_to(build))
        if rel not in bm.get('outputs', {}):
            raise ValueError(f'Build lacks checksums for {p}; rebuild with current run_build_core_fit.sh')
        with uproot.open(p) as f: a = f['TFit'].arrays(library='np')
        order = np.argsort(a['entry'])
        a = {k: v[order] for k, v in a.items()}
        if a['entry'].dtype.kind not in 'iu' or len(np.unique(a['entry'])) != len(a['entry']):
            raise ValueError('Event IDs must be unique integers')
        if not np.array_equal(a['entry'], mask['entry'][idx]):
            raise ValueError(f'{name}: TFit membership differs from balanced {sample} mask')
        for k in ('foil', 'ndel', 'xscol', 'yscol', 'core_keep', 'sample'):
            if not np.array_equal(a[k], mask[k][idx]): raise ValueError(f'{name}: {k} mismatch')
        for field, key in [('xs', 'xsieve'), ('ys', 'ysieve'), ('ztarT', 'zfoil')]:
            if not np.allclose(a[field], mask[key][idx], rtol=0, atol=1e-8):
                raise ValueError(f'{name}: {field} changed')
        meta = run_metadata(int(row['optics_id']), source/'metadata/optics.dat')
        if edges is not None and meta['edges'] != edges: raise ValueError('Incompatible delta boundaries')
        edges = meta['edges']
        if np.any(a['foil'] < 0) or np.any(a['foil'] >= len(meta['foils'])) or np.any(a['ndel'] < 0) or np.any(a['ndel'] >= len(meta['edges'])-1):
            raise ValueError('Foil/delta index outside saved metadata')
        if np.any(a['xscol'] < 0) or np.any(a['xscol'] >= spec.nx) or np.any(a['yscol'] < 0) or np.any(a['yscol'] >= spec.ny):
            raise ValueError('Hole outside geometry')
        lo, hi = np.array(meta['edges'])[a['ndel']], np.array(meta['edges'])[a['ndel']+1]
        if not (np.all(a['delta'] >= lo) and np.all(a['delta'] < hi)):
            raise ValueError('TFit delta disagrees with saved slice')
        for k, values in [('ztarT', np.array(meta['foils'])[a['foil']]),
                          ('xsT', np.array([spec.xs(int(i)) for i in a['xscol']])),
                          ('ysT', np.array([spec.ys(int(i)) for i in a['yscol']]))]:
            if not np.allclose(a[k], values, rtol=0, atol=1e-8): raise ValueError(f'{name}: truth geometry mismatch')
        a['quality'] = mask['quality'][idx]
        if not np.isin(a['quality'], [0, 1, 2]).all(): raise ValueError('Unknown quality code')
        a['excluded'] = np.array([(int(x), int(y)) in blocked for x, y in zip(a['xscol'], a['yscol'])])
        if sample == 'fit' and np.any((a['quality'] != 2) | (a['core_keep'] != 1) | a['excluded']):
            raise ValueError('Training contains noncore or blocked labels')
        a['rungroup'] = np.full(len(idx), name)
        a['angle'] = np.full(len(idx), np.deg2rad(meta['angle_deg']))
        a['delta_low'], a['delta_high'] = lo, hi
        fields = ('entry', 'foil', 'ndel', 'xscol', 'yscol', 'xs', 'ys', 'xsT', 'ysT', 'ztarT',
                  'xfp', 'xpfp', 'yfp', 'ypfp', 'xtar', 'xptarT', 'ytarT', 'yptarT',
                  'delta', 'quality', 'excluded', 'rungroup', 'angle', 'delta_low', 'delta_high')
        for k in fields:
            if k != 'rungroup' and not np.isfinite(a[k]).all(): raise ValueError(f'Nonfinite {name}/{k}; no silent dropping')
        pieces.append({k: a[k] for k in fields})
    if not pieces: raise ValueError(f'Empty {sample} pool')
    data = {k: np.concatenate([a[k] for a in pieces]) for k in pieces[0]}
    if sample == 'fit':
        actual = set(zip(data['rungroup'], map(int, data['entry'])))
        if actual != {(r['rungroup'], int(r['entry'])) for r in wanted}: raise ValueError('Selected-ID mismatch')
        _, counts = np.unique(data['ztarT'], return_counts=True)
        if len(set(counts)) != 1: raise ValueError('Physical foil totals are not strictly equal')
    return data, hashes


def missing_inputs(campaign, tag):
    source = campaign/'05c_core_sample'/tag
    needed = [source/'manifest.json', campaign/'06c_core_ntuple'/tag/'fit/build.json',
              campaign/'06c_core_ntuple'/tag/'holdout/build.json']
    return [str(p) for p in needed if not p.is_file()]


def fit(campaign, tag, name):
    out = campaign/'06d_elastic_net'/name
    if out.exists(): raise FileExistsError(f'Immutable output exists: {out}')
    a, hashes = load_sample(campaign, tag, 'fit')
    seed = campaign/'06c_core_ntuple'/tag/'fit/oldfit.dat'
    rows = matrix_rows(seed)
    cfg = dict(DEFAULTS)
    policy = campaign/'config/elastic.json'
    config_source = None
    if policy.exists():
        cfg.update(json.loads(policy.read_text()))
        config_source = dict(path=str(policy), sha256=digest(policy))
    validate(cfg)
    x, y, offset = problem(a, rows)
    ids = [f'{r}:{int(e)}' for r, e in zip(a['rungroup'], a['entry'])]
    strata = [f'{r}:{f}:{d}:{h}:{v}' for r, f, d, h, v in zip(a['rungroup'], a['foil'], a['ndel'], a['xscol'], a['yscol'])]
    folds = folds_for(ids, strata, cfg['folds'], cfg['seed'])
    cells = np.array([f'{z}:{d}' for z, d in zip(a['ztarT'], a['ndel'])])
    result = tune(x[:, 1:], y, folds, cells, cfg)
    out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.elastic_', dir=out.parent))
    try:
        for folder in ('matrices', 'plots', 'tsv', 'code'): (stage/folder).mkdir()
        shutil.copy2(seed, stage/'seed.dat')
        for model in ('enet', 'svd'):
            write_matrix(stage/f'matrices/{model}.dat', rows, basis(rows), result[model])
            reread = matrix_rows(stage/f'matrices/{model}.dat')
            pred = design(a, [t for c, t in reread]) @ np.array([c[:3] for c, t in reread])
            if not np.allclose(pred, x @ result[model]+offset, rtol=1e-10, atol=1e-11):
                raise ValueError('Exported matrix does not reproduce predictions')
        np.savez_compressed(stage/'model.npz', enet=result['enet'], svd=result['svd'],
                            terms=basis(rows), mean=result['state'][0], scale=result['state'][1],
                            target_mean=result['state'][2], target_scale=result['state'][3])
        write_tsv(stage/'tsv/folds.tsv', [dict(rungroup=r, entry=int(e), fold=int(f)) for r, e, f in zip(a['rungroup'], a['entry'], folds)])
        write_tsv(stage/'tsv/cv.tsv', result['cv'])
        write_tsv(stage/'tsv/cv_summary.tsv', result['summary'])
        term_rows = []
        for i, t in enumerate(basis(rows)[1:]):
            for j, target in enumerate(TARGETS):
                term_rows.append(dict(term=''.join(map(str, t)), target=target,
                    coefficient=float(result['enet'][i+1, j]), standardized=float(result['beta'][i, j]),
                    active=bool(result['beta'][i, j] != 0), fold_frequency=float(result['frequency'][i, j])))
        write_tsv(stage/'tsv/terms.tsv', term_rows)
        from elastic_diagnostics import training_plots, metrics
        training_plots(stage, result)
        write_tsv(stage/'tsv/training.tsv', metrics(a, y, {m: x @ result[m] for m in ('enet', 'svd')}, np.full(len(y), 'training'), cfg['qa_min'], offset))
        for file in ('fit_elastic.py', 'elastic_net.py', 'elastic_diagnostics.py', 'run_elastic.sh',
                     'spectrometer_config.py', 'spectrometer_profiles.def', 'core_sample.py', 'preallocated_svd.py'):
            shutil.copy2(PROJECT/file, stage/'code'/file)
        report = dict(schema='elastic_v1', campaign=campaign.name, sample=tag, config=cfg, config_source=config_source,
                      training=len(y), choices=result['choices'], svd=result['svd_info'], inputs=hashes,
                      versions={m: importlib.metadata.version(m) for m in ('numpy', 'scipy', 'scikit-learn', 'uproot', 'matplotlib')})
        (stage/'RESULTS.md').write_text(report_text(report))
        report['outputs'] = {str(p.relative_to(stage)): digest(p) for p in stage.rglob('*') if p.is_file()}
        (stage/'manifest.json').write_text(json.dumps(report, indent=2)+'\n')
        stage.rename(out)
    except BaseException:
        shutil.rmtree(stage)
        raise
    print(f'Fit complete: {out}. Protected holdouts have not been evaluated.')
    return out


def report_text(report):
    text = '# Elastic-net angular fit: '+report['campaign']+'\n\n'
    text += f"{report['training']:,} balanced training cores; {report['config']['folds']} training-only validation folds.\n\n"
    text += '| Target | Alpha | L1 ratio | Active slopes | CV MSE / scaled SVD | Alpha at grid edge |\n|---|---:|---:|---:|---:|---|\n'
    for r in report['choices']:
        ratio = r['mse']/r['svd_cv_mse'] if r['svd_cv_mse'] else float('nan')
        text += f"| {r['target']} | {r['alpha']:g} | {r['l1']:g} | {r['final_active']} | {ratio:.3f} | {r['alpha_edge']} |\n"
    text += '\nScores average squared error equally over populated physical-foil/delta cells within each fold. '
    text += 'The one-standard-error rule favors fewer terms; its error bar is a fold-variability heuristic, not a confidence interval.\n\n'
    text += 'Compare [validation curves](plots/cv.png), [term stability](plots/terms.png), '
    text += '[training strata](tsv/training.tsv), and [term candidates](tsv/terms.tsv). '
    text += 'The scaled SVD baseline uses identical events, basis, fixed terms and training-only transformations; it is not the historical ROOT normal-equation fit.\n\n'
    text += 'Beam search is a possible development step, not a demonstrated improvement. Use fold frequency and standardized coefficients as candidate rankings, not proof of physical importance. '
    text += 'Keep the constant and fixed xtar/delta terms. Use these saved training folds, recomputing scaling in each fold. '
    text += 'Do not rank terms or tune search settings on protected core/noncore results. Prefer completing any beam-search study before the one-time protected evaluation. '
    text += 'If protected diagnostics subsequently guide changes, they become development evidence and are no longer an untouched final test.\n\n'
    text += 'Do not advance a matrix based on pooled RMS alone: inspect edge-hole counts/residuals and all foil/delta cells. '
    text += 'Grid-edge selections, unstable term membership, or little CV advantage over SVD call for further training-only study. '
    text += 'Candidate matrices have not been installed in replay. Reconstructed quantities here keep the saved replay xtar input fixed; full replay/iteration validation remains necessary.\n'
    return text


def evaluate(campaign, name):
    out = campaign/'06d_elastic_net'/name
    report = json.loads((out/'manifest.json').read_text())
    hashes = {}
    for rel, sha in report['outputs'].items(): checked(out/rel, sha, hashes)
    for path, sha in report['inputs'].items(): checked(Path(path), sha, hashes)
    target = out/'evaluation'
    if target.exists(): raise FileExistsError('Protected evaluation already exists; do not retune against it')
    a, inputs = load_sample(campaign, report['sample'], 'holdout')
    train_ids = {(r['rungroup'], int(r['entry'])) for r in read_tsv(out/'tsv/folds.tsv')}
    if train_ids & set(zip(a['rungroup'], map(int, a['entry']))): raise ValueError('Training/holdout overlap')
    rows = matrix_rows(out/'seed.dat')
    x, y, offset = problem(a, rows)
    with np.load(out/'model.npz') as saved:
        if not np.array_equal(saved['terms'], basis(rows)): raise ValueError('Saved basis mismatch')
        predictions = {m: x @ saved[m] for m in ('enet', 'svd')}
    # Blocked labels are reported but cannot define valid target residuals.
    pool = np.where(a['excluded'], 'blocked', np.where(a['quality'] == 2, 'core', np.where(a['quality'] == 1, 'noncore', 'unsupported')))
    stage = Path(tempfile.mkdtemp(prefix='.evaluation_', dir=out))
    try:
        from elastic_diagnostics import evaluation_plots, metrics
        (stage/'plots').mkdir(); (stage/'tsv').mkdir()
        write_tsv(stage/'tsv/metrics.tsv', metrics(a, y, predictions, pool, report['config']['qa_min'], offset))
        evaluation_plots(stage, a, y, predictions, offset, pool, from_campaign(campaign), report['config']['qa_min'])
        write_tsv(stage/'tsv/ids.tsv', [dict(rungroup=r, entry=int(e), pool=p) for r, e, p in zip(a['rungroup'], a['entry'], pool)])
        counts = {p: int(sum(pool == p)) for p in ('core', 'noncore', 'unsupported', 'blocked')}
        (stage/'README.md').write_text('# Protected evaluation\n\n'+str(counts)+'\n\n'
            'Core and noncore are scored separately. Unsupported labels are a separate stress test; blocked labels have counts only. '
            'Noncore residuals assume the inherited hole/foil labels are valid; tails can include scattered or mislabeled events. '
            'Plots share axes/scales within comparisons. Empty holes have no residual estimate, not zero error. '
            'Sparse cells are marked by N and qa_low in metrics.tsv; QA never changes membership.\n\n'
            'Review residuals.png, foil_delta.png and all slice_* figures; exact bias, RMS and P95 absolute residuals are in tsv/metrics.tsv. '
            'Clouds show reconstructed sieve coordinates from the candidate matrix with saved replay xtar; they are not a new replay. '
            'The ztar diagnostic uses the HMS target identity and beam position inferred from the exported truth geometry. '
            'No improvement or beam-search readiness is certified automatically.\n')
        m = dict(fit_manifest=digest(out/'manifest.json'), inputs=inputs, counts=counts,
                 diagnostic_code=digest(PROJECT/'elastic_diagnostics.py'), adapter_code=digest(PROJECT/'fit_elastic.py'))
        m['outputs'] = {str(p.relative_to(stage)): digest(p) for p in stage.rglob('*') if p.is_file()}
        (stage/'manifest.json').write_text(json.dumps(m, indent=2)+'\n')
        stage.rename(target)
    except BaseException:
        shutil.rmtree(stage)
        raise
    print(f'Protected evaluation complete: {target}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign', type=Path)
    p.add_argument('tag', nargs='?', default='equal15')
    p.add_argument('name', nargs='?', default='enet')
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--check', action='store_true', help='Verify required exports; no outputs or fit')
    mode.add_argument('--evaluate', action='store_true', help='Evaluate an existing frozen fit once on protected pools')
    a = p.parse_args()
    for label in (a.tag, a.name):
        if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*', label): p.error('Invalid tag/name')
    campaign = (a.campaign if a.campaign.is_absolute() else PROJECT/a.campaign).resolve()
    try:
        if from_campaign(campaign).name != 'HMS': raise ValueError('HMS campaigns only in this adapter')
        if a.check:
            missing = missing_inputs(campaign, a.tag)
            if missing: raise FileNotFoundError('Required campaign exports are missing:\n'+'\n'.join(missing))
            for sample in ('fit', 'holdout'):
                data, _ = load_sample(campaign, a.tag, sample)
                print(f'OK {sample}: {len(data["entry"]):,} exact-membership events')
        elif a.evaluate:
            report = json.loads((campaign/'06d_elastic_net'/a.name/'manifest.json').read_text())
            if report['sample'] != a.tag: raise ValueError('Requested tag differs from saved fit')
            evaluate(campaign, a.name)
        else: fit(campaign, a.tag, a.name)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        p.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__': main()
