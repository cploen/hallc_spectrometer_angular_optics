#!/usr/bin/env python3
"""Export small, traceable production results for the NPS angular-fit slides."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

PROJECT = Path(__file__).resolve().parent
CAMPAIGNS = (
    ('HMS_5p878GeV', 'inner10', 16.500),
    ('HMS_6p117GeV', 'inner20', 12.385),
    ('HMS_6p667GeV_joint_min5', 'inner20', 12.490),
)


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export(output, make_sieve=False):
    policy = json.loads((PROJECT/'historical_offsets_2024.json').read_text())
    files = {}
    campaigns = []

    def add(path):
        if not path.is_file():
            raise FileNotFoundError(path)
        files[str(path.relative_to(PROJECT))] = path

    add(PROJECT/'historical_offsets_2024.json')
    for name, tag, angle in CAMPAIGNS:
        campaign = PROJECT/name
        beam = campaign/'06e_beam_search'/f'beam20_{tag}'
        standard = campaign/'07_diagnostics'/f'beam20_{tag}_standard_2024flags'
        manifest = json.loads((standard/'manifest.json').read_text())
        if manifest.get('status') != 'complete':
            raise ValueError(f'{standard}: diagnostics are not complete')
        offsets = policy['campaigns'][name.split('_joint')[0]]['offsets']
        table = rows(standard/'tsv/matrices.tsv')
        for model in ('old', 'svd', 'start', 'beam'):
            record, = [r for r in table if r['model'] == model]
            actual = [float(record[k]) for k in ('xptar_offset_mrad', 'ytar_offset_cm', 'yptar_offset_mrad')]
            expected = offsets if model == 'old' else [0., 0., 0.]
            if actual != expected:
                raise ValueError(f'{name}/{model}: incorrect external offsets {actual}')
            if sha(standard/f'matrices/{model}.dat') != record['sha256']:
                raise ValueError(f'{name}/{model}: diagnostic matrix checksum mismatch')
        if sha(beam/'matrices/beam.dat') != sha(standard/'matrices/beam.dat'):
            raise ValueError(f'{name}: standard diagnostics use a different final matrix')

        for path in (standard/'manifest.json', beam/'manifest.json'):
            add(path)
        for folder in (standard/'matrices', beam/'matrices'):
            for path in sorted(folder.glob('*.dat')):
                add(path)
        for filename in ('matrices.tsv', 'cell_comparisons.tsv', 'residuals.tsv',
                         'pooled_fold_scores.tsv', 'pooled_pooled.tsv', 'pooled_terms.tsv',
                         'beam_fold_scores.tsv', 'beam_conditioning.tsv',
                         'replay_closure.tsv', 'afterburner_iterations.tsv'):
            add(standard/'tsv'/filename)
        for filename in ('final_fit.tsv', 'terms.tsv', 'conditioning.tsv', 'fold_scores.tsv'):
            add(beam/'tsv'/filename)
        for filename in ('core_shoulder_change.png', 'core_shoulder_foil_delta.png',
                         'final_beam_core_shoulder_foil_delta.png',
                         'reduced_basis_folds.png', 'beam_folds.png'):
            add(standard/'plots'/filename)
        for path in sorted((campaign/'config').glob('rungroups_*_inputs.tsv')):
            add(path)
        for filename in ('sieve_mask.json', 'optics.dat'):
            add(campaign/'config'/filename)

        views = []
        replay_inputs = []
        angle_name = f'{angle:.3f}'.replace('.', 'p')
        for label, matrix, additions in (
            ('2024flags', standard/'matrices/old.dat', offsets),
            ('beam20', beam/'matrices/beam.dat', [0., 0., 0.]),
        ):
            view_name = f'sieve_{label}_theta{angle_name}_all_delta'
            view = campaign/'07_diagnostics'/view_name
            if not view.exists() and make_sieve:
                command = [sys.executable, str(PROJECT/'sieve_afterburner.py'), str(campaign), view_name,
                           '--matrix', str(matrix), '--offsets', *map(str, additions),
                           '--angle', str(angle), '--delta', '-10', '10',
                           '--cer', '6', '--cal', '0.65', '--foil-width', '2']
                subprocess.run(command, cwd=PROJECT, check=True)
            if not view.exists():
                raise FileNotFoundError(f'{view}: rerun this export with --make-sieve to generate the missing view')
            saved = json.loads((view/'manifest.json').read_text())
            settings = saved['settings']
            actual_offsets = settings['matrix_offsets']
            if (saved.get('status') != 'complete' or actual_offsets != additions
                    or abs(settings['angle']-angle) > 1e-6
                    or settings['delta'] != [-10., 10.]
                    or any(settings[k] != v for k, v in [('cer', 6.), ('cal', .65), ('foil_width', 2.)])):
                raise ValueError(f'{view}: selection or offset settings differ from the comparison')
            if sha(view/'matrices/matrix.dat') != sha(matrix):
                raise ValueError(f'{view}: sieve plot uses a different matrix')
            replay_inputs.append(saved['replays'])
            for filename in ('manifest.json', 'histograms.npz', 'tsv/counts.tsv', 'tsv/closure.tsv'):
                add(view/filename)
            views.append(str(view.relative_to(PROJECT)))
        old_counts, new_counts = [rows(PROJECT/v/'tsv/counts.tsv') for v in views]
        count_key = lambda r: (r['rungroup'], r['zfoil'])
        selected = lambda rs: {count_key(r): int(r['n']) for r in rs if r['model'] == 'matrix'}
        if selected(old_counts) != selected(new_counts):
            raise ValueError(f'{name}: before/after sieve event counts differ')
        if replay_inputs[0] != replay_inputs[1]:
            raise ValueError(f'{name}: before/after sieve replay inputs differ')
        campaigns.append(dict(campaign=name, tag=tag, representative_angle=angle,
                              standard=str(standard.relative_to(PROJECT)), sieve_views=views))
        print(f'{name}: offset policy, final matrix and matched sieve counts verified', flush=True)

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite {output}')
    metadata = dict(schema='nps_collaboration_export_v1', campaigns=campaigns,
                    historical_reference='2024 replay matrix + pinned historical flags offsets',
                    supported_population='Protected quality-2 cores plus all supported quality-1 shoulders',
                    files={key:dict(sha256=sha(path), bytes=path.stat().st_size) for key, path in files.items()})
    temporary = output.with_suffix('.partial')
    try:
        with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr('export_manifest.json', json.dumps(metadata, indent=2)+'\n')
            for key, path in sorted(files.items()):
                archive.write(path, key)
        temporary.rename(output)
    finally:
        temporary.unlink(missing_ok=True)
    print(f'Exported {len(files)} files: {output} ({output.stat().st_size/1024**2:.1f} MiB)')
    print('No ROOT event files or fits were exported or modified.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=PROJECT/'collaboration_exports/nps_angular_update_20261103.zip')
    parser.add_argument('--make-sieve', action='store_true', help='Generate missing representative sieve views from existing final matrices; never refit')
    args = parser.parse_args()
    try:
        export(args.output.resolve(), args.make_sieve)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f'ERROR: {exc}')
