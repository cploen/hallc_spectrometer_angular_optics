#!/usr/bin/env python3
"""Evaluate frozen beam/seed/full-basis matrices on exact reserved holdout IDs."""
import argparse
import json
from pathlib import Path
import re
import shutil
import tempfile

import numpy as np
from threadpoolctl import threadpool_limits

from core_sample import PROJECT, digest, read_tsv, write_tsv
from compare_matrices import read_matrix
from fit_elastic import checked, design, load_sample
from spectrometer_config import from_campaign

MODELS = ('svd', 'start', 'beam')


def prepare(campaign, tag, source):
    saved = campaign/'06e_beam_search'/source
    manifest = json.loads((saved/'manifest.json').read_text())
    if (manifest.get('schema'), manifest.get('status'), manifest.get('sample')) != ('beam_v1', 'complete', tag):
        raise ValueError('Requires a completed beam result for the requested allocation tag')
    hashes = {str(saved/'manifest.json'): digest(saved/'manifest.json')}
    for path, sha in manifest['inputs'].items():
        checked(Path(path), sha, hashes)
    required = ['tsv/folds.tsv', 'tsv/final_fit.tsv'] + [f'matrices/{m}.dat' for m in MODELS]
    for rel in required:
        checked(saved/rel, manifest['outputs'][rel], hashes)
    allocated = campaign/'05c_core_sample'/tag
    suffix = f'/05c_core_sample/{tag}/manifest.json'
    original = {sha for path, sha in manifest['inputs'].items() if path.endswith(suffix)}
    if original != {digest(allocated/'manifest.json')}:
        raise ValueError('Allocation differs from the beam training allocation')
    rows = read_tsv(saved/'tsv/folds.tsv')
    training = {(r['rungroup'], int(r['entry'])) for r in rows}
    intended = read_tsv(allocated/'tsv/selected_ids.tsv')
    expected = {(r['rungroup'], int(r['entry'])) for r in intended}
    if len(training) != len(rows) or len(expected) != len(intended) or training != expected:
        raise ValueError('Beam training IDs differ from allocation IDs')
    a, inputs = load_sample(campaign, tag, 'holdout')
    hashes.update(inputs)
    identities = list(zip(a['rungroup'], map(int, a['entry'])))
    if len(set(identities)) != len(identities) or training.intersection(identities):
        raise ValueError('Duplicate holdout IDs or training/holdout overlap')
    matrices = {m: read_matrix(saved/f'matrices/{m}.dat') for m in MODELS}
    # Fixed xtar/delta terms must agree across the three frozen angular candidates.
    fixed = [{t: c for c, t in matrices[m] if t[4]} for m in MODELS]
    if any(set(f) != set(fixed[0]) or any(not np.array_equal(f[t], fixed[0][t]) for t in f)
           for f in fixed[1:]):
        raise ValueError('Fixed xtar terms differ between beam reference matrices')
    delta = [{t: c[3] for c, t in matrices[m] if c[3] != 0} for m in MODELS]
    if any(v != delta[0] for v in delta[1:]):
        raise ValueError('Delta coefficients differ between beam reference matrices')
    pool = np.where(a['excluded'], 'blocked', np.where(a['quality'] == 2, 'core',
                    np.where(a['quality'] == 1, 'noncore', 'unsupported')))
    counts = {p: int(np.sum(pool == p)) for p in ('core', 'noncore', 'unsupported', 'blocked')}
    return saved, manifest, a, matrices, pool, counts, hashes


def predict(a, matrices, threads=8):
    terms = sorted({t for rows in matrices.values() for c, t in rows})
    lookup = {t: i for i, t in enumerate(terms)}
    coeff = np.zeros((len(terms), 3*len(MODELS)))
    for j, m in enumerate(MODELS):
        for c, t in matrices[m]:
            coeff[lookup[t], 3*j:3*j+3] = c[:3]
    values = np.empty((len(a['entry']), 3*len(MODELS)))
    with threadpool_limits(limits=threads):
        for first in range(0, len(values), 20000):
            piece = {k: v[first:first+20000] for k, v in a.items()}
            values[first:first+20000] = design(piece, terms) @ coeff
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite frozen-matrix prediction; no events silently discarded')
    return {m: values[:, 3*j:3*j+3] for j, m in enumerate(MODELS)}


def run(campaign, tag, name, source, threads=8, check=False):
    if threads < 1:
        raise ValueError('threads must be positive')
    target = campaign/'07_diagnostics'/name
    if not check and target.exists():
        raise FileExistsError(f'Output exists: {target}; choose a new name')
    saved, beam, a, matrices, pool, counts, hashes = prepare(campaign, tag, source)
    print(f'Protected holdout: {len(pool):,} exact IDs; pools: {counts}', flush=True)
    if check:
        print('Frozen beam/seed/SVD matrices and disjoint membership verified. No residuals evaluated; no outputs written.')
        return
    from elastic_diagnostics import evaluation_plots, metrics, residuals
    target.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.beam_holdout_', dir=target.parent))
    try:
        for folder in ('plots', 'tsv', 'matrices', 'code'):
            (stage/folder).mkdir()
        predictions = predict(a, matrices, threads)
        truth = np.column_stack((a['xptarT'], a['ytarT']/100., a['yptarT']))
        offset = np.zeros_like(truth)  # Complete transport matrices include constants and fixed terms.
        qa_min = beam.get('qa_min', 10)
        write_tsv(stage/'tsv/metrics.tsv', metrics(a, truth, predictions, pool, qa_min, offset))
        write_tsv(stage/'tsv/ids.tsv', [dict(rungroup=r, entry=int(e), pool=p)
                   for r, e, p in zip(a['rungroup'], a['entry'], pool)])
        rr = {m: residuals(a, truth, p, offset) for m, p in predictions.items()}
        np.savez_compressed(stage/'residuals.npz', rungroup=a['rungroup'], entry=a['entry'], pool=pool,
                            targets=np.array(['xptar', 'ytar', 'yptar', 'ztar']), **rr)
        evaluation_plots(stage, a, truth, predictions, offset, pool, from_campaign(campaign),
                         qa_min, candidate='beam', candidate_label='Beam')
        for m in MODELS:
            shutil.copy2(saved/f'matrices/{m}.dat', stage/f'matrices/{m}.dat')
        shutil.copy2(saved/'tsv/final_fit.tsv', stage/'tsv/training_conditioning.tsv')
        shutil.copy2(campaign/'05c_core_sample'/tag/'metadata/optics.dat', stage/'optics.dat')
        for file in ('beam_holdout.py', 'run_beam_holdout.sh', 'elastic_diagnostics.py',
                     'fit_elastic.py', 'elastic_net.py', 'preallocated_svd.py', 'mixed_sample.py',
                     'core_sample.py', 'spectrometer_config.py', 'spectrometer_profiles.def', 'compare_matrices.py'):
            shutil.copy2(PROJECT/file, stage/'code'/file)
        (stage/'README.md').write_text(
            '# Frozen beam protected-holdout diagnostics\n\n'
            f'Exact reserved holdout IDs: {len(pool):,}. Pool counts: {counts}.\n\n'
            'svd = full-basis SVD; start = selected starting basis; beam = chosen beam basis. '
            'All three are final full-training matrices. No coefficients are fitted here. '
            'Ratio maps compare Beam / SVD; residual distributions include all three. '
            'Sieve clouds compare saved replay, SVD and beam on identical events. '
            'Core, supported noncore and unsupported pools remain separate. Blocked labels have counts only '
            'in the summary tables and figures. Sparse-cell flags do not remove events.\n\n'
            'These diagnostics use the saved replay xtar, vertex and momentum inputs. The ztar residual is '
            'derived from the HMS target identity. Legacy afterburner PDFs instead iterate xtar with the '
            'saved reaction vertex; their auxiliary ztar/delta branches retain replay values. Neither is a new HCANA replay.\n\n'
            'Read plots/residuals.png, foil_delta.png, slice_*_holes.png and slice_*_clouds.png. '
            'tsv/metrics.tsv records bias, RMS, P95, counts and sparse-cell flags down to setting/hole. '
            'tsv/ids.tsv and residuals.npz retain exact identities.\n')
        report = dict(schema='beam_holdout_v1', status='complete', campaign=campaign.name, sample=tag,
                      source=source, counts=counts, inputs=hashes, threads=threads,
                      reconstruction='Complete frozen matrices; saved replay xtar; no external offsets',
                      qa_min=qa_min)
        # Refuse to certify inputs changed during diagnostics.
        for path, sha in hashes.items():
            checked(Path(path), sha, {})
        report['outputs'] = {str(p.relative_to(stage)): digest(p) for p in stage.rglob('*') if p.is_file()}
        (stage/'manifest.json').write_text(json.dumps(report, indent=2)+'\n')
        stage.rename(target)
    except BaseException:
        shutil.rmtree(stage)
        raise
    print(f'Protected beam diagnostics complete: {target}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign', type=Path)
    p.add_argument('tag')
    p.add_argument('name')
    p.add_argument('--source', required=True)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    if not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', v) for v in (a.tag, a.name, a.source)):
        p.error('Invalid tag/name/source')
    campaign = (a.campaign if a.campaign.is_absolute() else PROJECT/a.campaign).resolve()
    try:
        run(campaign, a.tag, a.name, a.source, a.threads, a.check)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        p.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
