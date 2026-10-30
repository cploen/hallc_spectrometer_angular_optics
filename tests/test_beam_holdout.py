"""Frozen mixed-sample holdout checks, exact large IDs, and rendered beam plots."""
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
import uproot

from test_elastic import fixture as core_fixture
from core_sample import digest, read_tsv, write_tsv
from fit_elastic import load_sample
from compare_matrices import read_matrix
from beam_holdout import prepare, predict, run


def write_root(path, tree, arrays):
    with uproot.recreate(path) as root:
        root.mktree(tree, arrays)


def fixture(root):
    campaign = core_fixture(root)
    allocation = campaign/'05c_core_sample/balanced'
    builds = campaign/'06c_core_ntuple/balanced'
    path = allocation/'root/CoreSample_rg01_test.root'
    with uproot.open(path) as f:
        mask = f['CoreSample'].arrays(library='np')
    for foil in range(5):
        idx = np.flatnonzero((mask['sample'] == 1) & (mask['foil'] == foil))
        mask['sample'][idx[90:]] = 0
        mask['quality'][idx[:18]] = 1
        mask['core_keep'][idx[:18]] = 0
    mask['core_score'] = np.where(mask['quality'] == 1, .45, np.where(mask['quality'] == 2, .7, np.nan))
    write_root(path, 'CoreSample', mask)
    selected = [dict(rungroup='rg01_test', entry=int(e)) for e in mask['entry'][mask['sample'] == 1]]
    write_tsv(allocation/'tsv/selected_ids.tsv', selected)
    m = json.loads((allocation/'manifest.json').read_text())
    m.update(schema='inner20_allocation_v1', foil_budget=90, totals=dict(fit=450, core=360, inner_shoulder=90),
             policy=dict(core_fraction=.8, inner_shoulder_fraction=.2, event_weight=1, inner_quality=1,
                         q_low_inclusive=.4, q_high_exclusive=.5))
    m['outputs'] = {str(p.relative_to(allocation)): digest(p) for p in allocation.rglob('*')
                    if p.is_file() and p.name != 'manifest.json'}
    (allocation/'manifest.json').write_text(json.dumps(m))
    for sample in ('fit', 'holdout'):
        build = builds/sample
        path = build/'root/Optics_999_-1_fit_tree_gmm.root'
        with uproot.open(path) as f:
            a = f['TFit'].arrays(library='np')
        if sample == 'fit':
            chosen = np.isin(a['entry'], mask['entry'][mask['sample'] == 1])
            a = {k: v[chosen] for k, v in a.items()}
        idx = np.searchsorted(mask['entry'], a['entry'])
        for field in ('quality', 'core_keep', 'core_score'):
            a[field] = mask[field][idx]
        write_root(path, 'TFit', a)
        if sample == 'fit':
            shutil.copy2(allocation/'manifest.json', build/'allocation_manifest.json')
            shutil.copy2(allocation/'tsv/selected_ids.tsv', build/'selected_ids.tsv')
        bm = dict(schema='mixed_preallocated_v1', sample=sample, sample_manifest=digest(allocation/'manifest.json'))
        bm['outputs'] = {str(p.relative_to(build)): digest(p) for p in build.rglob('*')
                         if p.is_file() and p.name != 'build.json'}
        (build/'build.json').write_text(json.dumps(bm))
    a, hashes = load_sample(campaign, 'balanced', 'fit')
    source = campaign/'06e_beam_search/beam20'
    (source/'matrices').mkdir(parents=True)
    (source/'tsv').mkdir()
    matrix = read_matrix(builds/'fit/oldfit.dat')
    for model, change in [('svd', 0.), ('start', .0001), ('beam', .0002)]:
        rows = [(c.copy(), t) for c, t in matrix]
        rows[0][0][0] += change
        (source/f'matrices/{model}.dat').write_text('\n'.join(' '.join(map(str, c))+' '+''.join(map(str,t)) for c,t in rows)+'\n')
    write_tsv(source/'tsv/folds.tsv', [dict(rungroup=r, entry=int(e), fold=i%3)
              for i,(r,e) in enumerate(zip(a['rungroup'], a['entry']))])
    write_tsv(source/'tsv/final_fit.tsv', [dict(target='ytar', model='beam', terms=4, rank=3)])
    bm = dict(schema='beam_v1', status='complete', sample='balanced', inputs=hashes, qa_min=10)
    bm['outputs'] = {str(p.relative_to(source)): digest(p) for p in source.rglob('*') if p.is_file()}
    (source/'manifest.json').write_text(json.dumps(bm))
    return campaign, source, bm


class BeamHoldoutTests(unittest.TestCase):
    def test_mixed_membership_predictions_and_real_plots(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign, source, _ = fixture(Path(temporary))
            _, _, a, matrices, pool, counts, _ = prepare(campaign, 'balanced', 'beam20')
            self.assertEqual(counts, dict(core=120, noncore=120, unsupported=120, blocked=90))
            self.assertGreater(int(a['entry'].min()), 2**54)
            predicted = predict(a, matrices, threads=1)
            for model in matrices:
                expected = np.zeros_like(predicted[model])
                for c,t in matrices[model]:
                    variables = [a['xfp']/100, a['xpfp'], a['yfp']/100, a['ypfp'], a['xtar']/100]
                    term = np.ones(len(pool))
                    for variable, power in zip(variables,t): term *= variable**power
                    expected += term[:,None]*c[:3]
                np.testing.assert_allclose(predicted[model], expected, rtol=1e-13, atol=1e-13)
            before = {p: digest(p) for p in campaign.rglob('*') if p.is_file()}
            run(campaign, 'balanced', 'holdout', 'beam20', threads=1, check=True)
            self.assertFalse((campaign/'07_diagnostics').exists())
            run(campaign, 'balanced', 'holdout', 'beam20', threads=1)
            out = campaign/'07_diagnostics/holdout'
            manifest = json.loads((out/'manifest.json').read_text())
            self.assertEqual(manifest['counts'], counts)
            for rel,sha in manifest['outputs'].items(): self.assertEqual(digest(out/rel),sha)
            ids = read_tsv(out/'tsv/ids.tsv')
            self.assertEqual({int(r['entry']) for r in ids}, set(map(int,a['entry'])))
            rows = read_tsv(out/'tsv/metrics.tsv')
            self.assertNotIn('blocked', {r['pool'] for r in rows})
            self.assertEqual({r['model'] for r in rows}, {'svd','start','beam'})
            self.assertEqual(len(list((out/'plots').glob('*.png'))), 8)
            self.assertEqual(before, {p: digest(p) for p in before})
            with self.assertRaises(FileExistsError): run(campaign, 'balanced', 'holdout', 'beam20')
            if os.environ.get('BEAM_HOLDOUT_PREVIEW'):
                target = Path(os.environ['BEAM_HOLDOUT_PREVIEW'])
                target.mkdir(parents=True, exist_ok=True)
                for name in ('residuals.png','foil_delta.png','slice_0_core_clouds.png'):
                    shutil.copy2(out/'plots'/name, target/name)

    def test_incomplete_changed_and_duplicate_training_inputs_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign, source, manifest = fixture(Path(temporary))
            manifest['status'] = 'incomplete'
            (source/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'completed beam'):
                prepare(campaign,'balanced','beam20')
            manifest['status'] = 'complete'
            (source/'manifest.json').write_text(json.dumps(manifest))
            original = (source/'matrices/beam.dat').read_bytes()
            (source/'matrices/beam.dat').write_bytes(original+b'\n')
            with self.assertRaisesRegex(ValueError, 'Frozen input changed'):
                prepare(campaign,'balanced','beam20')
            (source/'matrices/beam.dat').write_bytes(original)
            p = source/'tsv/folds.tsv'
            rows = read_tsv(p); rows.append(rows[0])
            write_tsv(p, rows)
            manifest['outputs']['tsv/folds.tsv'] = digest(p)
            (source/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'training IDs'):
                prepare(campaign,'balanced','beam20')


if __name__ == '__main__':
    unittest.main()
