"""Real ROOT handoff tests for frozen 90/10 and 80/20 training samples."""
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import core_sample as cs
from test_core_sample import fixture
from allocate_inner20 import run as allocate
from build_core_fit import build
from fit_elastic import load_sample
from mixed_sample import mixed_fraction, validate_training
from preallocated_svd import verify_build, run as standalone_svd
from reallocate_core import run as allocate_core


class AdmissionTests(unittest.TestCase):
    def test_boundaries_and_original_labels(self):
        a = dict(quality=np.array([2, 1, 1]), core_keep=np.array([1, 0, 0]),
                 core_score=np.array([.7, .4, np.nextafter(.5, 0)]),
                 balance_excluded=np.zeros(3, dtype=int))
        selected = np.ones(3, dtype=bool)
        validate_training(a, selected, True)
        with self.assertRaises(ValueError): validate_training(a, selected, False)
        for key, value in [('core_score', .5), ('core_score', .399), ('core_score', np.nan),
                           ('quality', 0), ('core_keep', 1), ('balance_excluded', 1)]:
            bad = {k: v.copy() for k, v in a.items()}
            bad[key][1] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                validate_training(bad, selected, True)
        with self.assertRaises(ValueError): mixed_fraction(dict(schema='unknown'))
        with self.assertRaises(ValueError): mixed_fraction(dict(schema='inner10_allocation_v1', policy={}))


@unittest.skipUnless(shutil.which('root'), 'ROOT executable required')
class RootExportTests(unittest.TestCase):
    def test_exact_membership_labels_holdout_and_elastic_loader(self):
        with tempfile.TemporaryDirectory(prefix='mixed_export_') as temporary:
            campaign = Path(temporary)/'HMS_fixture'
            fixture(campaign)
            metadata = campaign/'config/optics.dat'
            shutil.copy2(cs.PROJECT/'DATfiles/list_of_optics_run.dat', metadata)
            # A changed saved angle proves export uses frozen metadata rather than the global DB.
            text = metadata.read_text()
            text = text.replace('666701,666701,12.490,', '666701,666701,13.000,')
            metadata.write_text(text)
            (campaign/'config/core_sample.json').write_text(json.dumps(dict(min_events=10)))
            (campaign/'config/sieve_mask.json').write_text(json.dumps(dict(blocked=[[4, 4]])))
            (campaign/'config/oldfit.dat').write_text('! fixture\n ----\n 0 0 0 0 00000\n ----\n')
            version = importlib.metadata.version
            with patch('core_sample.plot_slice'), patch('core_sample.importlib.metadata.version',
                    side_effect=lambda n: 'fixture-no-plots' if n == 'matplotlib' else version(n)):
                cs.run(campaign, 'min10')
            parent = campaign/'05c_core_sample/min10'
            frozen = {p: cs.digest(p) for p in parent.rglob('*') if p.is_file()}
            with uproot.open(parent/'root/CoreSample_rg01_test.root') as root:
                original = root['CoreSample'].arrays(library='np')
            for percent in (10, 20):
                tag = f'inner{percent}'
                allocate(campaign, 'min10', tag, shoulder_percent=percent)
                source = campaign/'05c_core_sample'/tag
                manifest = json.loads((source/'manifest.json').read_text())
                with patch.dict(os.environ, HCANA='root'):
                    for sample in ('fit', 'holdout', 'surplus'):
                        build(campaign, tag, sample)
                        data, _ = load_sample(campaign, tag, sample)
                        code = dict(fit=1, holdout=2, surplus=3)[sample]
                        with uproot.open(source/'root/CoreSample_rg01_test.root') as root:
                            mask = root['CoreSample'].arrays(library='np')
                        expected = mask['sample'] == code
                        np.testing.assert_array_equal(data['entry'], mask['entry'][expected])
                        np.testing.assert_array_equal(data['quality'], mask['quality'][expected])
                        if sample == 'fit':
                            self.assertEqual(len(data['entry']), manifest['totals']['fit'])
                            self.assertEqual(int(np.sum(data['quality'] == 1))*100, percent*len(data['entry']))
                        if sample == 'holdout':
                            np.testing.assert_array_equal(data['entry'], original['entry'][original['sample'] == 2])
                    bm, ids = verify_build(campaign/'06c_core_ntuple'/tag/'fit')
                    self.assertEqual(bm['schema'], 'mixed_preallocated_v1')
                    self.assertEqual(len(ids), manifest['totals']['fit'])
                    with self.assertRaisesRegex(ValueError, 'standalone core-only'):
                        standalone_svd(campaign/'06c_core_ntuple'/tag/'fit', campaign/'unused_solver')
                # Bad intended IDs fail before ROOT is called, even with an updated TSV checksum.
                ids_path = source/'tsv/selected_ids.tsv'
                ids = cs.read_tsv(ids_path)
                ids[0]['entry'] = '-1'
                cs.write_tsv(ids_path, ids)
                manifest['outputs']['tsv/selected_ids.tsv'] = cs.digest(ids_path)
                (source/'manifest.json').write_text(json.dumps(manifest))
                shutil.rmtree(campaign/'06c_core_ntuple'/tag/'fit')
                with patch.dict(os.environ, HCANA='root'), self.assertRaisesRegex(ValueError, 'selected IDs disagree'):
                    build(campaign, tag, 'fit')
            # Existing equal15 core-only allocation/export/loader remains valid.
            with patch.dict(sys.modules, balance_diagnostics=SimpleNamespace(diagnostics=lambda *a: None)), \
                    patch('reallocate_core.importlib.metadata.version', side_effect=lambda n: 'fixture-no-plots' if n == 'matplotlib' else version(n)):
                allocate_core(campaign, 'min10', 'equal15')
            with patch.dict(os.environ, HCANA='root'):
                build(campaign, 'equal15', 'fit')
                data, _ = load_sample(campaign, 'equal15', 'fit')
            self.assertTrue(np.all(data['quality'] == 2))
            bm, _ = verify_build(campaign/'06c_core_ntuple/equal15/fit')
            self.assertEqual(bm['schema'], 'core_preallocated_v1')
            self.assertEqual(frozen, {p: cs.digest(p) for p in parent.rglob('*') if p.is_file()})


if __name__ == '__main__':
    unittest.main()
