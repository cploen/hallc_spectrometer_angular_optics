"""Target units, integer IDs, stale-truth rejection and real legacy ROOT output."""
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot

from core_sample import digest, write_tsv
from legacy_afterburner import REQUIRED, UPDATED, export_tree, prepare, run
from spectrometer_config import hms_xptar_truth

ANGLE = 12.49
XMIS = .1 * (2.37 - .086 * ANGLE + .0012 * ANGLE ** 2)


class LegacyAfterburnerTest(unittest.TestCase):
    def fixture(self, directory):
        campaign = Path(directory) / 'HMS_legacy_fixture'
        config = campaign / 'config'
        config.mkdir(parents=True)
        fit_dir = campaign / '06c_core_ntuple/equal15/fit/root'
        fit_dir.mkdir(parents=True)
        # This synthetic ID exists only in the campaign snapshot.
        (config / 'optics.dat').write_text('999999,999999,12.49,1,1,2,0\n0\n-10,10\n')
        write_tsv(config / 'rungroups_fixture_inputs.tsv', [dict(
            rungroup='fixture', optics_id=999999, hms_angle_deg=ANGLE)])
        n = 24
        arrays = {k: np.zeros(n, np.float64) for k in REQUIRED if k != 'entry'}
        arrays.update(entry=np.arange(n, dtype=np.int64) + 2 ** 54 + 1,
                      reactycalc=np.full(n, -XMIS), foil=np.zeros(n, np.int32),
                      xscol=np.full(n, 4, np.int32), yscol=np.full(n, 4, np.int32),
                      sample=np.ones(n, np.int32), core_keep=np.ones(n, np.int32))
        arrays['xptarT'] = hms_xptar_truth(arrays['xsT'], arrays['ztarT'], ANGLE,
                                        arrays['reactycalc'])
        source = fit_dir / 'Optics_999999_-1_fit_tree_gmm.root'
        self.write_tree(source, arrays)
        matrix = campaign / 'new.dat'
        matrix.write_text('.003 .00123 .002 0 00000\n')
        return campaign, fit_dir, source, matrix, arrays

    @staticmethod
    def write_tree(path, arrays):
        with uproot.recreate(path) as root:
            root.mktree('TFit', arrays)

    def test_export_units_constants_xtar_offsets_and_membership(self):
        with tempfile.TemporaryDirectory() as tmp:
            campaign, _, source, _, original = self.fixture(tmp)
            original['reactycalc'] = np.linspace(-.2, .2, len(original['entry']))
            original['ztar'] = np.linspace(-8, 8, len(original['entry']))
            original['ztarT'] = original['ztar'].copy()
            original['xptarT'] = hms_xptar_truth(original['xsT'], original['ztarT'],
                                               ANGLE, original['reactycalc'])
            self.write_tree(source, original)
            const = np.array([.003, .00123, .002, 0.])
            xtar_coeff = np.array([.01, .02, .03, 0.])
            for label, slope, offsets in [('constant', np.zeros(4), [0., 0., 0.]),
                                          ('xtar_offsets', xtar_coeff, [2., .3, -1.])]:
                with self.subTest(label=label):
                    rows = [(const, (0, 0, 0, 0, 0)), (slope, (0, 0, 0, 0, 1))]
                    target = campaign / (label + '.root')
                    report = export_tree(source, target, rows, ANGLE, offsets)
                    with uproot.open(target) as root:
                        exported = root['TFit'].arrays(library='np')
                        self.assertEqual(root['TFit'].classname, 'TTree')
                    self.assertEqual(report, dict(entries=len(original['entry']), iteration_limit=0))
                    # Both fixtures converge after their required first update.
                    x0 = -original['reactycalc'] - XMIS
                    additions = np.asarray(offsets) / [1000., 100., 1000.]
                    initial_xp = const[0] + slope[0] * x0 / 100 + additions[0]
                    xtar = x0 - initial_xp * original['ztar'] * np.cos(np.deg2rad(ANGLE))
                    values = const[:3] + slope[:3] * xtar[:, None] / 100 + additions
                    expected = dict(xtar=xtar, xptar=values[:, 0], ytar=100 * values[:, 1],
                                    yptar=values[:, 2], xs=xtar + 168 * values[:, 0],
                                    ys=100 * values[:, 1] + 168 * values[:, 2])
                    for key in UPDATED:
                        self.assertEqual(exported[key].dtype, np.dtype('float64'))
                        np.testing.assert_allclose(exported[key], expected[key], rtol=1e-13, atol=1e-13)
                        np.testing.assert_array_equal(exported['replay_' + key], original[key])
                    for key in set(original) - set(UPDATED):
                        np.testing.assert_array_equal(exported[key], original[key])
                        self.assertEqual(exported[key].dtype, original[key].dtype)
                    np.testing.assert_array_equal(exported['afterburner_iterations'], 1)

    def test_missing_branch_and_stale_truth_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            campaign, fit_dir, source, _, arrays = self.fixture(tmp)
            missing = dict(arrays)
            del missing['reactycalc']
            self.write_tree(source, missing)
            with self.assertRaisesRegex(ValueError, 'missing TFit branches.*reactycalc'):
                prepare(campaign, fit_dir)
            arrays['xptarT'][0] += .001
            self.write_tree(source, arrays)
            with self.assertRaisesRegex(ValueError, 'Stale or inconsistent HMS xptarT'):
                export_tree(source, campaign / 'stale.root',
                            [(np.zeros(4), (0, 0, 0, 0, 0))], ANGLE, [0, 0, 0])

    @unittest.skipUnless(shutil.which('root'), 'ROOT required for legacy integration')
    def test_real_root_products_manifest_dry_run_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            campaign, fit_dir, source, matrix, original = self.fixture(tmp)
            source_hash = digest(source)
            run(campaign, 'dry_run', matrix, fit_dir, dry_run=True)
            self.assertFalse((campaign / '07_diagnostics').exists())
            run(campaign, 'legacy', matrix, fit_dir)
            out = campaign / '07_diagnostics/legacy'
            manifest = json.loads((out / 'manifest.json').read_text())
            self.assertEqual(manifest['status'], 'complete')
            self.assertEqual(manifest['external_offsets'], [0., 0., 0.])
            self.assertEqual(manifest['groups'][0]['entries'], len(original['entry']))
            self.assertEqual(len(list(out.rglob('*.pdf'))), 7)
            self.assertEqual(len(list((out / 'residuals').rglob('*.root'))), 3)
            folder = out / 'residuals/fixture'
            checks = [('xptar', 'hXpDiff_0_DelCut_0_Xs_4', 3., 'fTsumwx'),
                      ('yptar', 'hYpDiff_0_DelCut_0_Ys_4', 2., 'fTsumwx'),
                      ('ytar', 'hYDiff_0_DelCut_0', .123, 'fTsumwx'),
                      ('ytar', 'hYDiff_0_DelCut_0_Xs_4', .123, 'fTsumwx'),
                      ('ytar', 'hYDiff_YTrue_0_DelCut_0', .123, 'fTsumwy')]
            for variable, name, expected, sum_member in checks:
                with uproot.open(folder / f'fixture_{variable}_residuals.root') as root:
                    hist = root[name]
                    self.assertEqual(hist.member('fEntries'), len(original['entry']))
                    self.assertAlmostEqual(hist.member(sum_member) / hist.member('fTsumw'), expected)
            self.assertEqual(digest(source), source_hash)
            with self.assertRaises(FileExistsError):
                run(campaign, 'legacy', matrix, fit_dir)


if __name__ == '__main__':
    unittest.main()
