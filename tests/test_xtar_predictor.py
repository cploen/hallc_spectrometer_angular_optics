"""Single-target units, direct-X equivalence, rank safety and separate export."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from beam_search import Folds, fit_all
from elastic_net import scaling, transform, macro_mse
from xtar_predictor import fit_training_design, export_model, predict


class SeparateXtarTests(unittest.TestCase):
    def test_single_target_qr_svd_matches_scaled_direct_x_in_mm(self):
        rng = np.random.default_rng(19)
        x = rng.normal(size=(180, 5)) * [0.2, 0.01, 0.02, 0.001, 2.] + [1, 2, 3, 4, 5]
        x[:, 4] = x[:, 0] + 1e-6*x[:, 4]
        y = (0.3 + x @ np.array([0.2, -0.1, 0.3, 0.5, 0.2]) + rng.normal(size=180)*0.007)[:, None]
        folds = np.arange(len(x)) % 3
        cells = np.array(['a']*110 + ['b']*50 + ['sparse']*20)
        for intercept in (True, False):
            prepared = Folds(x, y, folds, cells, 1e-12, 1, intercept, targets=('xtar',), scales=(10.,))
            for mask in ((), (0, 2), (0, 1, 2, 3, 4)):
                actual, _ = prepared.predict(0, mask)
                direct = np.empty(len(x))
                for f in range(3):
                    tr, va = folds != f, folds == f
                    state = scaling(x[tr], y[tr], intercept)
                    z, t = transform(x[tr], y[tr], state)
                    zv, _ = transform(x[va], y[va], state)
                    beta = np.linalg.lstsq(z[:, mask], t[:, 0], rcond=1e-12)[0]
                    direct[va] = zv[:, mask] @ beta*state[3][0] + state[2][0]
                np.testing.assert_allclose(actual, direct, atol=2e-8, rtol=2e-8)
                expected = float(macro_mse((direct-y[:, 0])*10., cells))
                self.assertAlmostEqual(prepared.evaluate(0, mask)['score'], expected, places=5)
            c, stats = fit_all(x, y, [(0, 2)], 1e-12, intercept, targets=('xtar',), scales=(10.,))
            self.assertEqual(c.shape, (6, 1))
            self.assertEqual(stats[0]['target'], 'xtar')
            if not intercept: self.assertEqual(c[0, 0], 0.)

    def test_response_shape_and_rank_rejected(self):
        x = np.tile(np.arange(30.)[:, None], (1, 2)); y = x[:, :1]
        with self.assertRaisesRegex(ValueError, 'shape'):
            Folds(x, y, np.arange(30)%3, np.zeros(30), 1e-12)
        with self.assertRaisesRegex(ValueError, 'scale'):
            Folds(x, y, np.arange(30)%3, np.zeros(30), 1e-12, targets=('xtar',), scales=(1000., 100.))
        with self.assertRaisesRegex(ValueError, 'rank deficient'):
            fit_all(x, y, [(0, 1)], 1e-12, targets=('xtar',), scales=(10.,))

    def test_training_selection_separate_export_and_prediction_identity(self):
        rng = np.random.default_rng(731)
        x = rng.normal(size=(150, 4))
        target = 0.12 + 0.4*x[:, 0] - 0.2*x[:, 2] + 0.003*rng.normal(size=150)
        folds = np.arange(150)%3; cells = np.array(['a']*100+['b']*50)
        with contextlib.redirect_stdout(io.StringIO()):
            fit = fit_training_design(x, target, folds, cells,
                en_config=dict(alphas=[.01, .001], l1=[.5], max_iter=2000),
                beam_config=dict(threads=1, steps=2))
        self.assertEqual(fit['beam_config']['beam'], 30)
        self.assertTrue(all(row['target']=='xtar' for row in fit['enet_cv']))
        features = [dict(name=f'x{i}', expression=f'x{i}', unit='m', source='synthetic measured fixture') for i in range(4)]
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'predictor.json'
            sentinel = Path(tmp)/'oldfit.dat'; sentinel.write_bytes(b'untouched transport\n')
            model = export_model(p, fit, features=features,
                                 provenance={'test':'synthetic'}, validity={'test':'synthetic only'})
            np.testing.assert_allclose(predict(model, x), np.c_[np.ones(len(x)), x] @ fit['coefficients_cm'])
            self.assertEqual(sentinel.read_bytes(), b'untouched transport\n')
            with self.assertRaises(FileExistsError):
                export_model(p, fit, features=features, provenance={'test':'synthetic'}, validity={'test':'synthetic only'})
            reloaded = json.loads(p.read_text())
            self.assertEqual(reloaded['output_unit'], 'cm')
            self.assertEqual(reloaded['reporting_scale'], 10.)
            # A future candidate/holdout distribution must not recenter the export.
            shifted = x+100
            np.testing.assert_allclose(predict(reloaded, shifted), fit['coefficients_cm'][0]+shifted@fit['coefficients_cm'][1:])


if __name__ == '__main__': unittest.main()
