"""Numerical controls for the portable continuous-delta engine (no ROOT)."""
import json
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from joint_surface import fit_axis, predict_axis, surface_derivative, surface_values
from joint_surface.engine import _coverage_records, _domain, _occupancy, _configuration, _coordinates


def sample(n, low=-10, high=10, count=7, seed=42):
    rng = np.random.RandomState(seed)
    delta = rng.uniform(low, high, n); j = rng.randint(count, size=n)
    r = rng.uniform(-1.4, 1.4, n)+.2*(delta-(low+high)/2)/((high-low)/2)
    D = (delta-(low+high)/2)/((high-low)/2)
    t = (j-(count-1)/2)/((count-1)/2)
    q = .40*(j-(count-1)/2)+.18*r+.055*r*r+.055*D+.025*t*D+rng.normal(0, .012, n)
    theta = np.deg2rad(135.)
    z = np.c_[q*np.cos(theta)-r*np.sin(theta), q*np.sin(theta)+r*np.cos(theta)]
    frame = dict(mean=[0., 0.], scale=[.02, 30.], theta_deg=135.)
    return z*frame['scale'], delta, j, r, frame


class JointSurfaceEngineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw, delta, cls.truth, cls.r, cls.frame = sample(16000)
        cls.raw, cls.delta = raw, delta
        cls.cfg = dict(phase_starts=4, phase_iterations=60, iterations=3, irls=4,
                       max_fit_events=7000, families=['bilinear', 'curved', 'shift', 'spacing'], seed=86237)
        cls.model, cls.meta = fit_axis(raw[:13000], delta[:13000], [-10, -8, -5, 0, 5, 10], cls.cfg, cls.frame)

    def test_discovery_and_independent_prediction(self):
        self.assertEqual(self.model['count'], 7)
        self.assertFalse(self.meta['reconstructed_sieve_inputs_used'])
        self.assertFalse(self.meta['reference_labels_used'])
        self.assertEqual(set(self.meta['training_indices']) & set(self.meta['model_choice_indices']), set())
        p = predict_axis(self.model, self.raw[13000:], self.delta[13000:])
        self.assertGreater(p['accepted'].mean(), .95)
        self.assertGreater(np.mean(p['labels'][p['accepted']] == self.truth[13000:][p['accepted']]), .995)
        self.assertGreater(min(x['minimum_order'] for x in self.model['diagnostics']['history']), 0)
        restored = json.loads(json.dumps(self.model, allow_nan=False))
        np.testing.assert_array_equal(predict_axis(restored, self.raw[13000:], self.delta[13000:])['labels'], p['labels'])

    def test_count_changes_only_coverage(self):
        before = json.dumps(self.model, sort_keys=True)
        a = predict_axis(self.model, self.raw[13000:], self.delta[13000:], coverage_min_events=10)
        b = predict_axis(self.model, self.raw[13000:], self.delta[13000:], coverage_min_events=1000)
        self.assertGreater(a['accepted'].sum(), b['accepted'].sum())
        for key in ['q', 'r', 'u', 'finite', 'occupancy', 'posterior', 'proposed_label']:
            np.testing.assert_array_equal(a[key], b[key])
        self.assertEqual(json.dumps(self.model, sort_keys=True), before)

    def test_missing_population_not_detected(self):
        # Geometry was learned with seven ridges. Remove one ridge from one
        # target reporting interval; it must not become an occupied detection.
        mask = ~((self.truth == 3) & (self.delta >= -10) & (self.delta < -8))
        cfg = _configuration(self.cfg)
        occ = _occupancy(self.model, self.raw[mask], self.delta[mask], cfg)
        self.assertNotIn(3, occ['0']['accepted_indices'])
        self.assertIn(3, occ['2']['accepted_indices'])
        self.assertEqual(self.model['count'], 7)

    def test_exact_delta_and_nonfinite_reasons(self):
        raw = np.repeat(self.raw[0:1], 4, axis=0); raw[3, 0] = np.nan
        p = predict_axis(self.model, raw, np.array([-10., 10., 22., 0.]))
        self.assertNotEqual(p['reason'][0], 'delta_range')
        self.assertEqual(p['reason'][1], 'delta_range')
        self.assertEqual(p['reason'][2], 'delta_range')
        self.assertEqual(p['reason'][3], 'nonfinite')

    def test_order_over_observed_domain(self):
        r = np.linspace(-1.3, 1.3, 31); delta = np.linspace(-9.9, 9.9, 41)
        rr, dd = np.meshgrid(r, delta)
        for j in np.linspace(-.5, self.model['count']-.5, 9):
            self.assertGreater(surface_derivative(self.model, j, rr, dd).min(), 0)
        self.assertTrue(np.all(surface_values(self.model, 1, rr, dd) < surface_values(self.model, 2, rr, dd)))

    def test_shms_range_and_eleven_is_maximum(self):
        raw, delta, truth, r, frame = sample(14000, low=-10, high=22, count=9, seed=64)
        model, meta = fit_axis(raw[:12000], delta[:12000], [-10, -8, -5, 0, 5, 10, 15, 22], self.cfg, frame)
        self.assertEqual(model['dcenter'], 6.)
        self.assertEqual(model['dscale'], 16.)
        self.assertEqual(model['count'], 9)  # Not forced to SHMS hardware maximum11.
        pred = predict_axis(model, raw[12000:], delta[12000:])
        upper = delta[12000:] >= 10
        self.assertGreater(pred['accepted'][upper].mean(), .9)
        self.assertGreater(np.mean(pred['labels'][pred['accepted']] == truth[12000:][pred['accepted']]), .99)
        self.assertTrue(any(row[0] > 10 for rows in model['coverage'].values() for row in rows))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            fit_axis(self.raw[:600], self.delta[:600], [-10, 0, 10], dict(unknown_key=True))
        with self.assertRaises(ValueError):
            fit_axis(self.raw[:600], self.delta[:600], [-10, -10, 10])
        with self.assertRaises(ValueError):
            fit_axis(self.raw[:600], self.delta[:600], [-5, 0, 10])


if __name__ == '__main__':
    unittest.main()
