"""Targeted numerical/cut-contract checks; these are software fixtures, not physics validation.

Run with a Python containing NumPy, SciPy, and Matplotlib. The companion
``geometric_root_smoke.py`` exercises ROOT and the downstream campaign stages.
"""
import sys
import json
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np
from matplotlib.path import Path as PolygonPath

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from geometric_selection import (fit_x, detect_x, predict_x, fit_y, predict_y,
                                 global_centers)
from geometric_selection.cuts import build_x_polygons, build_y_polygons


def x_fixture(slots=9, per_band=550, seed=411):
    """An ordered family with both fan divergence and weak curvature."""
    rng = np.random.RandomState(seed)
    geometry = dict(mean=[.012, 3.], scale=[.035, 28.], theta_deg=135.,
                    coefficients=[0., .11, .03, 1.7, .102, .034, 0., 0., 0.],
                    p_scale=2., j_degree=2, p_degree=2,
                    possible_indices=list(range(slots)), max_slots=slots,
                    diagnostics=dict(full_event_p_extent=[-2., 2.],
                                     p_support=[-2., 2.]))
    labels = np.repeat(np.arange(slots), per_band)
    p = rng.uniform(-1.9, 1.9, len(labels))
    curves = global_centers(geometry, p)
    q = curves[labels, np.arange(len(p))] + rng.normal(0., .014, len(p))
    t = np.deg2rad(geometry['theta_deg'])
    z = np.column_stack((np.cos(t)*q - np.sin(t)*p,
                         np.sin(t)*q + np.cos(t)*p))
    return z*np.asarray(geometry['scale'])+geometry['mean'], labels, geometry


def y_fixture(per_band=650, seed=702):
    """Five separated anisotropic clouds with unequal widths and small tilts."""
    rng = np.random.RandomState(seed)
    labels = np.repeat(np.arange(5), per_band)
    q = np.concatenate([rng.normal(i-2., .075+.008*i, per_band)
                        for i in range(5)])
    p = rng.normal(0., .55, len(labels)) + .035*(labels-2.)*q
    t = np.deg2rad(135.)
    raw = np.column_stack((np.cos(t)*q-np.sin(t)*p,
                           np.sin(t)*q+np.cos(t)*p))
    return raw*np.array([.015, 14.])+np.array([.002, -1.]), labels


def polygon_labels(polygons, raw):
    """Evaluate independent point-in-polygon membership, including components."""
    membership = {}
    for polygon in polygons:
        label = int(polygon['index'])
        inside = PolygonPath(np.asarray(polygon['points'])).contains_points(raw)
        membership[label] = membership.get(label, np.zeros(len(raw), bool)) | inside
    counts = np.zeros(len(raw), int)
    labels = np.full(len(raw), -1, int)
    for label, inside in membership.items():
        counts += inside
        labels[inside] = label
    return labels, counts


def agreement_up_to_permutation(predicted, truth):
    """Slot orientation/offset do not establish physical sieve labels."""
    return sum(np.bincount(truth[predicted == label]).max()
               for label in np.unique(predicted) if label >= 0)/len(truth)


class GeometricSelectionTests(unittest.TestCase):
    def test_curved_x_discovery(self):
        raw, truth, _ = x_fixture()
        geometry = fit_x(raw, ndel=2, angle=135., max_slots=9)
        occupied = detect_x(geometry, raw)
        prediction = predict_x(geometry, occupied, raw)
        self.assertEqual(len(occupied['accepted_indices']), 9)
        self.assertTrue(geometry['diagnostics']['no_split_anywhere'])
        self.assertTrue(geometry['theta_was_overridden'])
        self.assertGreater(agreement_up_to_permutation(prediction['labels'], truth), .99)

    def test_x_root_polygon_contract_and_missing_slot(self):
        raw, truth, geometry = x_fixture()
        occupied = detect_x(geometry, raw)
        self.assertEqual(occupied['accepted_indices'], list(range(9)))
        numeric = predict_x(geometry, occupied, raw)['labels']
        labels, counts = polygon_labels(build_x_polygons(geometry, occupied, raw), raw)
        self.assertEqual(int(counts.max()), 1, 'Cuts mix distinct provisional bands')
        self.assertGreater(np.mean(labels == numeric), .998)
        # An expected but absent interior slot must not be invented by the prior.
        keep = truth != 3
        sparse = detect_x(geometry, raw[keep])
        self.assertEqual(sparse['accepted_indices'], [0, 1, 2, 4, 5, 6, 7, 8])
        self.assertNotIn(3, predict_x(geometry, sparse, raw[keep])['labels'])

    def test_shms_eleven_slot_geometry(self):
        raw, truth, geometry = x_fixture(slots=11)
        occupied = detect_x(geometry, raw)
        self.assertEqual(occupied['accepted_indices'], list(range(11)))
        labels, counts = polygon_labels(build_x_polygons(geometry, occupied, raw), raw)
        self.assertLessEqual(int(counts.max()), 1)
        self.assertGreater(np.mean(labels == truth), .998)

    def test_manual_x_divider_changes_the_cut(self):
        raw, _, geometry = x_fixture()
        occupied = detect_x(geometry, raw)
        # Place a probe inside the first default cell, beyond a manual divider.
        q = float(global_centers(geometry, [0.], indices=[.35])[0, 0])
        t = np.deg2rad(geometry['theta_deg'])
        probe = np.array([[np.cos(t)*q, np.sin(t)*q]])*geometry['scale']+geometry['mean']
        default, _ = polygon_labels(build_x_polygons(geometry, occupied, raw), probe)
        dividers = np.arange(.5, 8., 1.); dividers[0] = .25
        manual, _ = polygon_labels(build_x_polygons(geometry, occupied, raw,
                                                  boundaries_u=dividers), probe)
        self.assertEqual(default.tolist(), [0])
        self.assertEqual(manual.tolist(), [1])

    def test_y_discovery_and_root_polygon_contract(self):
        raw, truth = y_fixture()
        model = fit_y(raw)
        self.assertEqual(model['band_count'], 5)
        numeric = predict_y(model, raw)
        self.assertGreater(agreement_up_to_permutation(numeric, truth), .995)
        labels, counts = polygon_labels(build_y_polygons(model, raw), raw)
        self.assertLessEqual(int(counts.max()), 1, 'Cuts mix distinct provisional bands')
        self.assertGreater(np.mean(labels == numeric), .995)

    def test_y_decision_region_preserves_two_holes(self):
        # Two compact clouds carve disjoint holes out of a broad competitor.
        # Rendering-path simplification must not change selection topology.
        model = dict(theta_deg=0., mean=[0., 0.], scale=[1., 1.],
                     means=[[-1., 0.], [1., 0.], [0., 0.]],
                     covariances=[[[.05, 0.], [0., .05]],
                                  [[.05, 0.], [0., .05]],
                                  [[4., 0.], [0., 4.]]],
                     weights=[.15, .15, .7], band_count=3)
        rng = np.random.RandomState(121)
        raw = np.r_[rng.uniform(-5., 5., (12000, 2)),
                    [[-1., 0.], [1., 0.], [0., 3.], [0., 0.]]]
        numeric = predict_y(model, raw)
        labels, counts = polygon_labels(build_y_polygons(model, raw, grid=801), raw)
        self.assertLessEqual(int(counts.max()), 1)
        self.assertGreater(np.mean(labels == numeric), .998)
        np.testing.assert_array_equal(labels[-4:], [0, 1, 2, 2])

    def test_frozen_models_allow_sparse_targets_and_manual_overrides(self):
        from run_geometric_selection import fit_slice, plain
        x, _, geometry = x_fixture()
        occupancy = detect_x(geometry, x)
        occupancy = {key: occupancy[key] for key in ('accepted_indices', 'mixture_model', 'candidates')}
        y, _ = y_fixture(); ymodel = fit_y(y)
        # Sparse targets inherit reviewed geometry rather than rediscovering it.
        data = dict(xpfp=x[:12, 0], xfp=x[:12, 1],
                    ypfp=y[:12, 0], yfp=y[:12, 1], delta=np.full(12, -2.))
        select = np.ones(12, bool)
        with tempfile.TemporaryDirectory() as name:
            base = Path(name)
            (base/'x.json').write_text(json.dumps(plain(dict(axis='x', geometry=geometry,
                                                            occupancy=occupancy))))
            (base/'y.json').write_text(json.dumps(plain(dict(axis='y', model=ymodel))))
            dividers = np.arange(.5, 8., 1.).tolist(); dividers[0] = .25
            with patch('run_geometric_selection.fit_x', side_effect=AssertionError('Unexpected X refit')), \
                 patch('run_geometric_selection.detect_x', side_effect=AssertionError('Unexpected occupancy refit')):
                _, model, prediction, _, _, _ = fit_slice(
                    'x', data, select, select, 0, [-5., 0.], 9,
                    dict(frozen_model='x.json', boundaries_u=dividers), base)
            self.assertEqual(model['target_events'], 12)
            np.testing.assert_array_equal(prediction['labels'], np.zeros(12, int))
            raw = y[:12]; lo, hi = raw.min(0)-1., raw.max(0)+1.
            points = [lo.tolist(), [hi[0], lo[1]], hi.tolist(), [lo[0], hi[1]], lo.tolist()]
            (base/'manual.json').write_text(json.dumps([dict(index=4, points=points)]))
            with patch('run_geometric_selection.fit_y', side_effect=AssertionError('Unexpected Y refit')):
                _, model, prediction, _, _, automatic = fit_slice(
                    'y', data, select, select, 0, [-5., 0.], 9,
                    dict(frozen_model='y.json', polygons_file='manual.json'), base)
            self.assertEqual(model['target_events'], 12)
            self.assertTrue(np.all(prediction['labels'] == 4))
            self.assertTrue(np.any(automatic != prediction['labels']))

    def test_bad_inputs_fail_explicitly(self):
        for raw in (np.zeros((99, 2)), np.ones((200, 2)),
                    np.full((200, 2), np.nan)):
            with self.assertRaises(ValueError):
                fit_x(raw)
            with self.assertRaises(ValueError):
                fit_y(raw)

    def test_manual_polygons_bypass_failed_discovery(self):
        from run_geometric_selection import fit_slice
        # Degenerate tiny sample cannot initialize either automatic engine.
        data = {key: np.zeros(3) for key in ('xpfp', 'xfp', 'ypfp', 'yfp', 'delta')}
        select = np.ones(3, bool)
        points = [[-1., -1.], [1., -1.], [1., 1.], [-1., 1.], [-1., -1.]]
        with tempfile.TemporaryDirectory() as name:
            base = Path(name)
            (base/'manual.json').write_text(json.dumps([dict(index=2, points=points)]))
            for axis in ('x', 'y'):
                with patch('run_geometric_selection.fit_x', side_effect=AssertionError('Unexpected X fit')), \
                     patch('run_geometric_selection.fit_y', side_effect=AssertionError('Unexpected Y fit')):
                    _, model, prediction, _, _, automatic = fit_slice(
                        axis, data, select, select, 0, [-5., 0.], 9,
                        dict(polygons_file='manual.json'), base)
                self.assertTrue(model['manual_partition'])
                np.testing.assert_array_equal(prediction['labels'], [2, 2, 2])
                np.testing.assert_array_equal(automatic, [-1, -1, -1])


if __name__ == '__main__':
    unittest.main()
