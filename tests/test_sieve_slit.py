"""Check density ranking, ties, and the disjoint population accounting."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

spec = importlib.util.spec_from_file_location('sieve_slit', Path(__file__).resolve().parents[1] / 'diagnostics/validation/sieve_slit/plot_ytar.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PopulationTests(unittest.TestCase):
    def test_local_ranking_and_ties(self):
        data = dict(quality=np.array([2, 2, 2, 2, 2, 1, 0]),
                    core_score=np.array([.8, .9, .9, .2, .1, .05, np.nan]),
                    entry=np.array([1, 9, 3, 4, 5, 6, 7]),
                    foil=np.zeros(7), ndel=np.array([0, 0, 0, 1, 1, 0, 0]),
                    xscol=np.zeros(7), yscol=np.zeros(7))
        masks = module.populations(data)
        self.assertEqual(np.flatnonzero(masks['dense_half']).tolist(), [1, 2, 3])
        self.assertTrue(np.array_equal(masks['out_of_core'], masks['shoulders'] | masks['unsupported']))
        self.assertTrue(np.all(sum(masks[k].astype(int) for k in ['dense_half', 'outer_core', 'shoulders', 'unsupported']) == 1))
        data['core_score'][:3] = .9
        self.assertEqual(np.flatnonzero(module.populations(data)['dense_half']).tolist(), [0, 2, 3])

    def test_reject_invalid_core_score(self):
        data = dict(quality=np.array([2]), core_score=np.array([np.nan]))
        with self.assertRaises(ValueError):
            module.populations(data)


if __name__ == '__main__':
    unittest.main()
