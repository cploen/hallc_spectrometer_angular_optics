"""Exact out-of-fold pooling, incomplete-case exclusion, and protected-pool isolation."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from core_sample import read_tsv
from elastic_convergence import run
from elastic_pool import Pool
from test_elastic import fixture
import fit_elastic


class PoolTests(unittest.TestCase):
    def test_pooled_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = fixture(Path(tmp))
            with patch('elastic_diagnostics.training_plots'):
                fit_elastic.fit(p,'balanced','enet')
            (p/'config/elastic_refit.json').write_text(json.dumps(dict(alphas=[.003],l1=[.5],budgets=[2,10000])))
            (p/'06c_core_ntuple/balanced/holdout/build.json').unlink()
            run(p,'balanced','pooled',refit=True,pooled=True)
            out = p/'06d_elastic_net/pooled'
            manifest = json.loads((out/'manifest.json').read_text())
            self.assertEqual(manifest['converged_cases'],9)
            self.assertEqual(manifest['pooled_cases'],3)
            self.assertEqual(manifest['config']['fold'],'all')
            values = np.load(out/'residuals.npz')
            self.assertEqual(set(values['fold']),{0,1,2})
            self.assertEqual(len(set(zip(values['rungroup'],values['entry']))),len(values['fold']))
            a,_ = fit_elastic.load_sample(p,'balanced','fit')
            np.testing.assert_array_equal(values['entry'],a['entry'])
            for r in read_tsv(out/'pooled.tsv'):
                key = f"{r['target']}_a{float(r['alpha']):g}_l{float(r['l1']):g}"
                residual = values[key][:,0 if r['model']=='full_svd' else 1]
                self.assertTrue(np.isfinite(residual).all())
                self.assertAlmostEqual(float(r['median']),np.median(np.abs(residual)))
                self.assertAlmostEqual(float(r['p90']),np.percentile(np.abs(residual),90))
                self.assertAlmostEqual(float(r['mse']),np.mean(residual**2))
                self.assertEqual(int(r['n']),len(a['entry']))
            self.assertTrue((out/'plots/plateau.png').exists())
            self.assertTrue((out/'plots/foil_delta.png').exists())
            self.assertFalse((out/'matrices').exists())

    def test_missing_fold_is_not_pooled(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = dict(ztarT=np.zeros(3),ndel=np.zeros(3,dtype=int),rungroup=np.array(['rg']*3),
                     entry=np.arange(3),xscol=np.zeros(3,dtype=int),yscol=np.zeros(3,dtype=int))
            pool = Pool(Path(tmp),a,np.zeros((3,3)),np.arange(3),np.array(['0:0']*3),[],10)
            pool.values[0,.003,.5] = np.array([[1.,2.],[2.,3.],[np.nan,np.nan]])
            self.assertEqual(pool.finish(dict(alphas=[.003],l1=[.5])),0)
            self.assertFalse((Path(tmp)/'pooled.tsv').exists())


if __name__ == '__main__': unittest.main()
