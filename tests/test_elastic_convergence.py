"""Continuation budgets, saved fold identity, and a training-only study round trip."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from sklearn.linear_model import ElasticNet
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from elastic_convergence import checkpoints, saved_folds, validate, DEFAULTS, run
from test_elastic import fixture
import fit_elastic


class ConvergenceTests(unittest.TestCase):
    def test_continuation_matches_fixed_objective(self):
        rng=np.random.default_rng(5)
        x=np.asfortranarray(rng.normal(size=(200,5)))
        x[:,1]=x[:,0]+rng.normal(0,.05,200)
        x-=x.mean(axis=0);x/=x.std(axis=0)
        y=x[:,0]-.5*x[:,2]+rng.normal(0,.1,200)
        y-=y.mean();y/=y.std()
        gram=np.ascontiguousarray(x.T@x)
        got=list(checkpoints(x,y,gram,.001,.5,[2,6],1e-14))
        self.assertEqual([r['iterations'] for b,r in got],[2,6])
        self.assertFalse(got[-1][1]['converged'])
        with patch('warnings.warn'):
            direct=ElasticNet(alpha=.001,l1_ratio=.5,fit_intercept=False,precompute=gram,max_iter=6,tol=1e-14).fit(x,y)
        np.testing.assert_allclose(got[-1][0],direct.coef_,atol=1e-12,rtol=1e-12)
        self.assertLessEqual(got[-1][1]['objective'],got[0][1]['objective']+1e-14)
        finished=list(checkpoints(x,y,gram,.1,.1,[10000,20000],1e-6))
        self.assertEqual(len(finished),1)
        self.assertTrue(finished[0][1]['converged'])
        self.assertLessEqual(finished[0][1]['gap_ratio'],1)
        with self.assertRaises(ValueError):validate(dict(DEFAULTS,budgets=[20,10]),3)

    def test_saved_ids_and_training_only_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=fixture(Path(tmp))
            with patch('elastic_diagnostics.training_plots'):
                parent=fit_elastic.fit(p,'balanced','enet')
            (p/'config/elastic_conv.json').write_text(json.dumps(dict(alphas=[.003],l1=[.5],budgets=[2,10000])))
            # No holdout export is needed or consulted by the convergence study.
            hp=p/'06c_core_ntuple/balanced/holdout/build.json'; hp.unlink()
            run(p,'balanced','conv')
            out=p/'06d_elastic_net/conv'
            m=json.loads((out/'manifest.json').read_text())
            self.assertEqual(m['status'],'complete');self.assertEqual(m['cases'],3)
            self.assertEqual(m['tol'],1e-6)
            self.assertTrue((out/'convergence.png').exists())
            self.assertFalse((out/'matrices').exists())
            rows=fit_elastic.read_tsv(out/'checkpoints.tsv')
            self.assertEqual({r['fold'] for r in rows},{'0'})
            a,_=fit_elastic.load_sample(p,'balanced','fit')
            folds=fit_elastic.read_tsv(parent/'tsv/folds.tsv')
            v=saved_folds(a,folds,3)
            np.testing.assert_array_equal(v,saved_folds(a,list(reversed(folds)),3))
            with self.assertRaises(ValueError):saved_folds(a,folds[:-1],3)
            with self.assertRaises(FileExistsError):run(p,'balanced','conv')


if __name__=='__main__':unittest.main()
