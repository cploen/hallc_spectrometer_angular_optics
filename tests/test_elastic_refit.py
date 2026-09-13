"""Selected-term direct SVD, fold isolation and comparison-file integration."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from elastic_refit import subset_svd, native
from elastic_convergence import study, run
from elastic_net import scaling, transform
from core_sample import read_tsv
from test_elastic import fixture
import fit_elastic


class RefitTests(unittest.TestCase):
    def test_subset_svd_and_units(self):
        rng=np.random.default_rng(10)
        x=rng.normal(size=(200,4));y=x @ np.array([2.,0,-3.,0])
        selected=np.array([True,False,True,False])
        b,rank,s=subset_svd(x,y,selected,1e-12)
        np.testing.assert_allclose(b,[2,0,-3,0],atol=1e-12)
        self.assertEqual(rank,2)
        empty,rank,s=subset_svd(x,y,np.zeros(4,dtype=bool),1e-12)
        self.assertEqual(rank,0);np.testing.assert_array_equal(empty,0)
        xx=x*np.array([1e-8,1,1e-4,1])+np.array([.003,2,0,4])
        yy=np.column_stack((y,y*.001+5,y*.1-1))
        state=scaling(xx,yy);z,t=transform(xx,yy,state)
        b,_,_=subset_svd(z,t[:,1],selected,1e-12)
        np.testing.assert_allclose(np.c_[np.ones(len(z)),xx] @ native(b,state,1), z @ b*state[3][1]+state[2][1],atol=1e-10)

    def test_validation_does_not_select_or_refit(self):
        rng=np.random.default_rng(15)
        x=rng.normal(size=(180,5));y=x[:,:3]+rng.normal(size=(180,3))*.1
        folds=np.arange(180)%3;cells=np.array(['cell']*180)
        cfg=dict(fold=0,l1=[.5],alphas=[.01],budgets=[1000])
        outputs=[]
        def fit(j,alpha,l1,beta,f):
            coeff,_,_=subset_svd(f['x'],f['y'][:,j],beta!=0,f['rcond'])
            outputs.append((j,beta.copy(),coeff))
        study(x,y,folds,cells,cfg,1e-6,1e-12,lambda row:None,fit)
        original=outputs.copy();outputs.clear()
        x2=x.copy();y2=y.copy();x2[folds==0]*=100;y2[folds==0]+=1000
        study(x2,y2,folds,cells,cfg,1e-6,1e-12,lambda row:None,fit)
        self.assertEqual(len(original),3)
        for first,second in zip(original,outputs):
            np.testing.assert_array_equal(first[1],second[1]);np.testing.assert_array_equal(first[2],second[2])

    def test_comparison_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=fixture(Path(tmp))
            with patch('elastic_diagnostics.training_plots'):
                fit_elastic.fit(p,'balanced','enet')
            (p/'config/elastic_refit.json').write_text(json.dumps(dict(alphas=[.003],l1=[.5],budgets=[2,10000])))
            (p/'06c_core_ntuple/balanced/holdout/build.json').unlink()
            run(p,'balanced','refit',refit=True)
            out=p/'06d_elastic_net/refit'
            m=json.loads((out/'manifest.json').read_text())
            self.assertEqual(m['schema'],'elastic_refit_v1')
            self.assertEqual(m['converged_cases'],3)
            rows=read_tsv(out/'comparison.tsv');self.assertEqual(len(rows),3)
            v=read_tsv(out/'validation.tsv')
            self.assertEqual({r['model'] for r in v},{'full_svd','enet','refit_svd'})
            self.assertIn('ztar',{r['target'] for r in v})
            self.assertIn('hole',{r['level'] for r in v})
            self.assertEqual({r['pool'] for r in v},{'validation'})
            self.assertTrue((out/'plots/comparison.png').exists())
            self.assertTrue(list((out/'plots').glob('*_xptar_holes.png')))
            terms=read_tsv(out/'terms.tsv')
            for r in terms:
                if r['selected']=='False':self.assertEqual(float(r['refit_svd']),0)
                if r['term']=='00000':self.assertEqual(r['selected'],'True')
            self.assertFalse((out/'matrices').exists())


if __name__=='__main__':unittest.main()
