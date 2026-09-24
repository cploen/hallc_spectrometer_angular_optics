"""Numerical and sample isolation checks for smooth Huber regression."""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.optimize import check_grad, least_squares
from scipy.special import pseudo_huber
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from huber_optics import FixedBasis,loss_gradient,robust_scale
from fit_huber import sample_ladder


class HuberTests(unittest.TestCase):
    def test_loss_and_gradient(self):
        rng=np.random.default_rng(3)
        q=rng.normal(size=(40,5));y=rng.normal(size=40);b=rng.normal(size=5)
        value,gradient=loss_gradient(b,q,y,1.5)
        self.assertAlmostEqual(value,float(np.mean(pseudo_huber(1.5,q@b-y))),places=13)
        self.assertLess(check_grad(lambda v:loss_gradient(v,q,y,1.5)[0],
                                   lambda v:loss_gradient(v,q,y,1.5)[1],b),1e-6)

    def test_squared_matches_direct_svd_and_scale_invariance(self):
        rng=np.random.default_rng(4)
        x=np.column_stack((np.ones(1000),rng.normal(size=(1000,4))))
        x[:,4]=x[:,1]+1e-4*x[:,4]
        x[:,2]*=1e-8
        y=x@np.array([.8,2.,4e7,.2,-.1])+rng.normal(0,.1,len(x))
        factor=FixedBasis(x)
        coef,_=factor.fit(y)
        reference=np.linalg.lstsq(x,y,rcond=1e-14)[0]
        np.testing.assert_allclose(x@coef,x@reference,atol=1e-7)
        h,_=factor.fit(y,.1)
        xp=x.copy();xp[:,2]*=1e8
        hp,_=FixedBasis(xp).fit(y,.1)
        np.testing.assert_allclose(x@h,xp@hp,atol=1e-7)

    def test_outliers_and_intercept(self):
        rng=np.random.default_rng(8)
        x=np.column_stack((np.ones(1000),rng.normal(size=1000)))
        truth=x@np.array([4.,2.]);y=truth+rng.normal(0,.15,len(x))
        y[::10]+=30
        factor=FixedBasis(x);sq,_=factor.fit(y);h,info=factor.fit(y,.15)
        self.assertLess(np.mean((x@h-truth)**2),.01*np.mean((x@sq-truth)**2))
        self.assertLess(abs(h[0]-4),.1)
        self.assertTrue(info['converged'])
        with self.assertRaises(ValueError):factor.fit(y,0)
        with self.assertRaises(ValueError):robust_scale(np.ones(20))

    def test_rank_deficient_against_independent_solver(self):
        rng=np.random.default_rng(21)
        x=np.column_stack((np.ones(300),rng.normal(size=(300,3))))
        x=np.column_stack((x,2*x[:,1]))
        y=3+x[:,2]*.7+rng.normal(0,.2,len(x));y[::13]+=8
        factor=FixedBasis(x)
        self.assertEqual(factor.rank,4)
        h,_=factor.fit(y,.2)
        reference=least_squares(lambda b:(x@b-y)/.2,np.zeros(5),
                                loss='soft_l1',f_scale=1.5,ftol=1e-12,gtol=1e-12,xtol=1e-12)
        np.testing.assert_allclose(x@h,x@reference.x,atol=2e-6)

    def test_nested_sample_and_protection(self):
        n=200
        a=dict(entry=np.arange(n),rungroup=np.full(n,'rg01'),sample=np.zeros(n,int),
               quality=np.tile([0,1,2,2],50),excluded=np.zeros(n,bool),balanced=np.zeros(n,bool))
        a['sample'][::7]=2;a['excluded'][::9]=True
        core=(a['quality']==2)&(a['sample']!=2)&~a['excluded']
        a['balanced'][np.flatnonzero(core)[:10]]=True
        ladder=sample_ladder(a)
        old=np.zeros(n,bool)
        for h in ladder.values():
            self.assertFalse(np.any(old&~h))
            self.assertFalse(np.any(h&((a['sample']==2)|a['excluded']|(a['quality']==0))))
            old=h
        np.testing.assert_array_equal(ladder['all_core'],core)


if __name__=='__main__':unittest.main()
