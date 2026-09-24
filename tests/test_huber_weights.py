"""Population weights outside robust loss, weighted-SVD parity and foil totals."""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.optimize import check_grad,least_squares
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from huber_optics import FixedBasis,loss_gradient
from fit_huber_balanced import foil_weights,weight_audit


class WeightedHuberTests(unittest.TestCase):
    def test_gradient(self):
        rng=np.random.default_rng(50);q=rng.normal(size=(100,4));y=rng.normal(size=100)
        b=rng.normal(size=4);w=np.exp(rng.normal(size=100))
        self.assertLess(check_grad(lambda v:loss_gradient(v,q,y,1.5,w)[0],
                                   lambda v:loss_gradient(v,q,y,1.5,w)[1],b),2e-6)

    def test_weighted_svd_and_duplicate_observation_equivalence(self):
        rng=np.random.default_rng(20);x=np.column_stack((np.ones(200),rng.normal(size=(200,3))))
        y=x@np.array([3.,.7,-.2,1.1])+rng.normal(0,.3,len(x));y[::11]+=8.
        w=rng.integers(1,5,len(x))
        f=FixedBasis(x,weights=w);sq,_=f.fit(y)
        direct=np.linalg.lstsq(x*np.sqrt(w[:,None]),y*np.sqrt(w),rcond=1e-12)[0]
        np.testing.assert_allclose(x@sq,x@direct,atol=1e-10)
        h,_=f.fit(y,.3)
        xd=np.repeat(x,w,axis=0);yd=np.repeat(y,w)
        ref=least_squares(lambda b:(xd@b-yd)/.3,np.zeros(4),loss='soft_l1',f_scale=1.5,
                          ftol=1e-12,gtol=1e-12,xtol=1e-12)
        np.testing.assert_allclose(x@h,x@ref.x,atol=2e-6)
        scaled,_=FixedBasis(x,weights=100*w).fit(y,.3)
        np.testing.assert_allclose(x@h,x@scaled,atol=1e-8)

    def test_equal_foil_base_and_robust_accounting(self):
        z=np.repeat([-8.,-3.,0.,3.,8.],[5,10,50,7,8]);w=foil_weights(z)
        self.assertAlmostEqual(w.mean(),1.)
        for v in np.unique(z): self.assertAlmostEqual(w[z==v].sum()/w.sum(),.2)
        q=np.where(np.arange(len(z))%3,2,1);r=np.where(q==2,.1,8.)
        rows=weight_audit(dict(zfoil=z,quality=q),w,r,1.,1.5,'test','xptar',True)
        foil=[v for v in rows if v['group_kind']=='foil']
        self.assertAlmostEqual(sum(v['irls_share'] for v in foil),1.)
        shoulder=next(v for v in rows if v['group_kind']=='quality' and v['group']=='shoulder')
        self.assertLess(shoulder['irls_share'],shoulder['base_share'])

    def test_invalid_weights(self):
        x=np.column_stack((np.ones(10),np.arange(10)))
        for w in (np.zeros(10),-np.ones(10),np.ones(9),np.full(10,np.nan)):
            with self.assertRaises(ValueError):FixedBasis(x,weights=w)
        with self.assertRaises(ValueError):foil_weights([])


if __name__=='__main__':unittest.main()
