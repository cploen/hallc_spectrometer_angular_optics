import unittest
import numpy as np
from xtar_correction import Reconstruction,Metrics,guards,scaled_fit,validate_geometry,features


class CorrectionTests(unittest.TestCase):
    def sample(self,n=1000):
        rng=np.random.default_rng(19)
        return dict(entry=np.arange(n),angle=np.full(n,np.deg2rad(15)),
                    reactycalc=rng.normal(0,.02,n),rz=rng.uniform(-10,10,n),
                    xfp=rng.normal(0,10,n),yfp=rng.normal(0,5,n),
                    xpfp=rng.normal(0,.02,n),ypfp=rng.normal(0,.02,n),
                    quality=np.full(n,2),ztarT=np.where(np.arange(n)<n/2,-10.,10.),
                    delta=np.where(np.arange(n)<n/2,-6.,6.))

    def test_zero_matches_existing_afterburner(self):
        from sieve_afterburner import reconstruct
        a=self.sample();a['ry']=a['reactycalc']
        rows=[(np.array([.001,.002,-.001,0]),(0,0,0,0,0)),
              (np.array([.2,.1,.01,0]),(1,0,0,0,0)),
              (np.array([.3,.01,.02,0]),(0,0,0,0,1))]
        r=Reconstruction(a,rows)
        _,_,act,v,x=reconstruct(a,rows,15,[0,0,0],return_target=True)
        expected=np.c_[v[:,0]*1000,v[:,2]*1000,v[:,1]*1000,x*10]
        np.testing.assert_allclose(r.base,expected,rtol=1e-14,atol=1e-12)
        np.testing.assert_array_equal(r.base_limited,act)

    def test_correction_does_not_accumulate(self):
        a=self.sample();rows=[(np.array([0.,0.,0.,0.]),(0,0,0,0,0))]
        r=Reconstruction(a,rows);v,_=r.predict(np.full(len(a['entry']),.03))
        np.testing.assert_allclose(v[:,3]/10-r.anchor,.03,rtol=0,atol=1e-15)

    def test_truth_geometry_and_mispointing(self):
        a=self.sample();a['xsT']=np.ones(len(a['entry']))*2
        anchor=-a['reactycalc']-.1*(2.37-.086*15+.0012*15**2)
        z=a['ztarT']*np.cos(a['angle']);a['xptarT']=(a['xsT']-anchor)/(168-z)
        a['xtarT']=anchor-a['xptarT']*z
        validate_geometry(a)
        a['xtarT']+=.01
        with self.assertRaises(AssertionError):validate_geometry(a)

    def test_cell_loss_cannot_hide_in_overall(self):
        a=self.sample();rng=np.random.default_rng(2);base=rng.normal(size=(1000,4));candidate=base.copy()
        candidate[:500,0]*=1.05;candidate[500:,0]*=.8
        g=guards(Metrics(a),base,candidate)
        self.assertFalse(g['passed']);self.assertTrue(g['width_violations'])

    def test_mean_shift_not_width_gain(self):
        a=self.sample();base=np.random.default_rng(8).normal(size=(1000,4));shift=base+10
        m=Metrics(a)
        self.assertAlmostEqual(m.score(base),m.score(shift))
        self.assertFalse(guards(m,base,shift)['passed'])

    def test_direct_ridge_scaling_and_zero_intercept(self):
        x=np.random.default_rng(6).normal(size=(300,3))*[1,.001,30];y=x@np.array([.2,30,-.02])
        b,info=scaled_fit(x,y,(0,1,2),'ridge',.03)
        z=x/np.sqrt(np.mean(x*x,axis=0))
        expected=np.linalg.lstsq(np.r_[z,np.eye(3)*3],np.r_[y,np.zeros(3)],rcond=1e-12)[0]
        np.testing.assert_allclose(b,expected/np.sqrt(np.mean(x*x,axis=0)))
        np.testing.assert_array_equal(np.zeros((1,3))@b,[0.])
        self.assertEqual(info['rank'],3)

    def test_features_ignore_truth_labels(self):
        a=self.sample();x,n=features(a);a['ztarT']*=2;a['delta']*=2;a['xtarT']=np.ones(len(x))*999
        xx,_=features(a);np.testing.assert_array_equal(x,xx);self.assertEqual(x.shape[1],20)


if __name__=='__main__':unittest.main()
