import unittest
import numpy as np
from test_xtar_correction import CorrectionTests
from xtar_correction import Metrics
from xtar_tolerance_scan import Moments,permitted

class ToleranceTests(unittest.TestCase):
    def test_moments_match_direct_statistics(self):
        a=CorrectionTests().sample();a['excluded']=np.zeros(1000,bool);a['quality'][::3]=1
        r=np.random.default_rng(12).normal(size=(1000,4))+[3,-1,.02,.1]
        fast=Moments(a);mean,sd=fast.stats(r);slow=Metrics(a).measure(r)
        for i,label in enumerate(fast.labels):
            if label in slow:
                self.assertEqual(fast.count[i],slow[label][0])
                np.testing.assert_allclose(mean[i],slow[label][1],atol=1e-12)
                np.testing.assert_allclose(sd[i],slow[label][2],atol=1e-12)
        self.assertAlmostEqual(fast.score((mean,sd)),Metrics(a).score(r))

    def test_cell_and_overall_are_independent(self):
        r=dict(physical=True,**{'global':np.array([.019,.01,0]),'cell':np.array([.049,.03,0]),'bias':np.zeros(3)})
        self.assertFalse(permitted(r,.02,.02));self.assertTrue(permitted(r,.05,.02))
        r['global'][0]=.04
        self.assertFalse(permitted(r,.05,.02));self.assertTrue(permitted(r,.10,.05))
        r['cell'][0]=.11;self.assertFalse(permitted(r,.10,.05))

    def test_bias_and_physical_guards_do_not_relax(self):
        r=dict(physical=True,**{'global':np.zeros(3),'cell':np.zeros(3),'bias':np.array([.021,0,0])})
        self.assertFalse(permitted(r,.10,.05));r['bias'][:]=0;r['physical']=False
        self.assertFalse(permitted(r,.10,.05))

if __name__=='__main__':unittest.main()
