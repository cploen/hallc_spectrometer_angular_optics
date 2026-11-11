import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from xtar_extended_heatmaps import metrics,MODELS,rms_percent_change


class ExtendedMetricsTests(unittest.TestCase):
    def fixture(self):
        n=7
        a=dict(rungroup=np.array(['r']*n),entry=np.arange(n),delta=np.array([-15,-12.1,-12,-10,14,15,3.]),
               ztarT=np.zeros(n),quality=np.array([2,2,2,2,1,0,0]),excluded=np.array([False]*6+[True]))
        bundles={m:dict(rungroup=a['rungroup'].copy(),entry=a['entry'].copy(),values=np.ones((n,4))*v)
                 for m,v in zip(MODELS,[.02,.01,.005])}
        return a,bundles

    def test_bins_bias_rms_signed_changes_and_population_counts(self):
        a,bundles=self.fixture()
        rows,counts=metrics(a,np.zeros((7,4)),bundles,physical_foils=[-10,0,10],supported=np.array([True]*5+[False,False]))
        def row(model,lo,pop='core',z=0):
            return next(r for r in rows if r['model']==model and r['delta_low']==lo and r['population']==pop and r['foil_cm']==z and r['target']=='xtar')
        r=row(MODELS[2],-15)
        self.assertEqual(r['n'],2)
        self.assertAlmostEqual(r['bias'],.05)
        self.assertAlmostEqual(r['rms_including_bias'],.05)  # includes pure bias
        self.assertAlmostEqual(r['rms_change_vs_oldfit'],-.15)
        self.assertAlmostEqual(r['rms_change_vs_refit'],-.05)
        self.assertEqual(row(MODELS[2],-12)['n'],1)
        self.assertEqual(row(MODELS[2],-10)['n'],1)
        self.assertEqual(row(MODELS[2],14,'supported_noncore')['n'],1)
        self.assertEqual(row(MODELS[2],14,z=-10)['rms_including_bias'],'')
        self.assertEqual(counts['outside_extended_delta'],1)
        self.assertEqual(counts['unsupported_noncore'],1)
        self.assertEqual(counts['blocked'],1)

    def test_previous_percent_change_direction_and_zero_reference(self):
        self.assertAlmostEqual(rms_percent_change(1.2,1.),20.)
        self.assertAlmostEqual(rms_percent_change(.8,1.),-20.)
        self.assertEqual(rms_percent_change(0.,0.),0.)
        self.assertTrue(np.isnan(rms_percent_change(1.,0.)))

    def test_id_and_prediction_failures_are_not_reselection(self):
        a,bundles=self.fixture();bundles[MODELS[2]]['entry']=a['entry'][::-1]
        with self.assertRaisesRegex(ValueError,'IDs/order'):
            metrics(a,np.zeros((7,4)),bundles,physical_foils=[0],supported=np.ones(7,dtype=bool))
        a,bundles=self.fixture();bundles[MODELS[2]]['values'][0,0]=np.nan
        with self.assertRaisesRegex(ValueError,'nonfinite'):
            metrics(a,np.zeros((7,4)),bundles,physical_foils=[0],supported=np.ones(7,dtype=bool))


if __name__=='__main__':unittest.main()
