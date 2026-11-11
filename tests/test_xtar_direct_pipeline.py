"""Frozen-ID, measured-input, anchor and offline transport integration checks."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from test_xtar_refit import geometric_fixture
from core_sample import digest,write_tsv
from fit_elastic import load_sample,matrix_rows
from run_xtar_predictor import run,measured_problem,polynomial_targets
from xtar_predictor import fit_training_design,export_model,predict


class DirectPipelineTests(unittest.TestCase):
    def test_no_truth_features_and_signed_mispoint_anchor(self):
        a=dict(xfp=np.array([10.,20.]),xpfp=np.array([.01,.02]),yfp=np.array([2.,3.]),
               ypfp=np.array([-.02,.01]),reactycalc=np.array([.1,-.1]),angle=np.deg2rad([15.,-30.]))
        x,anchor,features=measured_problem(a)
        self.assertEqual(x.shape,(2,20))
        degree=np.array([15.,30.]);mis=.1*(2.37-.086*degree+.0012*degree**2)
        np.testing.assert_allclose(anchor,-a['reactycalc']-mis)
        b=dict(a,foil=np.array([999,999]),xtarT=np.array([1000.,1000.]),
               ztarT=np.array([10000.,10000.]),xscol=np.array([9,9]),xptarT=np.array([100.,100.]))
        np.testing.assert_array_equal(measured_problem(b)[0],x)
        self.assertEqual(measured_problem(b)[2],features)
        np.testing.assert_allclose(x[:,0],a['xfp']/100*np.cos(a['angle']))

    def test_anchored_export_uses_single_target_and_requires_measured_offset(self):
        rng=np.random.default_rng(91);x=rng.normal(size=(90,3));anchor=rng.normal(size=90)*.01
        truth=anchor+.1+.2*x[:,0]+.002*rng.normal(size=90)
        with contextlib.redirect_stdout(io.StringIO()):
            fit=fit_training_design(x,truth,np.arange(90)%3,np.array(['a']*90),
                known_offset_cm=anchor,en_config=dict(alphas=[.01],l1=[.5],max_iter=2000),beam_config=dict(steps=1,threads=1))
        features=[dict(name=f'x{i}',expression=f'x{i}',unit='m',source='synthetic measured fixture') for i in range(3)]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'model.json'
            with self.assertRaisesRegex(ValueError,'anchor'):
                export_model(p,fit,features=features,provenance={'test':'fixture'},validity={'test':'fixture'})
            model=export_model(p,fit,features=features,provenance={'test':'fixture'},validity={'test':'fixture'},anchor={'unit':'cm','expression':'measured fixture'})
            with self.assertRaisesRegex(ValueError,'anchor'):predict(model,x)
            np.testing.assert_allclose(predict(model,x,anchor_cm=anchor),anchor+np.c_[np.ones(90),x]@fit['coefficients_cm'])

    def test_real_loader_fit_and_reserved_comparison_leave_matrices_fixed(self):
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            c=geometric_fixture(Path(tmp));seed=c/'06c_core_ntuple/balanced/fit/oldfit.dat'
            reference=c/'config/oldfit.dat';reference.write_bytes(seed.read_bytes())
            a,_=load_sample(c,'balanced','fit',xtar_coordinate='xtarT')
            folds=c/'folds.tsv';write_tsv(folds,[dict(rungroup=r,entry=int(e),fold=i%3) for i,(r,e) in enumerate(zip(a['rungroup'],a['entry']))])
            args=SimpleNamespace(campaign=c,tag='balanced',name='direct_test',baseline=seed,
                reference_sha256=digest(reference),baseline_sha256=digest(seed),folds_from=folds,
                threads=1,steps=2,fit=True,evaluate=True)
            before=seed.read_bytes();ref_before=reference.read_bytes()
            captured={}
            def save_metrics(out,rows,counts,**kwargs):
                captured['rows']=rows;captured['counts']=counts
                out.mkdir();(out/'metrics.json').write_text(json.dumps(rows))
            with patch('run_xtar_predictor.write_report',save_metrics):run(args)
            out=c/'06h_xtar_predictor/direct_test'
            manifest=json.loads((out/'manifest.json').read_text())
            self.assertEqual(manifest['status'],'complete');self.assertEqual(manifest['training'],480)
            self.assertEqual(manifest['reserved_counts']['total'],450)
            self.assertEqual(set(r['target'] for r in captured['rows']),{'xptar','yptar','ytar','xtar'})
            self.assertEqual(seed.read_bytes(),before);self.assertEqual(reference.read_bytes(),ref_before)
            self.assertEqual((out/'completed_refit_fixed.dat').read_bytes(),before)
            model=json.loads((out/'xtar_predictor.json').read_text())
            self.assertEqual(model['beam_config']['beam'],30);self.assertEqual(len(model['features']),20)
            with np.load(out/'reserved_predictions.npz') as predictions:
                self.assertEqual(len(predictions['entry']),450)
                self.assertFalse(set(map(int,predictions['entry']))&set(map(int,a['entry'])))
            # Independently check row-column and cm/m conventions, including xtar dependence.
            h=dict(a,xtar=np.full(len(a['entry']),2.))
            got=polynomial_targets(h,matrix_rows(seed),h['xtar'])
            np.testing.assert_allclose(got[:,0],.01*.02)
            np.testing.assert_allclose(got[:,1],.03*.02)
            np.testing.assert_allclose(got[:,2],100*.02*.02)
            np.testing.assert_array_equal(got[:,3],2.)


if __name__=='__main__':unittest.main()
