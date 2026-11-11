"""Conditional-fit intercept, coordinate contract and frozen export checks."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import uproot
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from test_elastic import fixture
from core_sample import digest, write_tsv
from elastic_net import DEFAULTS, scaling, transform, tune, macro_mse, SCALES
from beam_search import Folds, fit_all
from fit_elastic import load_sample, matrix_rows
import fit_xtar


def geometric_fixture(root):
    campaign=fixture(root)
    for sample in ('fit','holdout'):
        build=campaign/'06c_core_ntuple/balanced'/sample
        path=build/'root/Optics_999_-1_fit_tree_gmm.root'
        with uproot.open(path) as f:a=f['TFit'].arrays(library='np')
        a['xtarT']=a['xtar'].copy()
        with uproot.recreate(path) as f:f.mktree('TFit',a)
        manifest=json.loads((build/'build.json').read_text())
        manifest['outputs'][str(path.relative_to(build))]=digest(path)
        (build/'build.json').write_text(json.dumps(manifest))
    return campaign


class XtarTests(unittest.TestCase):
    def test_uncentered_cv_matches_direct_x_without_intercept(self):
        rng=np.random.default_rng(71)
        x=rng.normal(size=(120,4))+np.array([3,1,2,4])
        y=x@rng.normal(size=(4,3))+.1*rng.normal(size=(120,3))
        folds=np.arange(120)%3;cells=np.array(['a']*90+['b']*30)
        prepared=Folds(x,y,folds,cells,1e-12,1,fit_intercept=False)
        for selected in ((),(0,2),tuple(range(4))):
            prediction,stats=prepared.predict(0,selected)
            expected=np.zeros(len(y))
            for f in range(3):
                tr,va=folds!=f,folds==f
                state=scaling(x[tr],y[tr],fit_intercept=False)
                z,t=transform(x[tr],y[tr],state)
                if selected:
                    b=np.linalg.lstsq(z[:,selected],t[:,0],rcond=1e-12)[0]
                    expected[va]=(x[va][:,selected]/state[1][list(selected)])@b*state[3][0]
            np.testing.assert_allclose(prediction,expected,atol=1e-11)
            self.assertAlmostEqual(prepared.evaluate(0,selected)['score'],
                float(macro_mse((expected-y[:,0])*SCALES[0],cells)),places=5)
            self.assertEqual(stats[0]['terms'],len(selected))
        coefficients,_=fit_all(x,y,[tuple(range(4))]*3,1e-12,fit_intercept=False)
        np.testing.assert_array_equal(coefficients[0],0)
        cfg=dict(DEFAULTS,alphas=[.01],l1=[.5],max_iter=2000)
        result=tune(x,y,folds,cells,cfg,fit_intercept=False)
        np.testing.assert_array_equal(result['enet'][0],0)
        np.testing.assert_array_equal(result['svd'][0],0)

    def test_geometric_coordinate_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            campaign=geometric_fixture(Path(tmp))
            a,_=load_sample(campaign,'balanced','fit',xtar_coordinate='xtarT')
            np.testing.assert_allclose(a['xtarT'],a['xtar'])
            path=campaign/'06c_core_ntuple/balanced/fit/root/Optics_999_-1_fit_tree_gmm.root'
            with uproot.open(path) as f:values=f['TFit'].arrays(library='np')
            values['xtarT']+=.1
            with uproot.recreate(path) as f:f.mktree('TFit',values)
            bm=path.parent.parent/'build.json';m=json.loads(bm.read_text())
            m['outputs'][str(path.relative_to(bm.parent))]=digest(path);bm.write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError,'inconsistent geometric xtarT'):
                load_sample(campaign,'balanced','fit',xtar_coordinate='xtarT')

    def test_pipeline_preserves_frozen_blocks_and_ids(self):
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            campaign=geometric_fixture(Path(tmp))
            a,_=load_sample(campaign,'balanced','fit',xtar_coordinate='xtarT')
            folds=campaign/'folds.tsv'
            write_tsv(folds,[dict(rungroup=r,entry=int(e),fold=i%3)
                for i,(r,e) in enumerate(zip(a['rungroup'],a['entry']))])
            args=SimpleNamespace(campaign=campaign,tag='balanced',name='xtar_test',
                baseline=campaign/'06c_core_ntuple/balanced/fit/oldfit.dat',folds_from=folds,
                max_order=4,beam=2,steps=2,threads=1,check=False)
            def stub(out,a,predictions,label):
                (out/f'{label}_ids.json').write_text(json.dumps([int(i) for i in a['entry']]))
            with patch('elastic_diagnostics.training_plots'),patch('fit_xtar.diagnostics',stub):
                fit_xtar.run(args)
                fit_xtar.evaluate(args)
            out=campaign/'06g_xtar_refit/xtar_test'
            manifest=json.loads((out/'manifest.json').read_text())
            self.assertEqual(manifest['status'],'complete')
            self.assertFalse(manifest['fit_intercept'])
            old={t:c for c,t in matrix_rows(args.baseline)}
            for name in ('enet','full','start','beam'):
                for c,t in matrix_rows(out/f'matrices/{name}.dat'):
                    self.assertEqual(c[3],old[t][3])
                    if t[4]==0:np.testing.assert_array_equal(c,old[t])
            train=set(json.loads((out/'development_cv_ids.json').read_text()))
            hold=set(json.loads((out/'evaluation/reserved_holdout_ids.json').read_text()))
            self.assertEqual(len(train),480);self.assertEqual(len(hold),450)
            self.assertFalse(train&hold)


if __name__=='__main__':unittest.main()
