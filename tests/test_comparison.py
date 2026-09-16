"""Targeted checks for frozen evaluation, physical units and historical membership."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from compare_matrices import MODELS,read_matrix,predict,offset_rows,apply_offsets,historical_overlap,run,replot
from comparison_plots import report,stats,tail_percent,rms_change,grouped_statistics,TARGETS
from elastic_diagnostics import groups
from core_sample import read_tsv,digest


class ComparisonTest(unittest.TestCase):
    def test_grouped_statistics_parity_and_selection_reuse(self):
        # Interleaved pools/settings and sparse holes, including the N=10 boundary.
        rng=np.random.default_rng(667);n=120
        pool=np.resize(['protected_core','protected_noncore','surplus_core',
                        'protected_unsupported','blocked'],n)
        a=dict(ztarT=np.resize([-8.,0.,8.],n),ndel=np.arange(n)%2,
               rungroup=np.resize(['rg1','rg2'],n),xscol=np.arange(n)%4,yscol=np.arange(n)%3)
        rr={m:rng.normal(size=(n,4)) for m in MODELS}
        expected=[];ngroups=0
        for g,sel in groups(a,pool):
            ngroups+=1
            for j,t in enumerate(TARGETS):
                for m in MODELS:
                    expected.append(dict(g,model=m,target=t,qa_low=int(sum(sel))<10,
                                         **stats(rr[m][sel,j])))
        class CountSelections:
            def __init__(self,values):self.values=values;self.calls=0
            def __getitem__(self,index):
                # A single integer-index gather per group/model, never a mask per target.
                self.calls+=1
                assert isinstance(index,np.ndarray) and index.dtype.kind in 'iu'
                return self.values[index]
        counted={m:CountSelections(v) for m,v in rr.items()}
        self.assertEqual(grouped_statistics(a,pool,counted),expected)
        self.assertTrue(all(v.calls==ngroups for v in counted.values()))
        # Explicit single-hole groups at N=1, 9, 10, plus blocked-only input.
        for size in (1,9,10):
            data={k:v[:size]*0 if v.dtype.kind!='U' else np.full(size,'rg') for k,v in a.items()}
            rows=grouped_statistics(data,np.full(size,'protected_core'),{m:v[:size] for m,v in rr.items()})
            self.assertTrue(all(r['n']==size and r['qa_low']==(size<10) for r in rows))
        self.assertEqual(grouped_statistics(a,np.full(n,'blocked'),rr),[])

    def test_matrix_comments_units_and_offsets(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'matrix.dat'
            path.write_text('! 99 99 99 0 00000\n---\n1D-3 2D-3 3D-3 0 00000\n1 2 3 0 10000\n---\n')
            rows=read_matrix(path)
            a=dict(entry=np.arange(2),rungroup=np.array(['rg','rg']),xfp=np.array([0.,1.]),
                   xpfp=np.zeros(2),yfp=np.zeros(2),ypfp=np.zeros(2),xtar=np.zeros(2))
            pred=predict(a,{m:rows for m in MODELS},threads=1,chunk=1)
            np.testing.assert_allclose(pred['old'],[[.001,.002,.003],[.011,.022,.033]])
            policy=dict(offsets={m:[0,0,0] for m in MODELS},offset_sources={m:'test' for m in MODELS})
            policy['offsets']['old']=[.711802209,0,.584959117]
            before=pred['old'].copy();apply_offsets(pred,a,offset_rows(a['rungroup'],policy))
            np.testing.assert_allclose(pred['old']-before,np.tile([.000711802209,0,.000584959117],(2,1)))
            np.testing.assert_allclose(pred['beam'],before)
            self.assertAlmostEqual(stats(before[:,0])['spread'],stats(pred['old'][:,0])['spread'])
            policy['offsets']['old']=None
            with self.assertRaisesRegex(ValueError,'unresolved'):offset_rows(['rg'],policy)
            policy['offsets']['old']={'other':[0,0,0]}
            with self.assertRaisesRegex(ValueError,'old/rg'):offset_rows(['rg'],policy)

    def test_corrupt_row_requires_explicit_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'matrix.dat';bad='1e-310 2e-310 3e-310 0 96228528088595472097541840'
            path.write_text('0 0 0 0 00000\n'+bad+'\n')
            with self.assertRaisesRegex(ValueError,'Malformed'):read_matrix(path)
            self.assertEqual(len(read_matrix(path,[bad])),1)
            with self.assertRaisesRegex(ValueError,'not found'):read_matrix(path,[bad,'missing'])

    def test_exact_historical_membership(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'ids.npz';base=2**54
            np.savez(path,**{'999_entry':np.array([base+1],dtype=np.int64)})
            a=dict(entry=np.array([base,base+1,base+1]),rungroup=np.array(['rg','rg','unknown']))
            known,seen=historical_overlap(a,path,[dict(rungroup='rg',optics_id='999')])
            np.testing.assert_array_equal(known,[True,True,False]);np.testing.assert_array_equal(seen,[False,True,False])

    def test_check_does_not_evaluate(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch('compare_matrices.prepare',return_value={'pool':np.array(['protected_core'])}),patch('compare_matrices.predict',side_effect=AssertionError('must not evaluate')):
                run(Path(tmp),'test','compare','beam',{},check=True)
            self.assertFalse((Path(tmp)/'07_diagnostics').exists())

    def test_report_common_subset_and_bias(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);(out/'tsv').mkdir();(out/'plots').mkdir()
            n=24;a=dict(entry=np.arange(n),rungroup=np.full(n,'rg'),ztarT=np.repeat([-8.,0.,8.],8),
                        ndel=np.tile(np.repeat([0,1],4),3),xscol=np.ones(n,int),yscol=np.ones(n,int))
            pool=np.full(n,'protected_core');known=np.ones(n,bool);seen=np.arange(n)%2==0
            residual=np.tile(np.linspace(-1,1,n)[:,None],(1,4));rr={m:residual+(4-i)*.2 for i,m in enumerate(MODELS)}
            with patch('comparison_plots.plots'):report(out,a,pool,rr,known,seen,{'method':'fixture'})
            rows=read_tsv(out/'tsv/summary.tsv');r=next(r for r in rows if r['model']=='beam' and r['target']=='xptar' and r['subset']=='outside_gmm')
            self.assertEqual(int(r['n']),12)
            self.assertAlmostEqual(float(r['rms']),np.sqrt(np.mean(residual[~seen,0]**2)))
            for m in MODELS:
                r=next(r for r in rows if r['model']==m and r['target']=='xptar' and r['subset']=='all')
                self.assertAlmostEqual(float(r['rms'])**2,float(r['bias'])**2+float(r['spread'])**2)
            saved=np.load(out/'residuals.npz');np.testing.assert_array_equal(saved['entry'],a['entry'])

    def test_tail_percent_and_change(self):
        # Ties count at the threshold; zooming must not change the denominator.
        np.testing.assert_allclose(tail_percent(np.array([0.,1.,1.,4.]),np.array([0.,1.,2.,4.,5.])),[100,75,25,25,0])
        self.assertAlmostEqual(rms_change(1.8,2.),-10.)
        self.assertTrue(np.isnan(rms_change(1.,0.)))

    def test_replot_saved_results_without_evaluation(self):
        with tempfile.TemporaryDirectory() as tmp:
            campaign=Path(tmp);out=campaign/'07_diagnostics/compare';(out/'tsv').mkdir(parents=True)
            n=24;a=dict(entry=np.arange(n),rungroup=np.full(n,'rg'),ztarT=np.repeat([-8.,0.,8.],8),
                        ndel=np.tile(np.repeat([0,1],4),3),xscol=np.ones(n,int),yscol=np.ones(n,int))
            pool=np.full(n,'protected_core');known=np.ones(n,bool)
            residual=np.tile(np.linspace(-1,1,n)[:,None],(1,4));rr={m:residual*(1+i*.1) for i,m in enumerate(MODELS)}
            with patch('comparison_plots.plots'):report(out,a,pool,rr,known,~known,{'method':'fixture'})
            hashes={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}
            (out/'manifest.json').write_text(json.dumps(dict(schema='matrix_comparison_v1',status='complete',sample='test',outputs=hashes)))
            with patch('compare_matrices.prepare',side_effect=AssertionError('must not load source events')), \
                 patch('compare_matrices.predict',side_effect=AssertionError('must not evaluate')), \
                 patch('comparison_plots.report',side_effect=AssertionError('must not repeat statistics')):
                replot(campaign,'test','compare')
            for rel,sha in hashes.items():
                if rel!='MATRIX_COMPARISON.md':self.assertEqual(digest(out/rel),sha)
            self.assertEqual(len(list((out/'plots').glob('*.png'))),5)
            manifest=json.loads((out/'manifest.json').read_text())
            for rel,sha in manifest['outputs'].items():self.assertEqual(digest(out/rel),sha)
            (out/'tsv/residuals.tsv').write_text('modified')
            with self.assertRaises(ValueError):replot(campaign,'test','compare')


if __name__=='__main__':unittest.main()
