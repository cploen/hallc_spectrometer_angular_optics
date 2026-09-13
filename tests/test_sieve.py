"""Afterburner polynomial parity, fixed selection, and six-plot replay integration."""
import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import uproot
from core_sample import write_tsv
from compare_matrices import read_matrix
from fit_elastic import design
from sieve_afterburner import coefficients,evaluate,reconstruct,selected,run,FIELDS,DEFAULTS


class SieveTest(unittest.TestCase):
    def test_polynomial_and_projection(self):
        rng=np.random.default_rng(3);n=20
        a={k:rng.normal(size=n) for k in ('xfp','xpfp','yfp','ypfp','ry','rz')}
        a['entry']=np.arange(n);a['xtar']=rng.normal(size=n)
        rows=[(np.array([.1,.2,.3,0]),t) for t in [(1,0,0,0,0),(0,1,1,0,2),(0,0,0,0,0)]]
        expected=design(a,[t for c,t in rows])@np.array([c[:3] for c,t in rows])
        np.testing.assert_allclose(evaluate(coefficients(a,rows),a['xtar'],[0,0,0]),expected)
        constant=[(np.array([.01,.02,.03,0]),(0,0,0,0,0))]
        xy,it,active=reconstruct(a,constant,12.49,[0,0,0],xmis=.1)
        np.testing.assert_allclose(xy[:,0],-a['ry']-.1-.01*a['rz']*np.cos(np.deg2rad(12.49))+1.68)
        np.testing.assert_allclose(xy[:,1],7.04);self.assertTrue(np.all(it==1));self.assertFalse(active.any())

    def test_cuts(self):
        a=dict(cer=np.array([3,3,2,3,3]),cal=np.ones(5),delta=np.array([0,10,0,0,0]),rz=np.array([-8,-8,-8,0,8]))
        good,foil=selected(a,[-8,8],DEFAULTS)
        np.testing.assert_array_equal(good,[True,False,False,False,True]);self.assertEqual(foil[-1],1)

    def test_replay_to_six_plots(self):
        with tempfile.TemporaryDirectory() as tmp:
            campaign=Path(tmp)/'HMS_test';config=campaign/'config';config.mkdir(parents=True)
            mat='0.1 0 0 0 10000\n0 0.2 0 0 00100\n0 0 .1 0 00010\n.01 0 0 0 00001\n'
            (config/'oldfit.dat').write_text(mat);(campaign/'gmm.dat').write_text(mat)
            dest=campaign/'06e_beam_search/beam10/matrices';dest.mkdir(parents=True);(dest/'beam.dat').write_text(mat+'0.001 0 0 0 00000\n')
            (config/'comparison.json').write_text(json.dumps({'gmm':'gmm.dat'}))
            rng=np.random.default_rng(11);table=[]
            for rg,foils in [('rg0',[0]),('rg8',[-8,8])]:
                n=120;a={k:rng.normal(0,.1,n) for k in FIELDS};a.update(cer=np.full(n,5.),cal=np.ones(n),delta=np.zeros(n),rz=np.resize(foils,n)+rng.normal(0,.1,n))
                xy,_,_=reconstruct(a,read_matrix(config/'oldfit.dat'),12.49,[0,0,0]);a['xs']=xy[:,0];a['ys']=xy[:,1]
                path=campaign/f'{rg}.root'
                with uproot.recreate(path) as f:f.mktree('T',{'H.'+v:a[k] for k,v in FIELDS.items()})
                table.append(dict(rungroup=rg,optics_id=1 if rg=='rg0' else 2,hms_angle_deg=12.49,rootfile=str(path)))
            write_tsv(config/'rungroups_test_inputs.tsv',table)
            meta=lambda run,*args:dict(foils=[0] if run==1 else [-8,8])
            with patch('sieve_afterburner.run_metadata',side_effect=meta):
                run(campaign,'check',DEFAULTS,threads=1,check=True)
                self.assertFalse((campaign/'07_diagnostics/check').exists())
                run(campaign,'sieve',DEFAULTS,threads=1)
            out=campaign/'07_diagnostics/sieve';saved=json.loads((out/'manifest.json').read_text())
            self.assertTrue(saved['closure_passed']);self.assertEqual(saved['status'],'complete')
            self.assertEqual(len(list((out/'plots').glob('*.png'))),6)
            h=np.load(out/'histograms.npz')
            for rg,foils in [('rg0',[0]),('rg8',[-8,8])]:
                for z in foils:self.assertEqual(h[f'gmm_{rg}_{z}'].sum(),h[f'beam_{rg}_{z}'].sum())


if __name__=='__main__':unittest.main()
