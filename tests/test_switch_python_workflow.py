"""Exercise cleanup and target diagnostics with synthetic data on both arms."""
from array import array
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
import ROOT
REPO=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(REPO),str(REPO/'diagnostics/conditioning')]
from spectrometer_config import from_campaign
from preliminary_angular_conditioning import selected_indices
from angular_solver_comparison import reconstruct_vertex
ROOT.gROOT.SetBatch(True)

def run(command,cwd):
    r=subprocess.run(command,cwd=cwd,capture_output=True,text=True,timeout=90,env=dict(os.environ,MPLBACKEND='Agg',MPLCONFIGDIR='/tmp/hallc-test-mpl'))
    assert r.returncode==0,r.stdout+r.stderr
    return r

def main():
  with tempfile.TemporaryDirectory(prefix='hallc-python-switch-') as name:
    temp=Path(name)
    for arm in ('HMS','SHMS'):
      spec=from_campaign(arm+'_same_target');campaign=temp/(arm+'_same_target');campaign.mkdir();(campaign/'config').mkdir()
      (campaign/'config/rungroups_test_inputs.tsv').write_text('rungroup\toptics_id\tangle_deg\tfoils\tnruns\truns\trootfile\ntest\t99001\t12.5\t-10,0,10\t1\t99001\tunused.root\n')
      arrays=dict(delta=np.array([1.,1.,1.]),ztarT=np.array([-10.,0.,10.]),ysT=np.array([spec.ys(spec.ny-1)]*3))
      indices,counts=selected_indices(arrays,[-10,0,10],[-10,0,5,10,15,20,22],0,100,spec)
      assert list(indices)==[0,1,2] and list(counts)==[1,1,1]
      if arm=='SHMS':
        z=np.array([-10.,0.,10.]);yp=np.array([-.02,.01,.03]);a=math.radians(12.5);rx=.27
        y=-z*np.sin(a)+rx*np.cos(a)-spec.ymp-yp*(z*np.cos(a)+rx*np.sin(a))
        prediction=np.column_stack((np.zeros(3),y/100,yp))
        aux=dict(xtar_prediction=np.zeros((3,3)),angle=np.full(3,a),ymis=np.full(3,spec.ymp),xbeam=np.full(3,rx),spectrometer=arm)
        assert np.allclose(reconstruct_vertex(prediction,aux),z,atol=1e-12)
      for axis in ('y','x'):
        folder=campaign/('04a_candidate_trees_y' if axis=='y' else '04b_candidate_trees_x')/'root';folder.mkdir(parents=True)
        path=folder/(('Yscol' if axis=='y' else 'Xscol')+'Candidates_test.root')
        f=ROOT.TFile(str(path),'RECREATE');t=ROOT.TTree('TYCand' if axis=='y' else 'TXCand','Synthetic candidate')
        values={}
        columns=['entry','run','foil','ndel',axis+'scol']
        suffixes=['dc.x_fp','dc.xp_fp','dc.y_fp','dc.yp_fp','extcor.xsieve','extcor.ysieve','gtr.dp','gtr.y','gtr.th','gtr.ph']
        for field in columns+['delta_low','delta_high']+[spec.branch(s) for s in suffixes]:
          code='i' if field in columns else 'd';values[field]=array(code,[0]);t.Branch(field,values[field],field+('/I' if code=='i' else '/D'))
        rng=np.random.default_rng(13)
        for i in range(60):
          for key,value in [('entry',i),('run',99001),('foil',2),('ndel',6),(axis+'scol',spec.ny-1),('delta_low',15),('delta_high',20),(spec.branch('gtr.dp'),16)]: values[key][0]=value
          values[spec.branch('extcor.ysieve')][0]=spec.ys(spec.ny-1)+rng.normal(0,.03)
          values[spec.branch('extcor.xsieve')][0]=spec.xs(spec.nx-1)+rng.normal(0,.03)
          for suffix in suffixes[:4]:values[spec.branch(suffix)][0]=rng.normal(0,.01)
          t.Fill()
        t.Write();f.Close()
        script=REPO/f'gmm_cleanup_{axis}scol_candidates_v2.py'
        run([sys.executable,str(script),'--campaign',str(campaign),'--rungroup','test','--foil','2','--ndel','6','--n-init','1'],temp)
        outfolder=campaign/('05a_gmm_cleanup_y' if axis=='y' else '05b_gmm_cleanup_x')
        output=next((outfolder/'root').glob('*.root'));f=ROOT.TFile(str(output));t=f.Get('GMMClean')
        assert t.GetEntries()==60 and str(f.Get('hallc_spectrometer').GetTitle())==arm
        assert all(int(row.__getattr__(axis+'scol'))==spec.ny-1 for row in t)
        f.Close()
      # Plot the same triple target for both arms; includes all three foil positions.
      tsv=campaign/'peaks.tsv';tsv.write_text('run\topticsID\tnominal_foil_z_cm\tfit_mean_cm\tfit_sigma_cm\tmean_minus_nominal_cm\tentries\n'+''.join(f'99001\ttest\t{z}\t{z+.1}\t.3\t.1\t100\n' for z in (-10,0,10)))
      run([sys.executable,str(REPO/'diagnostics/validation/ztar/plot_ztar_stability.py'),str(tsv),'--outdir',str(campaign/'plots')],temp)
      import pandas as pd
      summary=pd.read_csv(campaign/'plots/peaks_pair_summary.tsv',sep='\t')
      assert 20 in summary['expected_separation_cm'].values
      assert len(pd.read_csv(campaign/'plots/peaks_fit_summary_errors.tsv',sep='\t'))==3
    print('PASS: both-arm cleanup with raw branch fallbacks; final sieve column; independent three-foil selection; SHMS vertex closure; target plots include -10/0/+10 for both arms')
if __name__=='__main__':main()
