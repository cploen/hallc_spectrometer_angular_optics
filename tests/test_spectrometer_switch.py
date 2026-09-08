"""Synthetic pipeline regression. No replay or physics validation is implied."""
from array import array
import json
import math
import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO))
from spectrometer_config import from_campaign, run_metadata, delta_index
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare(f'#include {json.dumps(str(REPO / "spectrometer_config.h"))}')


def check_profiles():
    for arm in ('HMS','SHMS'):
        py=from_campaign(f'/analysis/{arm}_anything/06a_fit_ntuple/root')
        cpp=ROOT.hallc.profileForName(arm)
        assert (py.prefix,py.nx,py.ny)==(str(cpp.prefix),cpp.nx,cpp.ny)
        for i in range(py.ny):
            assert abs(py.ys(i)-cpp.ys(i))<1e-14
            assert abs(py.xs(i)-cpp.xs(i))<1e-14
        for z in (-10,0,10):
            for dp in (-9,0,16,21):
                for rx,ry,xb in ((0,0,0),(.13,-.09,.07)):
                    angle=12.5; s=math.sin(math.radians(angle)); c=math.cos(math.radians(angle))
                    t=ROOT.hallc.targetTruth(cpp,angle,z,py.xs(py.nx-1),py.ys(0),dp,rx,ry,xb)
                    if arm=='SHMS':
                        zv=z*c+rx*s; xv=-ry-py.xmp; yv=-z*s+rx*c-py.ymp
                        assert abs(t.xtar+t.xptar*zv-xv)<1e-13
                        assert abs(t.ytar+t.yptar*zv-yv)<1e-13
                        assert abs(t.xtar+t.xptar*py.sieve_distance-py.xs(py.nx-1))<1e-13
                        assert abs(t.ytar+t.yptar*py.sieve_distance+cpp.hb(dp)-py.ys(0))<1e-13
                    else:
                        # Independent pre-switch HMS formula, including its beam convention.
                        ymis=.1*(.52-.012*angle+.002*angle*angle)
                        xmis=.1*(2.37-.086*angle+.0012*angle*angle)
                        yp=(py.ys(0)-(z*s-xb*c-ymis))/(168-z*c)
                        xp=py.xs(py.nx-1)/(168-z*c)
                        assert abs(t.yptar-yp)<1e-14
                        assert abs(t.xptar-xp)<1e-14
                        assert abs(t.ytar-(z*(s-yp*c)-xb*(c+yp*s)-ymis))<1e-13
                        assert abs(t.xtar-(-ry-xmis-xp*z*c))<1e-13
    for path in ('unknown','SHMSfoo','HMS_parent/SHMS_child'):
        try: from_campaign(path)
        except ValueError: pass
        else: raise AssertionError(path)


def make_replay(path):
    f=ROOT.TFile(str(path),'RECREATE'); t=ROOT.TTree('T','Synthetic dual-arm input')
    values={}
    suffixes=['cal.etottracknorm','gtr.dp','gtr.x','gtr.y','gtr.th','gtr.ph','dc.x_fp','dc.xp_fp','dc.y_fp','dc.yp_fp','react.x','react.y','react.z','extcor.xsieve','extcor.ysieve','rb.raster.fr_xbpm_tar','rb.raster.fr_ybpm_tar']
    for prefix,cer in (('H','cer'),('P','ngcer')):
        for suffix in suffixes+[cer+'.npeSum']:
            name=prefix+'.'+suffix; values[name]=array('d',[0]); t.Branch(name,values[name],name+'/D')
    rows=[]
    for foil,z in enumerate((-10,0,10)):
        for dp,nd in ((1,3),(16,6)):
            for col in (0,4,8,10):
                for prefix,cer,rx in (('H','cer',.13),('P','ngcer',.27)):
                    for suffix,value in ((cer+'.npeSum',7),('cal.etottracknorm',.9),('gtr.dp',dp),('react.x',rx),('react.y',-.09),('react.z',z),('rb.raster.fr_xbpm_tar',.07),('dc.x_fp',rx*10)):
                        values[prefix+'.'+suffix][0]=value
                rows.append((len(rows),foil,nd,col,dp));t.Fill()
    t.Write();f.Close();return rows


def make_masks(campaign,rows):
    for axis in ('Y','X'):
        folder=campaign/f'masks{axis}';folder.mkdir(parents=True)
        for foil in range(3):
            for nd,tag in ((3,'0_to_5'),(6,'15_to_20')):
                f=ROOT.TFile(str(folder/f'gmm_clean_test_foil{foil}_delta_{tag}_all{axis}.root'),'RECREATE')
                t=ROOT.TTree('GMMClean','Synthetic classification mask'); vals={}
                colname=axis.lower()+'scol'
                for name in ('entry','foil','ndel',colname,'gmm_keep'):
                    vals[name]=array('i',[0]);t.Branch(name,vals[name],name+'/I')
                score=array('d',[1.]);t.Branch('gmm_score',score,'gmm_score/D')
                for entry,nf,nidx,col,dp in rows:
                    if (nf,nidx)!=(foil,nd):continue
                    for name,value in zip(('entry','foil','ndel',colname,'gmm_keep'),(entry,nf,nidx,col,1)):vals[name][0]=value
                    t.Fill()
                t.Write();f.Close()


def run_macro(macro,expression,cwd,allow_error=False):
    result=subprocess.run(['root','-l','-b','-q','-e',f'#include {json.dumps(str(macro))}','-e',expression],cwd=cwd,text=True,capture_output=True,timeout=90)
    assert result.returncode==0 and (allow_error or 'error:' not in result.stderr.lower()),result.stdout+result.stderr
    return result


def check_pipeline():
    with tempfile.TemporaryDirectory(prefix='hallc-switch-') as name:
        temp=Path(name);(temp/'DATfiles').mkdir()
        metadata=temp/'DATfiles/list_of_optics_run.dat'
        metadata.write_text('99001,test,12.5,3,1,9,0.0\n-10,0,10\n-10,-8,-5,0,5,10,15,20,22\n')
        assert run_metadata(99001,metadata)['foils']==[-10,0,10]
        assert delta_index(run_metadata(99001,metadata)['edges'],15,20)==6
        replay=temp/'dual.root'; rows=make_replay(replay)
        for arm in ('HMS','SHMS'):
            campaign=temp/f'{arm}_same_target';make_masks(campaign,rows)
            args=[99001,-1,'test','',str(campaign/'masksY'),str(campaign/'masksX'),6.,.65,False,'',str(campaign/'06a_fit_ntuple'),str(replay),'test']
            expression='make_fit_ntuple_from_gmm('+','.join('kFALSE' if x is False else json.dumps(x) for x in args)+');'
            run_macro(REPO/'make_fit_ntuple_from_gmm.C',expression,temp)
            f=ROOT.TFile(str(campaign/'06a_fit_ntuple/root/Optics_99001_-1_fit_tree_gmm.root'))
            assert not f.IsZombie();t=f.Get('TFit');assert t
            expected=[r for r in rows if arm=='SHMS' or (r[3]<9 and r[4]<10)]
            assert t.GetEntries()==len(expected),(arm,t.GetEntries(),len(expected))
            assert str(f.Get('hallc_spectrometer').GetTitle())==arm
            cpp=ROOT.hallc.profileForName(arm)
            for i,(entry,foil,nd,col,dp) in enumerate(expected):
                t.GetEntry(i)
                assert (int(t.entry),int(t.foil),int(t.ndel),int(t.xscol),int(t.yscol))==(entry,foil,nd,col,col)
                assert t.ztarT==(-10,0,10)[foil]
                rx=.13 if arm=='HMS' else .27
                truth=ROOT.hallc.targetTruth(cpp,12.5,t.ztarT,cpp.xs(col),cpp.ys(col),dp,rx,-.09,.07)
                assert abs(t.reactxcalc-rx)<1e-14
                for name in ('xtar','ytar','xptar','yptar'):assert abs(getattr(t,name+'T')-getattr(truth,name))<1e-13
            f.Close()
            # A one-term toy fit exercises variable foil/slice allocation and output naming.
            # Its matrix is only a software fixture, never a calibration candidate.
            config=campaign/'config';config.mkdir()
            seed=config/'toy.dat';seed.write_text('! synthetic constant term\n ---------------------------------------------\n 0 0 0 0 00000\n ---------------------------------------------\n')
            groups=config/'rungroups_test_inputs.tsv';groups.write_text('rungroup\toptics_id\ntest\t99001\n')
            svdargs=['toy',-1,100,-1,str(campaign/'06a_fit_ntuple/root'),str(campaign/'06b_svd_fit'),str(seed),str(groups),str(metadata)]
            result=run_macro(REPO/'fit_opt_matrix_gmm.C','fit_opt_matrix_gmm('+','.join(json.dumps(x) for x in svdargs)+');',temp,allow_error=True)
            assert f'number to fit = {len(expected)} ' in result.stdout,result.stdout+result.stderr
            assert (campaign/f'06b_svd_fit/matrices/nps_{arm.lower()}_newfit_toy.dat').is_file()
            plot_metadata=temp/(arm+'_plot_metadata.dat')
            edges='-10,-8,-5,0,5,10' if arm=='HMS' else '-10,-8,-5,0,5,10,15,20,22'
            intervals=len(edges.split(','))-1
            plot_metadata.write_text(f'99001,test,12.5,3,1,{intervals+1},0.0\n-10,0,10\n{edges}\n99002,following,15,1,1,3,0.0\n0\n-10,0,10\n')
            for quantity in ('ytar','yptar','xptar'):
                plotargs=[99001,str(campaign/'06a_fit_ntuple/root/Optics_99001_-1_fit_tree_gmm.root'),str(plot_metadata),str(campaign/'residuals'),'toy']
                result=run_macro(REPO/f'diagnostics/fit/plot_{quantity}_residuals.C',f'plot_{quantity}_residuals('+','.join(json.dumps(x) for x in plotargs)+');',temp,allow_error=True)
                assert f'{intervals+1} delta boundaries; {intervals} intervals' in result.stdout,result.stdout+result.stderr
                plotted=ROOT.TFile(str(campaign/f'residuals/toy_{quantity}_residuals.root'))
                keys=[key.GetName() for key in plotted.GetListOfKeys()]
                used=[int(m.group(1)) for key in keys if (m:=re.search(r'_DelCut_(\d+)',key))]
                assert max(used)==intervals-1 and min(used)==0,keys
                plotted.Close()
            if arm=='HMS':
                # Compare all TFit columns against the unchanged baseline macro.
                baseline=temp/'baseline.C'
                baseline.write_text(subprocess.check_output(['git','show','057b908:make_fit_ntuple_from_gmm.C'],cwd=REPO,text=True))
                oldout=campaign/'baseline';args[10]=str(oldout)
                run_macro(baseline,'make_fit_ntuple_from_gmm('+','.join('kFALSE' if x is False else json.dumps(x) for x in args)+');',temp)
                fa=ROOT.TFile(str(campaign/'06a_fit_ntuple/root/Optics_99001_-1_fit_tree_gmm.root'));fb=ROOT.TFile(str(oldout/'root/Optics_99001_-1_fit_tree_gmm.root'))
                ta=fa.Get('TFit');tb=fb.Get('TFit');assert ta.GetEntries()==tb.GetEntries()
                for i in range(ta.GetEntries()):
                    ta.GetEntry(i);tb.GetEntry(i)
                    for branch in ta.GetListOfBranches():
                        key=branch.GetName();assert math.isclose(getattr(ta,key),getattr(tb,key),rel_tol=1e-13,abs_tol=1e-13),key
                fa.Close();fb.Close()
        # An explicit shifted sieve must not produce centered-geometry output.
        metadata.write_text(metadata.read_text().replace(',3,1,9',',3,2,9'))
        shifted=temp/'SHMS_shifted';make_masks(shifted,rows)
        args=[99001,-1,'test','',str(shifted/'masksY'),str(shifted/'masksX'),6.,.65,False,'',str(shifted/'06a_fit_ntuple'),str(replay),'test']
        result=run_macro(REPO/'make_fit_ntuple_from_gmm.C','make_fit_ntuple_from_gmm('+','.join('kFALSE' if x is False else json.dumps(x) for x in args)+');',temp,allow_error=True)
        assert 'centered sieve' in result.stderr and not (shifted/'06a_fit_ntuple/root').exists()

if __name__=='__main__':
    check_profiles();check_pipeline()
    print('PASS: shared profile parity; independent -10/0/+10 target on both arms; geometry closure; mask-to-TFit integration including SHMS column10 and delta16%; all HMS TFit columns match baseline; toy SVD accounting/output routing; residual interval counts; shifted sieve rejected')
