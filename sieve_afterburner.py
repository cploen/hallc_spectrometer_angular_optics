#!/usr/bin/env python3
"""Compare HMS matrices in sieve space using full replay events and fixed vertices."""
import argparse
import json
import re
import shutil
from pathlib import Path
import numpy as np
import uproot
from threadpoolctl import threadpool_limits
from core_sample import PROJECT,digest,read_tsv,write_tsv
from compare_matrices import read_matrix
from spectrometer_config import from_campaign,run_metadata

FIELDS={'xfp':'dc.x_fp','xpfp':'dc.xp_fp','yfp':'dc.y_fp','ypfp':'dc.yp_fp',
        'ry':'react.y','rz':'react.z','delta':'gtr.dp','cer':'cer.npeSum',
        'cal':'cal.etottracknorm','xs':'extcor.xsieve','ys':'extcor.ysieve'}
DEFAULTS=dict(angle=12.490,source='beam10',delta=[-10,10],cer=2.,cal=.65,
              foil_width=2.,replay_offsets=[0,0,0],gmm_offsets=[0,0,0],beam_offsets=[0,0,0],
              xmis=None,closure=.05)


def coefficients(a,rows):
    """Collapse four focal-plane variables once; retain a polynomial in xtar."""
    variables=[a['xfp']/100,a['xpfp'],a['yfp']/100,a['ypfp']]
    powers=[{k:v**k for k in {t[j] for c,t in rows}} for j,v in enumerate(variables)]
    result=np.zeros((max(t[4] for c,t in rows)+1,len(a['xfp']),3))
    for c,t in rows:
        if not np.any(c[:3]):continue
        term=np.ones(len(a['xfp']))
        for j in range(4):term*=powers[j][t[j]]
        result[t[4]]+=term[:,None]*c[:3]
    return result


def evaluate(poly,x,offset):
    value=np.zeros((len(x),3));xm=x[:,None]/100
    for row in poly[::-1]:value=value*xm+row
    return value+np.asarray(offset)/[1000,100,1000]


def reconstruct(a,rows,angle,offset,xmis=None,return_target=False):
    """HCANA ExtTarCor convention: saved vertex; >=1, <=5 updates, 2 mrad stop."""
    theta=np.deg2rad(angle)
    if xmis is None:
        degree=min(abs(angle),50);xmis=.1*(2.37-.086*degree+.0012*degree**2)
    x0=-a['ry']-xmis;poly=coefficients(a,rows)
    result=evaluate(poly,x0,offset);xtar=x0.copy();active=np.ones(len(x0),bool)
    iterations=np.zeros(len(x0),int)
    for iteration in range(5):
        proposal=x0-result[:,0]*a['rz']*np.cos(theta)
        candidate=evaluate(poly,proposal,offset)
        change=1000*np.abs(candidate[:,0]-result[:,0])
        xtar[active]=proposal[active];result[active]=candidate[active];iterations[active]+=1
        active &= change>2
        if not active.any():break
    xy=np.column_stack((xtar+168*result[:,0],100*result[:,1]+168*result[:,2]))
    if not np.isfinite(xy).all():raise ValueError('Nonfinite afterburner result; no events silently discarded')
    if return_target:
        return xy,iterations,active,result,xtar
    return xy,iterations,active


def selected(a,foils,cfg):
    good=(a['cer']>cfg['cer'])&(a['cal']>cfg['cal'])&(a['delta']>cfg['delta'][0])&(a['delta']<cfg['delta'][1])
    distance=np.abs(a['rz'][:,None]-np.asarray(foils)[None,:])
    index=distance.argmin(axis=1)
    good &= distance[np.arange(len(index)),index]<cfg['foil_width']
    return good,index


def prepare(campaign,cfg,gmm=None,beam=None):
    if from_campaign(campaign).name!='HMS':raise ValueError('This afterburner currently supports HMS only')
    tables=list((campaign/'config').glob('rungroups_*_inputs.tsv'))
    if len(tables)!=1:raise ValueError('Expected one campaign rungroup table')
    groups=[r for r in read_tsv(tables[0]) if abs(float(r['hms_angle_deg'])-cfg['angle'])<1e-6]
    if not groups:raise ValueError('No rungroups at the requested angle')
    policy=json.loads((campaign/'config/comparison.json').read_text())
    paths=dict(replay=campaign/'config/oldfit.dat',gmm=Path(gmm) if gmm else campaign/policy['gmm'],
               beam=Path(beam) if beam else campaign/'06e_beam_search'/cfg['source']/'matrices/beam.dat')
    rows={}
    for model,path in paths.items():
        excluded=policy.get('exclude_rows',{}).get(model,{}) if not (model=='gmm' and gmm) else {}
        if excluded and digest(path)!=policy.get(model+'_sha256'):
            raise ValueError('GMM row-exclusion checksum mismatch')
        rows[model]=read_matrix(path,excluded.get('rows',[]))
    for r in groups:
        r['foils']=run_metadata(int(r['optics_id']),PROJECT/'DATfiles/list_of_optics_run.dat')['foils']
        if len(r['foils'])>1 and 2*cfg['foil_width']>min(np.diff(sorted(r['foils']))):
            raise ValueError('Foil windows overlap; choose a smaller --foil-width')
        with uproot.open(r['rootfile']) as f:
            tree=f['T'];missing=set('H.'+v for v in FIELDS.values())-set(tree.keys())
            if missing:raise ValueError(f"{r['rungroup']}: missing branches {sorted(missing)}")
            r['entries']=tree.num_entries
    return groups,paths,rows,tables[0]


def run(campaign,name,cfg,threads=8,check=False,gmm=None,beam=None):
    out=campaign/'07_diagnostics'/name
    if not check and out.exists():raise FileExistsError(f'Output exists: {out}')
    groups,paths,matrices,table=prepare(campaign,cfg,gmm,beam)
    yedges=np.linspace(-7,7,281);xedges=np.linspace(-13,13,521)
    hist={};counts={};closure=[]
    manifest=dict(schema='sieve_afterburner_v1',status='running',settings=cfg,threads=threads,
                  inputs={str(p):digest(p) for p in [table,*paths.values(),campaign/'config/comparison.json',PROJECT/'DATfiles/list_of_optics_run.dat']},
                  replays=[dict(path=r['rootfile'],entries=r['entries'],size=Path(r['rootfile']).stat().st_size,
                                mtime_ns=Path(r['rootfile']).stat().st_mtime_ns) for r in groups],
                  reconstruction='Saved golden focal-plane tracks and react.y/react.z; HCANA extended-target iteration; fixed original delta and foil windows. No retracking or new reaction-vertex solve.')
    def save():(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if not check:
        for d in ['plots','tsv','matrices','code']:(out/d).mkdir(parents=True,exist_ok=True)
        save()
        for m,p in paths.items():shutil.copy2(p,out/f'matrices/{m}.dat')
        for file in ['sieve_afterburner.py','run_sieve.sh','compare_matrices.py','core_sample.py','spectrometer_config.py','spectrometer_profiles.def']:
            shutil.copy2(PROJECT/file,out/'code'/file)
    try:
        with threadpool_limits(limits=threads):
            for r in groups:
                rg=r['rungroup'];total=0;invalid=0;ncheck=0;sums=np.zeros(2);squares=np.zeros(2);maxdiff=np.zeros(2);unfinished={m:0 for m in matrices}
                for z in r['foils']:
                    key=(rg,z);counts[key]=0
                    for m in ('gmm','beam'):hist[m,rg,z]=np.zeros((280,520),dtype=np.int64)
                print(f"Read {rg}: {r['entries']:,} replay entries"+(' (check: first 50,000)' if check else ''),flush=True)
                with uproot.open(r['rootfile']) as f:
                    for raw in f['T'].iterate(['H.'+v for v in FIELDS.values()],step_size=50000,entry_stop=min(r['entries'],50000) if check else None,library='np'):
                        a={k:np.asarray(raw['H.'+v]) for k,v in FIELDS.items()}
                        if any(v.ndim!=1 or v.dtype.kind not in 'fiu' for v in a.values()):raise ValueError('Expected scalar numeric replay branches')
                        total+=len(a['rz']);finite=np.logical_and.reduce([np.isfinite(v)&(np.abs(v)<1e10) for v in a.values()])
                        invalid+=int(np.count_nonzero(~finite));a={k:v[finite] for k,v in a.items()}
                        keep,foil=selected(a,r['foils'],cfg);a={k:v[keep] for k,v in a.items()};foil=foil[keep]
                        if not len(foil):continue
                        for j,z in enumerate(r['foils']):counts[rg,z]+=int(np.count_nonzero(foil==j))
                        for m,rows in matrices.items():
                            xy,it,active=reconstruct(a,rows,float(r['hms_angle_deg']),cfg[m+'_offsets'],cfg['xmis'])
                            unfinished[m]+=int(np.count_nonzero(active))
                            if m=='replay':
                                diff=xy-np.column_stack((a['xs'],a['ys']));ncheck+=len(diff)
                                sums+=diff.sum(axis=0);squares+=(diff*diff).sum(axis=0);maxdiff=np.maximum(maxdiff,np.abs(diff).max(axis=0))
                            elif not check:
                                for j,z in enumerate(r['foils']):
                                    v=xy[foil==j];h,_,_=np.histogram2d(v[:,1],v[:,0],bins=(yedges,xedges));hist[m,rg,z]+=h.astype(np.int64)
                        print(f'  {total:,}/{r["entries"]:,} read; {ncheck:,} selected',flush=True)
                if not ncheck:raise ValueError(f'{rg}: no events pass cuts in inspected input')
                rms=np.sqrt(squares/ncheck)
                for j,axis in enumerate(('xsieve','ysieve')):
                    closure.append(dict(rungroup=rg,axis=axis,n=ncheck,bias_cm=sums[j]/ncheck,rms_cm=rms[j],max_cm=maxdiff[j],passed=bool(rms[j]<=cfg['closure']),invalid_entries=invalid,**{m+'_iteration_limit':v for m,v in unfinished.items()}))
                print(f"  Replay closure RMS: xs={rms[0]:.6g}, ys={rms[1]:.6g} cm; "+('PASS' if max(rms)<=cfg['closure'] else 'WARNING: inspect replay settings'),flush=True)
        if check:
            print('Check complete; no plots or files written. Closure uses only the first 50,000 entries per replay.');return
        manifest['closure_passed']=all(r['passed'] for r in closure)
        write_tsv(out/'tsv/closure.tsv',closure)
        report=[]
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        vmax=max(2,max(h.max() for h in hist.values()))
        for (m,rg,z),h in hist.items():
            label=('m' if z<0 else 'p')+f'{abs(z):g}' if z else '0'
            filename=f'{m}_{rg}_foil{label}'
            n=counts[rg,z];outside=n-int(h.sum())
            fig,ax=plt.subplots(figsize=(7,8),layout='constrained')
            im=ax.pcolormesh(yedges,xedges,np.ma.masked_equal(h.T,0),norm=LogNorm(1,vmax),cmap='viridis',rasterized=True)
            ax.set(xlabel='Y sieve (cm)',ylabel='X sieve (cm)',aspect='equal',xlim=(-7,7),ylim=(-13,13))
            ax.set_title(f'{m.upper()} | foil {z:+g} cm | θ={cfg["angle"]:g}°\n{cfg["delta"][0]:g}<δ<{cfg["delta"][1]:g}%; N={n:,}; outside frame={outside:,}')
            fig.colorbar(im,ax=ax,label='Events per bin (log scale)')
            fig.suptitle('Fixed-vertex afterburner'+(' — replay check needs review' if not manifest['closure_passed'] else ''),fontsize=11)
            fig.savefig(out/f'plots/{filename}.png',dpi=160);plt.close(fig)
            report.append(dict(model=m,rungroup=rg,zfoil=z,n=n,outside_frame=outside,file=filename+'.png'))
            print(f'Wrote {filename}.png',flush=True)
        write_tsv(out/'tsv/counts.tsv',report)
        np.savez_compressed(out/'histograms.npz',xedges=xedges,yedges=yedges,**{f'{m}_{rg}_{z:g}':h for (m,rg,z),h in hist.items()})
        manifest['status']='complete';manifest['outputs']={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and p.name!='manifest.json'};save()
        print(f'Complete: {out}',flush=True)
    except BaseException as exc:
        if not check:manifest.update(status='incomplete',error=str(exc));save()
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path);p.add_argument('name',nargs='?',default='sieve')
    p.add_argument('--angle',type=float);p.add_argument('--source');p.add_argument('--delta',nargs=2,type=float)
    p.add_argument('--cer',type=float);p.add_argument('--cal',type=float);p.add_argument('--foil-width',type=float)
    p.add_argument('--gmm',type=Path);p.add_argument('--beam',type=Path);p.add_argument('--threads',type=int,default=8);p.add_argument('--check',action='store_true')
    a=p.parse_args();campaign=a.campaign.resolve();cfg=DEFAULTS.copy()
    try:
        config=campaign/'config/sieve.json'
        if config.exists():
            supplied=json.loads(config.read_text())
            if set(supplied)-set(cfg):raise ValueError('Unknown sieve configuration key')
            cfg.update(supplied)
        for k in ('angle','source','delta','cer','cal','foil_width'):
            if getattr(a,k) is not None:cfg[k]=getattr(a,k)
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',a.name) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',cfg['source']):raise ValueError('Invalid output/source name')
        if a.threads<1 or cfg['foil_width']<=0 or cfg['closure']<=0 or len(cfg['delta'])!=2 or cfg['delta'][0]>=cfg['delta'][1]:raise ValueError('Invalid threads, window or threshold')
        for m in ('replay','gmm','beam'):
            v=cfg[m+'_offsets']
            if len(v)!=3 or not np.isfinite(v).all():raise ValueError('Offsets must be [xptar mrad, ytar cm, yptar mrad]')
        run(campaign,a.name,cfg,a.threads,a.check,a.gmm,a.beam)
    except (OSError,ValueError,KeyError,RuntimeError) as exc:p.exit(1,f'ERROR: {exc}\n')


if __name__=='__main__':main()
