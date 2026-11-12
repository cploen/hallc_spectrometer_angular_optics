#!/usr/bin/env python3
"""Rescore saved xtar candidates; freeze choices before reserved evaluation.

No new search, selection changes, transport edits, or replay installation.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys
import numpy as np
from xtar_correction import (EDGES,TARGETS,UNITS,OLD_HASH,REFIT_HASH,
    Reconstruction,Metrics,subset,features,truth,scaled_fit,attach_vertex,validate_geometry)

CASES={'cell2_overall2':(.02,.02),'cell5_overall2':(.05,.02),'cell10_overall5':(.10,.05)}


class Moments:
    """Fixed event groups, reduced with bincount instead of repeated mask scans."""
    def __init__(self,a,minimum=100):
        self.minimum=minimum
        eligible=np.isin(a['quality'],[1,2])&~a['excluded']
        nd=np.searchsorted(EDGES,a['delta'],side='right')-1
        cells=np.array([f'{z:g}:{d}' for z,d in zip(a['ztarT'],nd)])
        names,inv=np.unique(cells,return_inverse=True);nc=len(names)
        idx=[];codes=[];labels=['overall','core','shoulder']
        labels += [f'{pop}|{cell}' for pop in ('overall','core','shoulder') for cell in names]
        for j,m in enumerate((eligible,eligible&(a['quality']==2),eligible&(a['quality']==1))):
            ii=np.flatnonzero(m)
            idx.extend((ii,ii));codes.extend((np.full(len(ii),j),3+j*nc+inv[ii]))
        self.idx=np.concatenate(idx);self.codes=np.concatenate(codes)
        self.labels=labels;self.count=np.bincount(self.codes,minlength=len(labels))
        self.good=self.count>=minimum;self.cell=np.arange(len(labels))>=3
        self.score_cells=np.arange(3,3+nc)

    def stats(self,r):
        sums=np.column_stack([np.bincount(self.codes,weights=r[self.idx,j],minlength=len(self.labels)) for j in range(4)])
        squares=np.column_stack([np.bincount(self.codes,weights=r[self.idx,j]**2,minlength=len(self.labels)) for j in range(4)])
        mean=sums/np.maximum(self.count[:,None],1)
        sd=np.sqrt(np.maximum(0,squares/np.maximum(self.count[:,None],1)-mean**2))
        return mean,sd

    def score(self,stats):
        i=self.score_cells
        return float(np.sum(self.count[i]*stats[1][i,3]**2)/np.sum(self.count[i]))

    def changes(self,b,c):
        ratio=c[1][:,:3]/np.maximum(b[1][:,:3],1e-12)-1
        bias=(np.abs(c[0][:,:3])-np.abs(b[0][:,:3]))/np.maximum(b[1][:,:3],1e-12)
        def worst(v,mask):
            return np.max(v[mask],axis=0) if mask.any() else np.full(3,-1.)
        return dict(overall=worst(ratio,self.good&~self.cell),
                    cell=worst(ratio,self.good&self.cell),bias=worst(bias,self.good))


def permitted(r,cell,overall):
    return bool(r['physical'] and max(r['global'])<=overall+1e-10 and
                max(r['cell'])<=cell+1e-10 and max(r['bias'])<=.02+1e-10)


def read_table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def write_table(p,rows):
    if not rows:return
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


def scalar_record(r):
    v={k:r[k] for k in ('id','name','kind','alpha','l1','amplitude','terms','score','mse','condition','physical','correction_max_cm')}
    v['mask']=','.join(map(str,r['mask']))
    for key in ('global','cell','bias'):
        for j,t in enumerate(TARGETS[:3]):v[f'{key}_{t}']=float(r[key][j])
    for case,(c,o) in CASES.items():v[case]=permitted(r,c,o)
    return v


def run(args):
    sys.path.insert(0,str(args.checkout))
    from fit_elastic import load_sample,matrix_rows,checked
    from core_sample import read_tsv,digest
    from threadpoolctl import threadpool_limits
    original=args.experiment/'results';out=args.out
    out.mkdir(parents=True,exist_ok=False);(out/'plots').mkdir();(out/'models').mkdir()
    shutil.copy2(__file__,out/'xtar_tolerance_scan.py')
    old_manifest=json.loads((original/'manifest.json').read_text())
    if old_manifest['status']!='complete':raise ValueError('Original run is incomplete')
    campaign=Path(old_manifest['campaign']);hashes={}
    refit=original/'matrices/refit.dat';old=original/'matrices/oldfit.dat'
    checked(refit,REFIT_HASH,hashes);checked(old,OLD_HASH,hashes)
    source=original/'search_candidates.tsv';hashes[str(source)]=digest(source)
    # Freeze the exact original implementation as well as its matrix/data IDs.
    frozen=original/'code/xtar_correction.py'
    checked(Path(__file__).with_name('xtar_correction.py'),digest(frozen),hashes)
    manifest=dict(status='running',cases=CASES,bias_limit=.02,minimum=100,
        mode='Rescoring the saved candidate set; no new beam exploration.',
        source_candidates=str(source),constants='All stored constants retained; zero external additions.',
        heldout='Previously inspected reserved sample is development evidence.',inputs=hashes)
    def save():(out/'manifest.json').write_text(json.dumps(manifest,indent=2,default=str)+'\n')
    save()
    try:
      with threadpool_limits(limits=1):
        a,hh=load_sample(campaign,'inner20','fit',xtar_coordinate='xtarT');hashes.update(hh)
        a=attach_vertex(campaign,'inner20','fit',a,read_tsv);validate_geometry(a)
        x,labels=features(a);tt=truth(a);rows=matrix_rows(refit);reco=Reconstruction(a,rows)
        residual0=reco.base-tt;y=a['xtarT']-reco.base[:,3]/10
        foldrows=read_table(original/'folds.tsv');lookup={(r['rungroup'],int(r['entry'])):int(r['fold']) for r in foldrows}
        folds=np.array([lookup[str(r),int(e)] for r,e in zip(a['rungroup'],a['entry'])])
        if len(lookup)!=len(x):raise ValueError('Fold ID mismatch')
        metric=Moments(a);bstat=metric.stats(residual0);baseline=metric.score(bstat)
        np.testing.assert_allclose(baseline,old_manifest['selection']['baseline_variance'],rtol=1e-10)
        # Verify the fast statistics against the established direct std path.
        slow=Metrics(a).measure(residual0)
        for i,label in enumerate(metric.labels):
            if label in slow:
                np.testing.assert_allclose(bstat[0][i],slow[label][1],rtol=1e-10,atol=1e-10)
                np.testing.assert_allclose(bstat[1][i],slow[label][2],rtol=1e-9,atol=1e-10)
        parts=[]
        for f in sorted(set(folds)):
            tr=folds!=f;va=folds==f;mm=Moments(subset(a,va))
            child=object.__new__(Reconstruction);child.anchor=reco.anchor[va];child.z=reco.z[va];child.poly=reco.poly[:,va]
            parts.append((x[tr],y[tr],x[va],va,mm,mm.stats(residual0[va]),child))
        candidates=read_table(source);records=[];fits={};max_score_difference=0.;eligibility_mismatches=0
        print(f'Rescoring {len(candidates)} frozen candidates on {len(x):,} training rows; three tolerance cases; one shared solve cache.',flush=True)
        for k,raw in enumerate(candidates):
            mask=tuple(int(s) for s in raw['mask'].split(',') if s)
            r=dict(id=k,name=raw['name'],kind=raw['kind'],alpha=float(raw['alpha']),
                   l1=float(raw['l1']),amplitude=float(raw['amplitude']),mask=mask)
            rr=np.empty_like(residual0);physical=True;term=0;condition=1.;gmax=0.
            change={key:np.full(3,-np.inf) for key in ('overall','cell','bias')}
            try:
                for f,(xt,yt,xv,va,mm,bs,child) in enumerate(parts):
                    key=(f,r['kind'],mask,r['alpha'],r['l1'])
                    if key not in fits:
                        try:fits[key]=scaled_fit(xt,yt,mask,r['kind'],r['alpha'],r['l1'])
                        except (ValueError,Warning,np.linalg.LinAlgError):fits[key]=None
                    if fits[key] is None:raise ValueError('Invalid scaled design')
                    coeff,info=fits[key];b=coeff*r['amplitude'];g=xv@b
                    term=max(term,int(np.count_nonzero(b)));condition=max(condition,info['condition'])
                    gmax=max(gmax,float(np.max(np.abs(g))))
                    p,limited=child.predict(g);rr[va]=p-tt[va]
                    physical &= bool(not np.any(limited&~reco.base_limited[va]))
                    cc=mm.changes(bs,mm.stats(rr[va]))
                    for key in change:change[key]=np.maximum(change[key],cc[key])
                st=metric.stats(rr);cc=metric.changes(bstat,st)
                for key in change:change[key]=np.maximum(change[key],cc[key])
                score=metric.score(st);mse=float(np.mean(rr[:,3]**2))
                physical &= bool(term<=6 and gmax<=.2 and mse<=float(np.mean(residual0[:,3]**2))*(1+1e-10))
                if float(raw['score'])<baseline*100:
                    difference=abs(score-float(raw['score']))/max(baseline,1e-20)
                    max_score_difference=max(max_score_difference,difference)
            except (ValueError,Warning,np.linalg.LinAlgError):
                physical=False;score=baseline*1e6;mse=1e6
                change={key:np.full(3,1e99) for key in change}
            r.update(physical=bool(physical),score=score,mse=mse,terms=term,condition=condition,
                     correction_max_cm=gmax,**{'global':change['overall'],'cell':change['cell'],'bias':change['bias']})
            if permitted(r,.02,.02)!=(raw['valid']=='True'):eligibility_mismatches+=1
            records.append(r)
            if (k+1)%200==0 or k+1==len(candidates):
                print(f'{k+1}/{len(candidates)} rescored; eligible '+str({c:sum(permitted(v,*lim) for v in records) for c,lim in CASES.items()}),flush=True)
                write_table(out/'candidate_scores.tsv',[scalar_record(v) for v in records])
        if max_score_difference>1e-8 or eligibility_mismatches:
            raise ValueError(f'Original comparison parity failed: score difference={max_score_difference}; eligibility mismatches={eligibility_mismatches}')
        print('Original 2%/2% eligibility and scores reproduced.',flush=True)
        choices={};coeffs={};selection=[]
        for case,limits in CASES.items():
            for family in ('EN','Beam'):
                pool=[r for r in records if r['name'] in (family,'unchanged') and permitted(r,*limits)]
                best=min(pool,key=lambda r:(r['score'],r['terms'],r['condition']))
                simple=min([r for r in pool if r['score']<=best['score']*1.005],key=lambda r:(r['terms'],r['score']))
                name=case+'_'+family;choices[name]=best
                b,info=scaled_fit(x,y,best['mask'],best['kind'],best['alpha'],best['l1']);b*=best['amplitude']
                coeffs[name]=b
                model=dict(features=labels,coefficients_cm=b.tolist(),intercept_cm=0.,spec=scalar_record(best),
                    final_fit=info,convention='Existing iterative xtar plus a nonaccumulating correction; frozen transport.',
                    status='Informational strongest eligible candidate, not installed')
                (out/'models'/f'{name}.json').write_text(json.dumps(model,indent=2)+'\n')
                selection.append(dict(case=case,family=family,id=best['id'],terms=best['terms'],
                    validation_xtar_width_gain_pct=100*(1-np.sqrt(best['score']/baseline)),
                    worst_global_loss_pct=100*max(best['global']),worst_cell_loss_pct=100*max(best['cell']),
                    simplicity_id=simple['id'],simplicity_terms=simple['terms'],eligible=len(pool)))
        write_table(out/'selection.tsv',selection)
        # Choices/coefficients are now frozen. Never retune after this point.
        manifest.update(selection=selection,max_score_parity_error=max_score_difference,
                        original_eligibility_mismatches=eligibility_mismatches,inputs=hashes);save()
        h,hh=load_sample(campaign,'inner20','holdout',xtar_coordinate='xtarT');hashes.update(hh)
        if set(zip(a['rungroup'],a['entry']))&set(zip(h['rungroup'],h['entry'])):raise ValueError('Training/reserved overlap')
        h=attach_vertex(campaign,'inner20','holdout',h,read_tsv);validate_geometry(h)
        hx,_=features(h);ht=truth(h);hr=Reconstruction(h,rows)
        with np.load(original/'reserved_predictions.npz') as saved:
            np.testing.assert_array_equal(saved['entry'],h['entry']);np.testing.assert_array_equal(saved['rungroup'],h['rungroup'])
            np.testing.assert_allclose(saved['refit'],hr.base,rtol=1e-12,atol=1e-10)
            old_prediction=saved['oldfit'].copy()
        hm=Moments(h);hbase=hr.base-ht;hb=hm.stats(hbase)
        pred={'oldfit':old_prediction,'refit':hr.base};evaluation=[];details=[]
        for name,b in coeffs.items():
            g=hx@b;p,lim=hr.predict(g);pred[name]=p;st=hm.stats(p-ht);loss=hm.changes(hb,st)
            case=name.rsplit('_',1)[0];cell_limit,global_limit=CASES[case]
            physical=bool(np.max(np.abs(g))<=.2 and not np.any(lim&~hr.base_limited) and np.count_nonzero(b)<=6)
            result={'physical':physical,'global':loss['overall'],'cell':loss['cell'],'bias':loss['bias']}
            row=dict(model=name,terms=int(np.count_nonzero(b)),
                xtar_width_gain_pct=100*(1-np.sqrt(hm.score(st)/hm.score(hb))),
                xtar_mse_gain_pct=100*(1-np.mean((p[:,3]-ht[:,3])**2)/np.mean(hbase[:,3]**2)),
                reserved_guards_pass=permitted(result,cell_limit,global_limit))
            for j,t in enumerate(TARGETS[:3]):
                row[t+'_overall_change_pct']=100*(st[1][0,j]/hb[1][0,j]-1)
                row[t+'_worst_population_change_pct']=100*loss['overall'][j]
                row[t+'_worst_cell_change_pct']=100*loss['cell'][j]
                row[t+'_worst_bias_increase_in_sigma']=loss['bias'][j]
            evaluation.append(row)
        for name,p in pred.items():
            mean,sd=hm.stats(p-ht)
            for i,label in enumerate(hm.labels):
                if not hm.count[i]:continue
                for j,t in enumerate(TARGETS):
                    details.append(dict(model=name,group=label,target=t,unit=UNITS[j],n=int(hm.count[i]),
                                        sd=float(sd[i,j]),bias=float(mean[i,j]),sparse=bool(hm.count[i]<100)))
        write_table(out/'reserved_summary.tsv',evaluation);write_table(out/'reserved_metrics.tsv',details)
        np.savez_compressed(out/'reserved_predictions.npz',truth=ht,entry=h['entry'],rungroup=h['rungroup'],
            delta=h['delta'],zfoil=h['ztarT'],quality=h['quality'],excluded=h['excluded'],**pred)
        plots(out,records,baseline,selection,evaluation,details,h)
        for p,sha in hashes.items():checked(Path(p),sha,{})
        manifest.update(status='complete',inputs=hashes,training=len(x),reserved=len(hx),
            conclusion='Informational tradeoff only; transport and adopted correction unchanged.',
            limitation='Saved candidate set was explored with the original strict gates. New relaxed-gate searches could find additional bases.')
        save();print('COMPLETE',out,flush=True)
        print((out/'selection.tsv').read_text(),flush=True);print((out/'reserved_summary.tsv').read_text(),flush=True)
    except BaseException as exc:
        manifest.update(status='incomplete',error=str(exc));save();raise


def plots(out,records,baseline,selection,evaluation,details,a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(15,10),layout='constrained')
    valid=[r for r in records if r['physical'] and max(r['bias'])<=.02+1e-10]
    gain=np.array([100*(1-np.sqrt(r['score']/baseline)) for r in valid])
    for ax,key,title in [(axes[0,0],'cell','Worst populated cell'),(axes[0,1],'global','Worst overall/population')]:
        ax.scatter([100*max(r[key]) for r in valid],gain,s=5,alpha=.25)
        ax.set(xlabel=title+' downstream SD increase (%)',ylabel='Validation xtar SD reduction (%)',title=title+' versus xtar gain')
        for limit in ((2,5,10) if key=='cell' else (2,5)):ax.axvline(limit,color='.5',ls=':',lw=.8)
        ax.set_xlim(-1,15);ax.grid(alpha=.2)
    for family,style in [('EN','o--'),('Beam','s-')]:
        ss=[r for r in selection if r['family']==family];ee=[next(e for e in evaluation if e['model']==r['case']+'_'+family) for r in ss]
        axes[1,0].plot(range(3),[r['validation_xtar_width_gain_pct'] for r in ss],style,label=family+' validation')
        axes[1,0].plot(range(3),[r['xtar_width_gain_pct'] for r in ee],style,alpha=.55,label=family+' reserved')
    axes[1,0].set(xticks=range(3),xticklabels=['2% cell / 2% overall','5% cell / 2% overall','10% cell / 5% overall'],
                   ylabel='xtar SD reduction (%)',title='Strongest eligible corrections; no simplicity suppression')
    axes[1,0].legend(fontsize=8);axes[1,0].grid(alpha=.2)
    beams=[r for r in evaluation if r['model'].endswith('_Beam')]
    for j,t in enumerate(TARGETS[:3]):axes[1,1].plot(range(3),[r[t+'_overall_change_pct'] for r in beams],'.-',label=t+' overall')
    axes[1,1].plot(range(3),[r['xtar_width_gain_pct']*-1 for r in beams],'s-',label='xtar within-cell')
    axes[1,1].set(xticks=range(3),xticklabels=['2 / 2','5 / 2','10 / 5'],ylabel='Reserved SD change (%); negative = narrower',title='Beam correction: reserved gain and downstream cost')
    axes[1,1].legend(fontsize=8);axes[1,1].grid(alpha=.2)
    fig.suptitle('HMS 2018 xtar correction — tolerance sensitivity\nSame saved candidates, run-group folds, constants and iterative reconstruction; bias guard unchanged',fontsize=14)
    fig.savefig(out/'plots/tolerance_summary.png',dpi=170);plt.close(fig)
    lookup={(r['model'],r['group'],r['target']):r for r in details}
    foils=sorted(np.unique(a['ztarT']));nc=len(EDGES)-1
    pairs=[('refit','oldfit')]+[(c+'_'+f,'refit') for c in CASES for f in ('EN','Beam')]
    for population in ('core','shoulder'):
      for metric in ('sd','bias'):
        fig,axs=plt.subplots(len(pairs),4,figsize=(32,25),layout='constrained')
        for row,(model,reference) in enumerate(pairs):
          for j,target in enumerate(TARGETS):
            ax=axs[row,j];v=np.full((len(foils),nc),np.nan);counts=np.zeros(v.shape,int)
            for i,z in enumerate(foils):
              for d in range(nc):
                group=f'{population}|{z:g}:{d}';r=lookup.get((model,group,target));b=lookup.get((reference,group,target))
                if r and b:
                    counts[i,d]=r['n'];v[i,d]=100*(r['sd']/max(b['sd'],1e-12)-1) if metric=='sd' else r['bias']-b['bias']
            finite=np.abs(v[np.isfinite(v)]);scale=20 if metric=='sd' else max(float(finite.max()) if len(finite) else 0,1e-5)
            im=ax.imshow(np.ma.masked_invalid(v),aspect='auto',cmap='RdBu_r',vmin=-scale,vmax=scale)
            ax.set(title=f'{target}: {model}\nminus {reference} | '+('SD change (%)' if metric=='sd' else f'signed bias change ({UNITS[j]})'),
                   yticks=range(len(foils)),yticklabels=foils,xticks=range(nc),
                   xticklabels=[f'[{l:g},{h:g})' for l,h in zip(EDGES[:-1],EDGES[1:])],ylabel='Foil (cm)')
            ax.tick_params(axis='x',rotation=40,labelsize=8)
            for i,d in np.ndindex(v.shape):
              if counts[i,d]:
                val=f'{v[i,d]:+.1f}%' if metric=='sd' else f'{v[i,d]:+.3g}'
                ax.text(d,i,val+f'\nN={counts[i,d]}'+('*' if counts[i,d]<100 else ''),ha='center',va='center',fontsize=6.1,
                        color='white' if abs(v[i,d])>scale*.5 else 'black')
            fig.colorbar(im,ax=ax,shrink=.65,pad=.01)
        fig.suptitle(f'{population}: strongest eligible EN and beam proposals by tolerance\n'+
            ('Blue = narrower, red = wider; color saturates ±20%, text uncapped. ' if metric=='sd' else 'Signed means are separate from width gains. ')+
            'Identical reserved events, stored constants retained, zero external additions. * N<100 uncertified.',fontsize=14)
        fig.savefig(out/'plots'/f'{population}_{metric}_changes.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--experiment',type=Path,required=True);p.add_argument('--checkout',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    run(p.parse_args())
