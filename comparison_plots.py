"""Frozen-matrix attribution tables, central residual zooms, and separate tail views."""
import json
import numpy as np
from core_sample import write_tsv
from elastic_diagnostics import groups,LABELS
from elastic_net import macro_mse

MODELS=('old','gmm','core','enet','beam')
NAMES=dict(old='Old replay + offsets',gmm='GMM fit + offsets',core='Full core SVD',
           enet='EN-selected SVD',beam='Beam')
TARGETS=('xptar','ytar','yptar','ztar')
POOLS=('protected_core','protected_noncore','surplus_core')


def stats(v):
    return dict(n=len(v),bias=float(v.mean()),spread=float(v.std()),rms=float(np.sqrt(np.mean(v*v))),
                median=float(np.median(np.abs(v))),p90=float(np.percentile(np.abs(v),90)))


def report(out,a,pool,rr,known,seen,history):
    cells=np.array([f'{z}:{d}' for z,d in zip(a['ztarT'],a['ndel'])])
    summary,details,overlap,changes=[],[],[],[]
    for p in np.unique(pool):
        mask=pool==p
        overlap.append(dict(pool=p,n=int(sum(mask)),gmm_seen=int(sum(mask & known & seen)),
                            gmm_unseen=int(sum(mask & known & ~seen)),gmm_unknown=int(sum(mask & ~known)),
                            old_seen='unknown'))
        if p=='blocked':continue
        for subset,sel in (('all',mask),('outside_gmm',mask & known & ~seen)):
            if not sel.any():continue
            for j,t in enumerate(TARGETS):
                for m in MODELS:
                    v=rr[m][sel,j]
                    summary.append(dict(pool=p,subset=subset,model=m,target=t,
                                        cell_mse=float(macro_mse(v,cells[sel])),**stats(v)))
    for g,sel in groups(a,pool):
        for j,t in enumerate(TARGETS):
            for m in MODELS:
                details.append(dict(g,model=m,target=t,qa_low=int(sum(sel))<10,**stats(rr[m][sel,j])))
    for r in summary:
        index=MODELS.index(r['model'])
        for label,base in [('starting','old')]+([('previous',MODELS[index-1])] if index else []):
            b=next(v for v in summary if (v['pool'],v['subset'],v['target'],v['model'])==
                   (r['pool'],r['subset'],r['target'],base))
            for metric in ('bias','spread','rms','median','p90','cell_mse'):
                changes.append(dict(pool=r['pool'],subset=r['subset'],target=r['target'],model=r['model'],
                                    reference=base,comparison=label,metric=metric,value=r[metric],
                                    reference_value=b[metric],change=r[metric]-b[metric],
                                    percent_change=100*(r[metric]/b[metric]-1) if b[metric] and metric!='bias' else ''))
    for f,rows in [('summary',summary),('residuals',details),('overlap',overlap),('changes',changes)]:
        write_tsv(out/f'tsv/{f}.tsv',rows)
    np.savez_compressed(out/'residuals.npz',rungroup=a['rungroup'],entry=a['entry'],pool=pool,
                        zfoil=a['ztarT'],ndel=a['ndel'],xscol=a['xscol'],yscol=a['yscol'],
                        gmm_known=known,gmm_seen=seen,targets=np.array(TARGETS),**rr)
    plots(out,a,pool,rr,details)
    text='# Frozen matrix comparison\n\n'
    policy=json.loads((out/'manifest.json').read_text()).get('policy',{}) if (out/'manifest.json').exists() else {}
    for m,item in policy.get('exclude_rows',{}).items():
        text+=f"**Matrix-file limitation ({m}):** {item['reason']} Exact excluded row(s) and the original file checksum are recorded in manifest.json.\n\n"
    text+='Order: old replay plus external offsets → GMM fit plus external offsets → full-basis core SVD → elastic-net-selected SVD refit → beam.\n\n'
    for p in POOLS:
        rows=[r for r in summary if r['pool']==p and r['subset']=='all']
        if not rows:continue
        text+=f'## {p}\n\n| Matrix | xptar RMS (mrad) | ytar RMS (cm) | yptar RMS (mrad) | ztar RMS (cm) |\n|---|---:|---:|---:|---:|\n'
        for m in MODELS:
            text+='| '+NAMES[m]+' | '+' | '.join(f"{next(r for r in rows if r['model']==m and r['target']==t)['rms']:.5g}" for t in TARGETS)+' |\n'
        text+='\n'
    text+=('## What these comparisons establish\n\n'
           'Every matrix is evaluated on the same saved focal-plane coordinates and xtar input, with its '
           'entire angular polynomial and explicit external corrections. No coefficient or offset is fitted '
           'on these evaluation events. Matrices and offset sources are frozen in this output. This is a '
           'common-input reconstruction comparison, not a full iterative HCANA replay with new event selection.\n\n'
           'Old → GMM quantifies the historical change. GMM → core combines event selection, balancing and '
           'solver changes; it does not isolate tighter cores alone. Core → EN-selected SVD → beam holds '
           'the training events, fixed terms and direct-X solver constant and measures the basis-change '
           'effect. EN-selected SVD is the all-training refit of the chosen compact seed basis, not the '
           'penalized elastic-net prediction. Its identity is the beam source matrices/start.dat.\n\n'
           'External offsets must come from historical settings, not from centering these residuals. '
           'tsv/offsets.tsv lists additions in mrad/cm. New core, EN and beam fits already contain freely '
           'fitted intercepts, so no historical correction is added to those by default. Bias, spread and '
           'RMS are reported separately: RMS² = bias² + spread². A lower RMS caused only by removing a '
           'constant offset is not an improvement in mean-subtracted width. Spread is not automatically '
           'an experimental fitted Gaussian resolution. Local bias/spread must be examined as well as '
           'pooled values, where opposing cell biases can cancel.\n\n'
           'Protected cores and protected noncore events are separate. Surplus cores are additional '
           'development coverage, not an independent final test. Unsupported labels are tabulated separately; '
           'blocked positions have counts only. Shoulder residuals depend on the validity of inherited labels.\n\n'
           'tsv/overlap.tsv identifies events used by the historical GMM fit, based on reconstructed '
           'membership rather than an original solver ID log. tsv/summary.tsv and changes.tsv also report '
           'the common subset outside that reconstructed GMM training membership. Unknown membership is '
           'not counted as unseen. Training membership for the old replay matrix remains unknown. Thus '
           'the all-events historical comparison is descriptive; new-fit held-out performance must not '
           'be confused with an independently held-out test of every historical matrix.\n\n'
           'plots/*_center.png zoom to ± the beam P90 absolute residual, with linear vertical axes and '
           'identical bins for every matrix. Densities are normalized to the entire pool, not just the '
           'visible events; percentages in the panel list how much each curve shows. No residual is '
           'recentered for the plot. plots/*_tails.png shows the fraction at or beyond each absolute '
           'residual over the full range. plots/*_foil_delta.png shows absolute RMS for all five matrices, '
           'with shared scales and N. Full bias, spread, median and P90 down to setting/hole are in '
           'tsv/residuals.tsv; qa_low flags N<10 and never changes the sample.\n\n'
           'tsv/changes.tsv reports both previous-stage and starting-matrix differences; negative width '
           'changes mean improvement. Mean cell MSE gives equal weight to populated foil/delta cells; '
           'other pooled statistics weight events equally. tsv/training_conditioning.tsv preserves the '
           'source training-design condition numbers for the modern fits only. These are not condition '
           'numbers of the held-out residuals, and historical effective rank cannot be inferred just '
           'by counting nonzero polynomial coefficients.\n\n'
           'Once these protected results guide additional tuning, subsequent performance on them is '
           'development evidence. This report makes no automatic claim of improvement.\n\n'
           'Historical membership method: '+history.get('method','unspecified')+'\n')
    (out/'MATRIX_COMPARISON.md').write_text(text)


def plots(out,a,pool,rr,details):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    colors=('0.35','C1','C0','C4','C2');styles=(':','--','-','-.','-')
    for p in POOLS:
        sel=pool==p
        if not sel.any():continue
        center,axes=plt.subplots(1,4,figsize=(18,4.5),constrained_layout=True)
        tails,tax=plt.subplots(1,4,figsize=(18,4.5),constrained_layout=True)
        for j,(ax,tx) in enumerate(zip(axes,tax)):
            limit=max(float(np.percentile(np.abs(rr['beam'][sel,j]),90)),1e-9)
            bins=np.linspace(-limit,limit,81);mid=(bins[:-1]+bins[1:])/2
            fractions=[]
            for m,color,style in zip(MODELS,colors,styles):
                v=rr[m][sel,j];count,_=np.histogram(v,bins)
                ax.plot(mid,count/(len(v)*np.diff(bins)),color=color,ls=style,label=NAMES[m])
                fractions.append(f'{NAMES[m]}: {100*np.mean(np.abs(v)<=limit):.1f}% visible')
                # Exact empirical survival at sampled unique values, including tied residuals.
                sorted_v=np.sort(np.abs(v));absolute,first=np.unique(sorted_v,return_index=True)
                take=np.unique(np.linspace(0,len(absolute)-1,min(700,len(absolute))).astype(int))
                positive=absolute[take]>0;take=take[positive]
                tx.plot(absolute[take],(len(v)-first[take])/len(v),color=color,ls=style,label=NAMES[m])
            ax.set(xlabel=LABELS[j],ylabel='Density per physical unit',title='Central signed residual',xlim=(-limit,limit))
            ax.text(.02,.98,'\n'.join(fractions),va='top',transform=ax.transAxes,fontsize=6)
            tx.set(xlabel='Absolute '+LABELS[j],ylabel='Fraction at or beyond',xscale='log',yscale='log',title='Full-range tails')
        center.legend(*axes[-1].get_legend_handles_labels(), loc='outside lower center', ncol=5, fontsize=8)
        tax[-1].legend(fontsize=7)
        center.suptitle(f'{p}: N={sum(sel):,}; central window = ±beam P90; no recentering')
        tails.suptitle(f'{p}: all events included in tail probabilities')
        center.savefig(out/f'plots/{p}_center.png',dpi=140);plt.close(center)
        tails.savefig(out/f'plots/{p}_tails.png',dpi=140);plt.close(tails)
        local=[r for r in details if r['pool']==p and r['level']=='foil_delta']
        foils=sorted({r['zfoil'] for r in local});deltas=sorted({r['ndel'] for r in local})
        lookup={(r['model'],r['target'],r['zfoil'],r['ndel']):r for r in local}
        fig,grid=plt.subplots(5,4,figsize=(17,15),constrained_layout=True)
        for j,t in enumerate(TARGETS):
            norm=Normalize(0,max((r['rms'] for r in local if r['target']==t),default=1) or 1)
            for i,m in enumerate(MODELS):
                ax=grid[i,j];values=np.full((len(foils),len(deltas)),np.nan)
                for zi,z in enumerate(foils):
                    for di,d in enumerate(deltas):
                        r=lookup.get((m,t,z,d))
                        if r:
                            values[zi,di]=r['rms']
                            ax.text(di,zi,f"{r['rms']:.3f}\nN={r['n']}",ha='center',va='center',fontsize=5,
                                    color='white' if norm(r['rms'])<.5 else 'black')
                im=ax.imshow(values,norm=norm,cmap='viridis',aspect='auto')
                ax.set(title=f'{NAMES[m]}: {LABELS[j]}',xticks=range(len(deltas)),xticklabels=deltas,
                       yticks=range(len(foils)),yticklabels=foils,xlabel='Delta slice',ylabel='Foil (cm)')
            fig.colorbar(im,ax=grid[:,j].tolist(),label='RMS residual',shrink=.6)
        fig.suptitle(f'{p}: absolute RMS on identical events; color scale shared vertically')
        fig.savefig(out/f'plots/{p}_foil_delta.png',dpi=140);plt.close(fig)
