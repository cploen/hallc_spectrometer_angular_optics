"""Frozen-matrix attribution tables, central residual zooms, and separate tail views."""
import json
import numpy as np
from core_sample import write_tsv
from elastic_diagnostics import groups,LABELS
from elastic_net import macro_mse

MODELS=('old','gmm','core','enet','beam')
NAMES=dict(old='Old replay + offsets',gmm='GMM SVD',core='Full core SVD',
           enet='EN-selected SVD',beam='Beam')
TARGETS=('xptar','ytar','yptar','ztar')
POOLS=('protected_core','protected_noncore','surplus_core')
PAIRS=(('gmm','old'),('core','gmm'),('enet','core'),('beam','enet'),('beam','core'),('beam','gmm'))
PLOT_GUIDE='''# Reading the comparison plots

- `*_center.png`: signed residuals near zero. Each bin shows a percentage of
  **all events in that pool**, with common bins for every matrix. The window is
  ±beam P90; the legend in each panel reports the percentage visible. No curve
  is recentered. A narrower, taller peak is not by itself evidence of less bias.
- `*_tails.png`: percentage of events with an absolute residual at least as
  large as the horizontal threshold. For example, 1% at 3 mrad means 1% of the
  pool has |residual| ≥3 mrad. Lower is better. The horizontal axis is linear,
  from zero to the largest P99.9 among the five matrices; the vertical axis is
  logarithmic, labelled in percent with 10%, 1%, and 0.1% guides.
- `*_extremes.png`: a separate full-range tail view, with logarithmic axes.
  This retains the rare extremes beyond the main tail view. Both tail figures
  use the entire pool as denominator; zooming never renormalizes the data.
- `*_foil_delta.png`: absolute RMS, all five matrices, all populated physical
  foils and delta slices. Each target uses a common color scale across matrices.
- `*_change.png`: RMS change, 100 × (new/reference − 1), for each foil/delta
  cell. Blue/negative means smaller RMS; red/positive means larger RMS. Colors
  saturate at ±20%; printed values remain exact to the displayed precision.
  Zero reference RMS or absent cells have no defined percentage and show “—”.

Angular residuals are at the target in **mrad**; ytar and derived ztar are in
**cm**. N is the number of evaluation events, identical for both matrices in
each comparison. An asterisk marks N<10 for caution; it never excludes events.
RMS includes bias and is not a fitted Gaussian resolution. Consult bias and
mean-subtracted spread in the unchanged TSVs as well. The historical old-matrix
offset convention remains provisional; changing the plots does not resolve it.

Protected cores and protected noncore events are separate. Surplus cores are
development coverage. These figures change presentation only, not matrices,
offsets, event membership, or residuals.
'''


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
    plots(out,pool,rr,details)
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
           'Read plots/README.md for the central residual, tail, extreme-tail and foil/delta figures. '
           'Central bins and tail probabilities are percentages of the entire pool. Percent-change '
           'heat maps complement absolute RMS maps. Full bias, spread, median and P90 down to '
           'setting/hole are in tsv/residuals.tsv; qa_low flags N<10 and never changes the sample.\n\n'
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



def tail_percent(sorted_absolute, thresholds):
    """Percentage at or beyond each threshold, including ties and all pool events."""
    return 100*(len(sorted_absolute)-np.searchsorted(sorted_absolute,thresholds,side='left'))/len(sorted_absolute)


def rms_change(value, reference):
    return 100*(value/reference-1) if reference>0 else np.nan


def heatmaps(out,p,details):
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    local=[r for r in details if r['pool']==p and r['level']=='foil_delta']
    if not local:return
    foils=sorted({r['zfoil'] for r in local});deltas=sorted({r['ndel'] for r in local})
    lookup={(r['model'],r['target'],r['zfoil'],r['ndel']):r for r in local}
    for change in (False,True):
        rows=PAIRS if change else [(m,None) for m in MODELS]
        fig,grid=plt.subplots(len(rows),4,figsize=(17,3*len(rows)),constrained_layout=True)
        for j,t in enumerate(TARGETS):
            norm=Normalize(-20,20) if change else Normalize(0,max(r['rms'] for r in local if r['target']==t) or 1)
            cmap=plt.get_cmap('RdBu_r' if change else 'viridis').with_extremes(bad='0.9')
            for i,(m,base) in enumerate(rows):
                ax=grid[i,j];values=np.full((len(foils),len(deltas)),np.nan)
                for zi,z in enumerate(foils):
                    for di,d in enumerate(deltas):
                        r=lookup.get((m,t,z,d));b=lookup.get((base,t,z,d))
                        v=rms_change(r['rms'],b['rms']) if change and r and b else (r['rms'] if r and not change else np.nan)
                        values[zi,di]=v
                        label=f'{v:+.1f}%' if change else f'{v:.3f}'
                        if not np.isfinite(v):label='—'
                        if r:label+=f"\nN={r['n']}{'*' if r['n']<10 else ''}"
                        rgba=cmap(norm(v)) if np.isfinite(v) else (.9,.9,.9,1)
                        luminance=np.dot(rgba[:3],[.2126,.7152,.0722])
                        ax.text(di,zi,label,ha='center',va='center',fontsize=6,color='black' if luminance>.5 else 'white')
                im=ax.imshow(values,norm=norm,cmap=cmap,aspect='auto')
                title=f'{NAMES[m]} vs {NAMES[base]}' if change else NAMES[m]
                ax.set_title(f'{title}\n{LABELS[j]}',fontsize=10)
                ax.set(xticks=range(len(deltas)),xticklabels=deltas,
                       yticks=range(len(foils)),yticklabels=foils,xlabel='Delta slice',ylabel='Foil (cm)')
            if not change:fig.colorbar(im,ax=grid[:,j].tolist(),label='RMS residual',shrink=.6)
        if change:fig.colorbar(im,ax=grid.ravel().tolist(),label='RMS change (%)',extend='both',shrink=.6)
        title=('RMS change: blue = smaller; red = larger; colors saturate at ±20%' if change else
               'absolute RMS on identical events; color scale shared vertically')
        fig.suptitle(f'{p}: {title}\n* N<10; no events excluded')
        suffix='change' if change else 'foil_delta'
        fig.savefig(out/f'plots/{p}_{suffix}.png',dpi=140);plt.close(fig)


def plots(out,pool,rr,details):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, FixedLocator, NullFormatter
    colors=('0.35','C1','C0','C4','C2');styles=(':','--','-','-.','-')
    (out/'plots').mkdir(exist_ok=True)
    (out/'plots/README.md').write_text(PLOT_GUIDE)
    for p in POOLS:
        sel=pool==p
        if not sel.any():continue
        n=int(np.count_nonzero(sel));values={m:rr[m][sel] for m in MODELS}
        print(f'Plot {p}: {n:,} events',flush=True)
        center,axes=plt.subplots(1,4,figsize=(18,5),constrained_layout=True)
        tails,tax=plt.subplots(1,4,figsize=(18,5),constrained_layout=True)
        extremes,eax=plt.subplots(1,4,figsize=(18,5),constrained_layout=True)
        for j,(ax,tx,ex) in enumerate(zip(axes,tax,eax)):
            limit=max(float(np.percentile(np.abs(values['beam'][:,j]),90)),1e-9)
            tail_limit=max(max(float(np.percentile(np.abs(v[:,j]),99.9)) for v in values.values()),1e-9)
            bins=np.linspace(-limit,limit,81)
            thresholds=np.linspace(0,tail_limit,701)
            fractions=[]
            for m,color,style in zip(MODELS,colors,styles):
                v=values[m][:,j];count,_=np.histogram(v,bins)
                ax.stairs(100*count/len(v),bins,color=color,ls=style,label=NAMES[m])
                fractions.append(f'{NAMES[m]}: {100*np.mean(np.abs(v)<=limit):.1f}% visible')
                sorted_v=np.sort(np.abs(v))
                tx.plot(thresholds,tail_percent(sorted_v,thresholds),color=color,ls=style,label=NAMES[m])
                # More samples in the rare tail; probabilities remain exact at each threshold.
                take=np.unique(np.rint(len(v)-np.geomspace(len(v),1,min(700,len(v)))).astype(int))
                absolute=np.unique(sorted_v[take]);absolute=absolute[absolute>0]
                if len(absolute):
                    ex.step(absolute,tail_percent(sorted_v,absolute),where='pre',color=color,ls=style,label=NAMES[m])
            ax.set(xlabel='Residual '+LABELS[j],ylabel='Events per bin (%)',title='Central signed residual',xlim=(-limit,limit))
            ax.text(.02,.98,'\n'.join(fractions),va='top',transform=ax.transAxes,fontsize=6)
            for view in (tx,ex):
                view.set(xlabel='Absolute residual '+LABELS[j],ylabel='Events at or beyond threshold (%)',yscale='log')
                for y in (10,1,.1):view.axhline(y,color='0.7',ls=':',lw=.8,zorder=0)
                ticks=[v for v in (100,10,1,.1,.01,.001,.0001) if v>=min(.1,50/n)]
                view.yaxis.set_major_locator(FixedLocator(ticks))
                view.yaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{v:g}%'))
                view.yaxis.set_minor_formatter(NullFormatter())
            tx.set(xlim=(0,tail_limit),ylim=(.1,105),title='Main tails: linear residual axis')
            ex.set(xscale='log',ylim=(min(.05,50/n),105),title='Rare extremes: full residual range')
        for fig,axs,suffix,title in (
            (center,axes,'center','central window = ±beam P90; no recentering'),
            (tails,tax,'tails','linear range to largest P99.9; rarer tails in the extremes figure'),
            (extremes,eax,'extremes','full-range tails; both axes logarithmic')):
            fig.legend(*axs[0].get_legend_handles_labels(),loc='outside lower center',ncol=5,fontsize=8)
            fig.suptitle(f'{p}: N={n:,}; {title}\nPercentages use all events in this pool')
            fig.savefig(out/f'plots/{p}_{suffix}.png',dpi=140);plt.close(fig)
        heatmaps(out,p,details)
