"""Common-event diagnostics for equal-foil Huber and squared-error fits."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core_sample import read_tsv,write_tsv
from elastic_net import TARGETS

NEW=tuple(f'{p}_{s}_equalfoil' for p in ('all_core','all_supported') for s in ('squared','huber'))
NAMES=dict(balanced_squared_full='Original balanced SVD',
           all_core_squared_full='Unweighted cores, squared',all_core_huber_full='Unweighted cores, Huber',
           all_supported_squared_full='Unweighted cores + shoulders, squared',
           all_supported_huber_full='Unweighted cores + shoulders, Huber',
           all_core_squared_equalfoil='Equal-foil cores, squared',all_core_huber_equalfoil='Equal-foil cores, Huber',
           all_supported_squared_equalfoil='Equal-foil cores + shoulders, squared',
           all_supported_huber_equalfoil='Equal-foil cores + shoulders, Huber')
BASE='balanced_squared_full'


def products(out,stats,intervals,weights,manifest):
    table={(r['model'],r['pool'],r['group_kind'],r['group'],r['target']):r for r in stats}
    holes=read_tsv(out/'tsv/holes.tsv')
    grouped={}
    for r in holes:
        if int(r['n'])>=10:
            grouped.setdefault((r['model'],r['pool'],r['target']),{})[r['hole']]=r
    models=list(dict.fromkeys(r['model'] for r in stats));summary=[]
    for pool in ('core','shoulder','supported'):
        for target in TARGETS:
            cell=[r for r in stats if r['model']==BASE and r['pool']==pool and r['target']==target and r['group_kind']=='foil_delta']
            reference=np.array([r['rms'] for r in cell])
            base_h=grouped[BASE,pool,target]
            for model in models:
                overall=table[model,pool,'overall','all',target]
                current=np.array([table[model,pool,'foil_delta',c['group'],target]['rms'] for c in cell])
                h=grouped[model,pool,target]
                if h.keys()!=base_h.keys(): raise ValueError('Within-hole comparison membership changed')
                n=np.array([int(v['n']) for v in h.values()]);s=np.array([float(v['spread']) for v in h.values()])
                hole_bias=np.array([float(v['bias']) for v in h.values()])
                hole_rms=np.array([float(v['rms']) for v in h.values()])
                base_s=np.array([float(base_h[k]['spread']) for k in h])
                gains=100*(1-current/reference)
                baseline=table[BASE,pool,'overall','all',target]
                width=float(np.sqrt(np.average(s*s,weights=n)))
                width_ref=float(np.sqrt(np.average(base_s*base_s,weights=n)))
                summary.append(dict(model=model,pool=pool,target=target,n=overall['n'],
                    rms=overall['rms'],bias=overall['bias'],spread=overall['spread'],width68=overall['width68'],p95=overall['p95'],
                    pooled_rms_gain=100*(1-overall['rms']/baseline['rms']),
                    macro_rms=float(np.sqrt(np.mean(current**2))),macro_gain=100*(1-np.sqrt(np.mean(current**2)/np.mean(reference**2))),
                    improved_cells=int(sum(gains>0)),cells=len(cell),worst_cell_gain=float(min(gains)),
                    within_hole_spread=width,within_hole_gain=100*(1-width/width_ref),
                    hole_centroid_rms=float(np.sqrt(np.average(hole_bias**2,weights=n))),
                    same_holes_rms=float(np.sqrt(np.average(hole_rms**2,weights=n))),
                    equal_hole_spread=float(np.sqrt(np.mean(s*s))),
                    improved_holes=int(sum(s<base_s)),holes=len(h),within_hole_events=int(sum(n))))
    write_tsv(out/'tsv/summary.tsv',summary)
    index={(r['model'],r['pool'],r['target']):r for r in summary}
    # Four prespecified comparisons against the same original balanced baseline.
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    labels=['Cores\nsquared','Cores\nHuber','Cores +\nshoulders\nsquared','Cores +\nshoulders\nHuber']
    for i,pool in enumerate(('core','shoulder')):
        for j,(metric,title) in enumerate((('pooled_rms_gain','Pooled RMS'),('macro_gain','Equal foil/delta-cell RMS'),('within_hole_gain','Within-hole standard deviation'))):
            ax=axes[i,j]
            for target,color in zip(TARGETS,('#2874a6','#aa4838','#4b8b53')):
                ax.plot(range(4),[index[m,pool,target][metric] for m in NEW],'o-',label=target,color=color)
            ax.axhline(0,color='.5',lw=.8);ax.grid(axis='y',alpha=.2)
            ax.set_xticks(range(4),labels,fontsize=8);ax.set_ylabel('Reduction from original balanced SVD (%)')
            ax.set_title(f'{title}\nReserved {pool}s',fontsize=11)
            if i==0 and j==0: ax.legend(fontsize=9)
    fig.suptitle('Equal total base weight per foil, retaining every training event\nPositive = smaller errors; full 210-term basis and fixed Huber threshold',fontsize=13)
    fig.savefig(out/'plots/comparison.png',dpi=170);plt.close(fig)
    heatmap(out,table,stats,'all_supported_huber_equalfoil',BASE,'foil_delta.png',
            'Equal-foil Huber with cores + shoulders versus original balanced SVD')
    heatmap(out,table,stats,'all_supported_huber_equalfoil','all_supported_huber_full','weighting_effect.png',
            'Effect of foil weighting: equal-foil versus unweighted Huber, identical training events')
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,pool in zip(axes,('all_core','all_supported')):
        for target,color in zip(TARGETS,('#2874a6','#aa4838','#4b8b53')):
            r=[r for r in weights if r['model']==pool+'_huber_equalfoil' and r['target']==target and r['group_kind']=='foil']
            ax.plot([float(v['group']) for v in r],[100*v['irls_share'] for v in r],'o-',label=target,color=color)
        ax.axhline(20,color='.4',ls='--',label='Base weight: 20% each')
        ax.set(xlabel='Physical foil z (cm)',ylabel='Share of total IRLS weight (%)',title=pool.replace('_',' '),xticks=[-8,-3,0,3,8])
        ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Foil contributions after Huber residual weighting\nBase weight × residual reduction factor; not a full leverage/influence measure',fontsize=12)
    fig.savefig(out/'plots/foil_weights.png',dpi=170);plt.close(fig)
    text='# Equal-foil weighting with smooth Huber\n\n'
    text+='Four prespecified fits were completed on the same development and reserved events as the unweighted trial. '
    text+='All cores: 310,513 training events. Cores plus supported shoulders: 396,320. '
    text+='Evaluation: 77,864 reserved cores and 22,532 reserved shoulders. No event was discarded to balance a foil.\n\n'
    text+='## Weighting definition\n\n'
    text+='Each event receives `w_i = N / (5 N_foil)`, so every physical foil has exactly 20% of the total base weight. '
    text+='Minimize `sum(w_i rho(r_i / sigma)) / sum(w_i)`; the weight is outside the loss, so it does not change the residual threshold. '
    text+='Huber uses `c = 1.5` and the previous balanced-training MAD scales: '+', '.join(f'{v:.7g}' for v in manifest['residual_scales_physical'])+' (mrad, cm, mrad). '
    text+='All 210 polynomial terms, fitted constants, inherited xtar/delta coefficients, and geometry remain fixed in definition. '
    text+='Weights are constructed from training counts only. Delta slices and holes are not separately balanced in this first follow-up. '
    text+='The original balanced sample also applied slice/hole caps; equal-foil weights alone do not recreate that allocation.\n\n'
    selected=(BASE,'all_supported_huber_full',*NEW)
    for pool in ('core','shoulder'):
        text+=f'## Reserved {pool}s: RMS\n\n| Fit | xptar (mrad) | ytar (cm) | yptar (mrad) |\n|---|---:|---:|---:|\n'
        for m in selected:
            text+='| '+NAMES[m]+' | '+' | '.join(f'{index[m,pool,t]["rms"]:.6f}' for t in TARGETS)+' |\n'
        text+='\n'
    text+='## Changes relative to the original balanced SVD\n\nPositive means smaller error. Equal-cell RMS is the square root of the mean of 25 foil/delta-cell MSEs. '
    text+='Within-hole width is the square root of the event-weighted mean within-hole variance, after subtracting each '
    text+='run/foil/delta/hole mean; only common groups with at least 10 reserved events are used. It is a descriptive width, not fitted Gaussian resolution.\n\n'
    for pool in ('core','shoulder'):
        text+=f'### {pool}\n\n| Fit | Target | Pooled RMS gain (%) | Equal-cell RMS gain (%) | Improved cells | Within-hole width gain (%) |\n|---|---|---:|---:|---:|---:|\n'
        for m in NEW:
            for t in TARGETS:
                r=index[m,pool,t]
                text+=f'| {NAMES[m]} | {t} | {r["pooled_rms_gain"]:.3f} | {r["macro_gain"]:.3f} | {r["improved_cells"]}/25 | {r["within_hole_gain"]:.3f} |\n'
        text+='\n'
    text+='## Residual weights after foil balancing\n\n'
    text+='For smooth Huber the IRLS factor is `h_i = 1/sqrt(1+(r_i/(c sigma))²)`. '
    text+='The following shares sum `w_i h_i`; they describe weights in the estimating equation, not full coefficient influence (which also depends on the polynomial design). '
    text+='No second foil normalization is applied after Huber weighting.\n\n'
    for pool in ('all_core','all_supported'):
        text+=f'### {pool}\n\n| Foil (cm) | Base share | xptar IRLS share | ytar IRLS share | yptar IRLS share |\n|---|---:|---:|---:|---:|\n'
        for z in (-8,-3,0,3,8):
            values=[next(r for r in weights if r['model']==pool+'_huber_equalfoil' and r['target']==t and r['group_kind']=='foil' and float(r['group'])==z) for t in TARGETS]
            text+=f'| {z:+g} | 20.00% | '+' | '.join(f'{100*r["irls_share"]:.2f}%' for r in values)+' |\n'
        text+='\n'
    text+='### Shoulder contributions in the full supported sample\n\n| Target | Raw event share | Base-weight share | IRLS-weight share | Mean residual factor |\n|---|---:|---:|---:|---:|\n'
    for t in TARGETS:
        r=next(r for r in weights if r['model']=='all_supported_huber_equalfoil' and r['target']==t and r['group_kind']=='quality' and r['group']=='shoulder')
        text+=f'| {t} | {100*r["raw_share"]:.2f}% | {100*r["base_share"]:.2f}% | {100*r["irls_share"]:.2f}% | {r["mean_robust_factor"]:.3f} |\n'
    text+='\n## Files and limits\n\n'
    text+='[Comparison](plots/comparison.png) · [Foil/delta changes](plots/foil_delta.png) · [Effect of adding weights](plots/weighting_effect.png) · [Foil weights](plots/foil_weights.png).\n\n'
    text+='`tsv/summary.tsv` contains pooled, equal-cell and within-hole statistics for all nine reference/candidate fits. '
    text+='`tsv/metrics.tsv` includes biases, spreads and tails; `tsv/holes.tsv` preserves sparse groups; '
    text+='`tsv/weight_contributions.tsv` contains totals, foil shares, core/shoulder breakdowns and Kish effective counts. '
    text+='`training_weights.npz` stores exact IDs and base weights. `tsv/paired_intervals.tsv` provides conditional paired cluster-bootstrap '
    text+='intervals for pooled RMS changes against the original balanced fit; it does not include training or systematic uncertainty.\n\n'
    text+='All four candidates were specified before viewing their reserved scores. Source inputs and residual identities were verified, '
    text+='all optimizations converged, and exported matrices reproduce predictions. This follows earlier examination of the same reserved sample '
    text+='and is not fresh independent validation. Shoulder labels, target geometry and saved xtar remain inherited. '
    text+='No replay matrix was installed; these are fixed-input reconstruction comparisons.\n'
    (out/'RESULTS.md').write_text(text)


def heatmap(out,table,stats,model,base,filename,title):
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    foils=(-8,-3,0,3,8);ds=range(5)
    for i,pool in enumerate(('core','shoulder')):
        for j,t in enumerate(TARGETS):
            values=np.array([[100*(1-table[model,pool,'foil_delta',f'{z}:{d}',t]['rms']/table[base,pool,'foil_delta',f'{z}:{d}',t]['rms']) for d in ds] for z in foils])
            lim=max(1.,float(np.max(np.abs(values))))
            ax=axes[i,j];im=ax.imshow(values,cmap='RdBu',vmin=-lim,vmax=lim,aspect='auto')
            for (r,c),v in np.ndenumerate(values):
                ax.text(c,r,f'{v:+.1f}',ha='center',va='center',fontsize=9,color='white' if abs(v)>.6*lim else 'black')
            ax.set(xticks=list(ds),yticks=range(5),yticklabels=foils,xlabel='Frozen delta slice',ylabel='Foil z (cm)',title=f'{t} · {pool}s')
            fig.colorbar(im,ax=ax,label='RMS reduction (%)',shrink=.8)
    fig.suptitle(title+'\nPositive = lower RMS on identical reserved events',fontsize=12)
    fig.savefig(out/'plots'/filename,dpi=170);plt.close(fig)
