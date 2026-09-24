"""Reserved-sample plots for the frozen Huber expansion experiment."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from elastic_net import TARGETS

STAGES=('balanced','half_surplus','all_core','half_shoulders','all_supported')
LABELS=('Balanced','+½ surplus','All cores','+½ shoulders','All supported')
UNITS=('mrad','cm','mrad')


def lookup(stats):
    return {(r['model'],r['pool'],r['group_kind'],r['group'],r['target']):r for r in stats}


def plots(out,stats,counts):
    table=lookup(stats)
    fig,axes=plt.subplots(2,3,figsize=(12.5,7),layout='constrained')
    for row,pool in enumerate(('core','shoulder')):
        for j,target in enumerate(TARGETS):
            ax=axes[row,j]
            baseline=table['balanced_squared_full',pool,'overall','all',target]['rms']
            for loss,color in [('squared','#477aa3'),('huber','#b34b36')]:
                values=[100*(1-table[f'{s}_{loss}_full',pool,'overall','all',target]['rms']/baseline) for s in STAGES]
                ax.plot(range(5),values,'o-',color=color,label='Squared error' if loss=='squared' else 'Smooth Huber')
            ax.axhline(0,color='.5',lw=.8)
            ax.set_xticks(range(5),LABELS,rotation=25,ha='right',fontsize=9)
            ax.set_title(f'{target} · reserved {pool}s')
            ax.set_ylabel('RMS reduction (%)')
            ax.grid(axis='y',alpha=.2)
            if row==0 and j==0: ax.legend(fontsize=9)
    fig.suptitle('More training data: same full polynomial, same reserved events\nPositive = lower RMS than balanced squared-error fit',fontsize=13)
    fig.savefig(out/'plots/expansion.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for i,pool in enumerate(('core','shoulder')):
        for j,target in enumerate(TARGETS):
            ax=axes[i,j]
            cells=sorted({r['group'] for r in stats if r['group_kind']=='foil_delta' and r['pool']==pool},
                         key=lambda s:tuple(map(float,s.split(':'))))
            foils=sorted({float(c.split(':')[0]) for c in cells})
            slices=sorted({int(c.split(':')[1]) for c in cells})
            values=np.full((len(foils),len(slices)),np.nan)
            for c in cells:
                z,d=map(float,c.split(':'))
                b=table['balanced_squared_full',pool,'foil_delta',c,target]['rms']
                v=table['all_supported_huber_full',pool,'foil_delta',c,target]['rms']
                values[foils.index(z),slices.index(int(d))]=100*(1-v/b)
            lim=max(1.,float(np.nanmax(np.abs(values))))
            im=ax.imshow(values,cmap='RdBu',vmin=-lim,vmax=lim,aspect='auto')
            for (r,c),v in np.ndenumerate(values):
                if np.isfinite(v): ax.text(c,r,f'{v:+.1f}',ha='center',va='center',fontsize=9,
                                          color='white' if abs(v)>.6*lim else 'black')
            ax.set_xticks(range(len(slices)),slices);ax.set_yticks(range(len(foils)),[f'{z:g}' for z in foils])
            ax.set_xlabel('Frozen delta slice');ax.set_ylabel('Foil z (cm)')
            ax.set_title(f'{target} · {pool}s');fig.colorbar(im,ax=ax,label='RMS reduction (%)',shrink=.8)
    fig.suptitle('All supported events + smooth Huber versus balanced squared error\nSame full basis; each foil/delta cell evaluated separately',fontsize=13)
    fig.savefig(out/'plots/foil_delta.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(11.5,6.5),layout='constrained')
    for i,pool in enumerate(('core','shoulder')):
        for j,target in enumerate(TARGETS):
            ax=axes[i,j]
            for kind,label,color in [('full','Full basis','#477aa3'),('enet_basis','EN-selected basis','#7c6a9c'),('beam_basis','Beam basis','#418c6b')]:
                base=table['balanced_squared_full',pool,'overall','all',target]['rms']
                models=[f'balanced_squared_{kind}',f'all_core_huber_{kind}',f'all_supported_huber_{kind}']
                vals=[100*(1-table[m,pool,'overall','all',target]['rms']/base) for m in models]
                ax.plot(range(3),vals,'o-',label=label,color=color)
            ax.axhline(0,color='.5',lw=.8);ax.grid(axis='y',alpha=.2)
            ax.set_xticks(range(3),['Balanced\nsquared','All cores\nHuber','All supported\nHuber'])
            ax.set_title(f'{target} · {pool}s');ax.set_ylabel('RMS reduction (%)')
            if i==0 and j==0: ax.legend(fontsize=8)
    fig.suptitle('Reuse of frozen elastic-net-selected and beam-search bases\nNo new term selection on the reserved sample',fontsize=13)
    fig.savefig(out/'plots/bases.png',dpi=180);plt.close(fig)


def report(stats,intervals,counts,manifest):
    table=lookup(stats)
    text='# Smooth Huber and training-sample expansion\n\n'
    text+='Real HMS 6.667 GeV frozen event data. Candidate matrices are exported but not installed in replay.\n\n'
    text+='## Experimental definition\n\n'
    text+='The regression loss is the smooth (pseudo-)Huber loss, '
    text+='`rho(u) = c² (sqrt(1+(u/c)²)-1)`, with `u = residual / sigma` and `c = 1.5`. '
    text+='`sigma` is 1.4826 times the residual MAD from the balanced-training squared-error fit, '
    text+='frozen for every later fit. Constants are fitted; xtar-dependent and delta coefficients '
    text+='are inherited from the archived seed. No foil, slice or hole population weights are applied. '
    text+='The regression loss follows [SciPy pseudo_huber](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pseudo_huber.html); '
    text+='the classification loss called modified_huber is not used.\n\n'
    text+='The main experiment retains all 210 terms (including the constant). '
    text+='Half increments use a reproducible event-ID hash and retain the entire preceding sample. '
    text+='Additional cores mostly increase statistics within existing accepted regions; shoulders also broaden the occupied sample. '
    text+='The expansion stays inside the original accepted labels and delta range; it cannot recover events absent from the frozen trees.\n\n'
    text+='| Training sample | Events |\n|---|---:|\n'
    for k in STAGES: text+=f'| {k} | {counts[k]:,} |\n'
    text+='\nResidual scales (xptar mrad, ytar cm, yptar mrad): '+', '.join(f'{v:.6g}' for v in manifest['residual_scales_physical'])+'.\n\n'
    text+='## Common reserved-sample RMS\n\n'
    for pool in ('core','shoulder'):
        text+=f'### Reserved {pool}s\n\n| Training | Loss | xptar (mrad) | ytar (cm) | yptar (mrad) |\n|---|---|---:|---:|---:|\n'
        for k in STAGES:
            for loss in ('squared','huber'):
                vals=[table[f'{k}_{loss}_full',pool,'overall','all',t]['rms'] for t in TARGETS]
                text+=f'| {k} | {loss} | '+ ' | '.join(f'{v:.6f}' for v in vals)+' |\n'
        text+='\n'
    text+='## Expanded Huber versus balanced squared error\n\n'
    text+='Positive gain means lower RMS. The intervals below are paired 500-replicate bootstrap intervals '
    text+='over run/foil/delta/hole clusters, stratified by physical foil/delta. They describe evaluation-sample '
    text+='variation conditional on these fitted models, not training uncertainty or systematic calibration error.\n\n'
    text+='| Pool | Target | RMS reduction (%) | Conditional 95% interval | Cell-macro RMS reduction (%) | Improved cells | Worst cell change (%) |\n|---|---|---:|---|---:|---:|---:|\n'
    for pool in ('core','shoulder'):
        for target in TARGETS:
            r=next(r for r in intervals if r['model']=='all_supported_huber_full' and r['pool']==pool and r['target']==target)
            cells=[v for v in stats if v['model']=='balanced_squared_full' and v['pool']==pool and v['target']==target and v['group_kind']=='foil_delta']
            old=np.array([v['rms'] for v in cells]);new=np.array([table['all_supported_huber_full',pool,'foil_delta',v['group'],target]['rms'] for v in cells])
            macro=100*(1-np.sqrt(np.mean(new**2)/np.mean(old**2)))
            gain=100*(1-new/old)
            text+=f'| {pool} | {target} | {r["rms_gain_percent"]:.3f} | [{r["low95"]:.3f}, {r["high95"]:.3f}] | {macro:.3f} | {sum(gain>0)}/{len(gain)} | {min(gain):.3f} |\n'
    text+='\n## Elastic net and beam-search controls\n\n'
    text+='The existing elastic-net-selected and beam10 term supports are reused and refitted on all cores '
    text+='and on cores plus shoulders with both losses. The elastic-net-selected basis here has unpenalized '
    text+='coefficients; it is not a new penalized elastic-net fit. No new beam search or hyperparameter '
    text+='optimization is claimed. See `plots/bases.png` and the complete `tsv/metrics.tsv`.\n\n'
    text+='## Checks and limitations\n\n'
    text+='Frozen input hashes and balanced event IDs were verified against the archived beam run. '
    text+='The locally rebuilt balanced SVD predictions reproduce the saved SVD matrix; maximum absolute '
    text+='differences in the reported physical units are '+str(manifest['balanced_prediction_max_abs_difference'])+'. '
    text+='All Huber optimizations pass an explicit gradient check; candidate matrix files reproduce in-memory predictions.\n\n'
    text+='All variants are fixed before examining the common reserved residuals. The reserved sample was already '
    text+='examined during earlier work, so this is diagnostic evidence, not a new independent final validation. '
    text+='Density labels themselves were estimated from development events. Unsupported and blocked labels are '
    text+='excluded. Shoulder residuals remain conditional on their inherited labels. Original geometry and saved '
    text+='xtar are held fixed; full iterative replay and independent physical calibration are untested. '
    text+='RMS includes bias and is not a Gaussian resolution. Bias, standard deviation, central 68% half-width, '
    text+='95th absolute residual, and sparse-hole statistics are retained in the TSV files.\n\n'
    text+='Files: `plots/expansion.png`, `plots/foil_delta.png`, `plots/bases.png`, `tsv/coverage.tsv`, '
    text+='`tsv/metrics.tsv`, `tsv/holes.tsv`, `tsv/paired_intervals.tsv`, `tsv/convergence.tsv`, '
    text+='`tsv/membership.tsv`, candidate `matrices/*.dat`, and the input/output manifest.\n'
    return text
