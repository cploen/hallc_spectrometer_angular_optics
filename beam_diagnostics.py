"""Small beam-search portfolio: search path, physical residuals, and foil/delta coverage."""
import numpy as np
from core_sample import write_tsv,read_tsv
from elastic_diagnostics import residuals,groups,LABELS
from elastic_net import TARGETS,macro_mse

NAMES = dict(svd='Full basis',start='Seed',beam='Beam')


def stats(v):
    return dict(n=len(v),bias=float(v.mean()),spread=float(v.std()),rms=float(np.sqrt(np.mean(v*v))),
                median=float(np.median(np.abs(v))),p90=float(np.percentile(np.abs(v),90)))


def report(out,a,truth,predictions,offset,folds,results,cells,qa_min):
    rr = {m:residuals(a,truth,p,offset) for m,p in predictions.items()}
    targets = (*TARGETS,'ztar')
    details,coverage,summary,byfold = [],[],[],[]
    for m,values in rr.items():
        for j,t in enumerate(targets):
            summary.append(dict(model=m,target=t,cell_mse=float(macro_mse(values[:,j],cells)),**stats(values[:,j])))
            for f in np.unique(folds):
                v = values[folds==f,j]
                byfold.append(dict(model=m,target=t,fold=int(f),cell_mse=float(macro_mse(v,cells[folds==f])),**stats(v)))
    for g,mask in groups(a,np.full(len(folds),'development')):
        n = int(sum(mask));coverage.append(dict(g,n=n,qa_low=n<qa_min))
        for m,values in rr.items():
            for j,t in enumerate(targets):
                details.append(dict(g,model=m,target=t,qa_low=n<qa_min,**stats(values[mask,j])))
    write_tsv(out/'tsv/summary.tsv',summary)
    write_tsv(out/'tsv/fold_scores.tsv',byfold)
    write_tsv(out/'tsv/residuals.tsv',details)
    write_tsv(out/'tsv/coverage.tsv',coverage)
    np.savez_compressed(out/'residuals.npz',rungroup=a['rungroup'],entry=a['entry'],fold=folds,
                        targets=np.array(targets),**rr)
    plots(out,rr,summary,details,results)
    text = '# Beam-search development results\n\n'
    text += '| Target | Seed terms | Beam terms | Seed MSE | Beam MSE | Full-basis MSE | Stop |\n|---|---:|---:|---:|---:|---:|---|\n'
    for j,t in enumerate(TARGETS):
        s = {r['model']:r['cell_mse'] for r in summary if r['target']==t}
        r = results[j]
        text += f"| {t} | {r['start']['terms']} | {r['chosen']['terms']} | {s['start']:.6g} | {s['beam']:.6g} | {s['svd']:.6g} | {r['stop']} |\n"
    text += ('\nMSE is averaged equally over populated physical foil/delta cells. Its units are mrad² '
             'for angles and cm² for ytar. Median and P90 describe event-pooled absolute residuals. '
             'A hole with six events does not receive the same weight as an entire well-populated cell, '
             'and sparse-hole flags never veto a candidate. No events are omitted from scoring.\n\n'
             'plots/search.png shows best score reached, terms, and the worst fold condition number; '
             'stars mark the final choice, which can differ from the lowest-MSE candidate because of '
             'the configured sparsity allowance. The full-basis reference always uses these same events '
             'and fold scaling, not the historical ROOT normal-equation solve on a different sample. '
             'The constant is counted as a term but excluded from the centered slope-block condition number.\n\n'
             'plots/residuals.png shows signed residual distributions. plots/foil_delta.png shows absolute '
             'RMS in mrad or cm, with a shared color scale for the three models within each target. '
             'RMS includes bias; spread in tsv/residuals.tsv is standard deviation after subtracting the '
             'mean, not a fitted Gaussian detector resolution. ztar uses the saved-xtar reconstruction '
             'diagnostic and is not a new replay.\n\n'
             'tsv/residuals.tsv gives bias, spread, RMS, median, P90 and N down to individual setting/hole. '
             'qa_low marks fewer than ten events and is informational. tsv/summary.tsv and fold_scores.tsv '
             'give overall and individual-fold scores; conditioning.tsv gives ranks and conditioning by '
             'fold, final_fit.tsv the final full-training solves. residuals.npz preserves exact event IDs '
             'and physical-unit residuals for all three models, ordered xptar, ytar, yptar, ztar.\n\n'
             '**These are adaptive development scores.** Candidate term choices use the same three folds '
             'that score the search, and the original elastic-net seeds also depended on these training '
             'events. Each prediction uses coefficients fitted without that event, but the basis selection '
             'is not independent of the scored events. There are no confidence intervals or claims of '
             'independent generalization. Protected core and noncore pools remain closed.\n\n'
             'Each target starts from all distinct saved bases for its chosen elastic-net setting. '
             'The Seed reference is the best of those fixed bases under this shared-fold scoring, '
             'not the earlier pooled result that used a different seed in each fold. No averaging or '
             'intersection of the seed masks is used.\n\n'
             'matrices/beam.dat is the chosen basis refitted on all balanced training events. '
             'matrices/start.dat and svd.dat provide corresponding seed and full-basis refits. '
             'seeds/fold*.dat preserve the actual input seed matrices; seed.dat preserves the fixed '
             'transport terms. Delta coefficients and xtar-dependent terms are unchanged. '
             'No matrix has been installed in replay or evaluated on protected events.\n')
    (out/'RESULTS.md').write_text(text)


def plots(out,rr,summary,details,results):
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    conditioning = read_tsv(out/'tsv/conditioning.tsv')
    fig,axes = plt.subplots(3,3,figsize=(13,10),constrained_layout=True)
    for j,t in enumerate(TARGETS):
        result = results[j];history = result['history'];chosen = result['chosen']
        full = next(r for r in summary if r['target']==t and r['model']=='svd')['cell_mse']
        full_cond = max(float(r['condition']) for r in conditioning if r['target']==t and r['model']=='svd')
        for i,(key,label) in enumerate((('score','Mean foil/delta MSE'),('terms','Terms including constant'),('condition','Worst fold condition'))):
            ax = axes[i,j];ax.plot([r['step'] for r in history],[r[key] for r in history],'.-',label='Best MSE reached')
            ax.scatter(history[-1]['step'],chosen[key],marker='*',s=140,color='C3',label='Chosen basis')
            if i == 0: ax.axhline(full,ls='--',color='black',label='Full basis')
            if i == 2: ax.set_yscale('log');ax.axhline(full_cond,ls='--',color='black')
            ax.set(xlabel='Search round',ylabel=label,title=LABELS[j]);ax.legend(fontsize=8)
    fig.suptitle('Beam-search development scores; full basis is a reference, not a pass/fail threshold')
    fig.savefig(out/'plots/search.png',dpi=140);plt.close(fig)
    fig,axes = plt.subplots(1,4,figsize=(17,4),constrained_layout=True)
    for j,ax in enumerate(axes):
        bins = np.histogram_bin_edges(np.concatenate([v[:,j] for v in rr.values()]),bins=150)
        for m,v in rr.items(): ax.hist(v[:,j],bins=bins,histtype='step',label=NAMES[m],log=True)
        ax.set(xlabel=LABELS[j],ylabel='Events',title='Signed reconstruction residual');ax.legend(fontsize=8)
    fig.savefig(out/'plots/residuals.png',dpi=140);plt.close(fig)
    cells = [r for r in details if r['level']=='foil_delta']
    foils = sorted({r['zfoil'] for r in cells});deltas = sorted({r['ndel'] for r in cells})
    lookup = {(r['model'],r['target'],r['zfoil'],r['ndel']):r for r in cells}
    fig,axes = plt.subplots(3,4,figsize=(17,10),constrained_layout=True)
    for j,t in enumerate((*TARGETS,'ztar')):
        vmax = max(r['rms'] for r in cells if r['target']==t)
        norm = Normalize(0,vmax or 1.)
        for i,m in enumerate(NAMES):
            ax = axes[i,j];values = np.full((len(foils),len(deltas)),np.nan)
            for zi,z in enumerate(foils):
                for di,d in enumerate(deltas):
                    r = lookup.get((m,t,z,d))
                    if r:
                        values[zi,di] = r['rms']
                        ax.text(di,zi,f"{r['rms']:.3f}\nN={r['n']}",ha='center',va='center',fontsize=6,
                                color='white' if r['rms']<.5*vmax else 'black')
            im = ax.imshow(values,cmap='viridis',norm=norm,aspect='auto')
            ax.set(title=f'{NAMES[m]}: {LABELS[j]}',xticks=range(len(deltas)),xticklabels=deltas,
                   yticks=range(len(foils)),yticklabels=foils,xlabel='Delta slice',ylabel='Foil (cm)')
        fig.colorbar(im,ax=axes[:,j].tolist(),label='RMS residual',shrink=.6)
    fig.suptitle('Absolute development RMS; color scales shared vertically. RMS includes systematic offsets.')
    fig.savefig(out/'plots/foil_delta.png',dpi=140);plt.close(fig)
