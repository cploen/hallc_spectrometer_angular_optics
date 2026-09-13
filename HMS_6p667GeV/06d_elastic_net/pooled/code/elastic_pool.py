"""Training-only, out-of-fold diagnostics for a small set of beam-search seeds."""
import numpy as np
from core_sample import write_tsv
from elastic_net import TARGETS, SCALES, macro_mse
from elastic_refit import subset_svd, native
from elastic_diagnostics import groups, LABELS


def scores(residual, cells):
    """Quantiles pool events; cell MSE preserves the existing foil/delta weighting."""
    return dict(n=len(residual), median=float(np.median(np.abs(residual))),
                p90=float(np.percentile(np.abs(residual), 90)),
                mse=float(np.mean(residual**2)), cell_mse=float(macro_mse(residual, cells)))


class Pool:
    def __init__(self, out, data, truth, folds, cells, terms, qa_min):
        self.out, self.data, self.truth = out, data, truth
        self.folds, self.cells, self.terms, self.qa_min = folds, cells, terms, qa_min
        self.fold = None
        self.rows, self.term_rows, self.values, self.coefficients = [], [], {}, {}

    def __call__(self, j, alpha, l1, enet, fit):
        selected = enet != 0
        beta, rank, singular = subset_svd(fit['x'], fit['y'][:,j], selected, fit['rcond'])
        mask, state = fit['mask'], fit['state']
        if not np.array_equal(mask, self.folds == self.fold):
            raise ValueError('Validation mask does not match saved fold')
        key = (j, alpha, l1)
        values = self.values.setdefault(key, np.full((len(self.folds),2), np.nan))
        if np.isfinite(values[mask]).any(): raise ValueError('Repeated validation events in pooled case')
        name = f'{TARGETS[j]}_a{alpha:g}_l{l1:g}_f{self.fold}'
        self.coefficients[name] = native(beta, state, j)
        for i,term in enumerate(self.terms):
            self.term_rows.append(dict(target=TARGETS[j], alpha=alpha, l1=l1, fold=self.fold,
                                       term=''.join(map(str,term)), selected=bool(i == 0 or selected[i-1])))
        for k,b in enumerate((fit['svd'][:,j], beta)):
            prediction = (fit['validation'] @ b)*state[3][j]+state[2][j]
            values[mask,k] = (prediction-self.truth[mask,j])*SCALES[j]
            s = fit['singular'] if k == 0 else singular
            condition = float(s[0]/s[-1]) if len(s) and s[-1] > 0 else np.inf
            self.rows.append(dict(target=TARGETS[j], alpha=alpha, l1=l1, fold=self.fold,
                                  model=('full_svd','refit_svd')[k],
                                  terms=len(self.terms) if k == 0 else int(sum(selected))+1,
                                  rank=int(np.sum(s > fit['rcond']*s[0])) if k == 0 else rank,
                                  condition=condition, **scores(values[mask,k],self.cells[mask])))
        write_tsv(self.out/'fold_scores.tsv',self.rows)
        write_tsv(self.out/'terms.tsv',self.term_rows)
        np.savez_compressed(self.out/'coefficients.npz',terms=self.terms,**self.coefficients)

    def finish(self, cfg):
        expected = [(j,a,l) for j in range(len(TARGETS)) for l in cfg['l1'] for a in cfg['alphas']]
        complete = [key for key in expected if key in self.values and np.isfinite(self.values[key]).all()]
        summary, details, archive = [], [], {}
        # Each residual is an untouched validation prediction from its assigned fold.
        for j,alpha,l1 in complete:
            values = self.values[j,alpha,l1]
            archive[f'{TARGETS[j]}_a{alpha:g}_l{l1:g}'] = values
            for k,model in enumerate(('full_svd','refit_svd')):
                rows = [r for r in self.rows if (r['target'],r['alpha'],r['l1'],r['model']) ==
                        (TARGETS[j],alpha,l1,model)]
                summary.append(dict(target=TARGETS[j],alpha=alpha,l1=l1,model=model,
                                    terms=float(np.mean([r['terms'] for r in rows])),
                                    terms_min=min(r['terms'] for r in rows), terms_max=max(r['terms'] for r in rows),
                                    condition=float(np.median([r['condition'] for r in rows])),
                                    condition_min=min(r['condition'] for r in rows),
                                    condition_max=max(r['condition'] for r in rows),
                                    **scores(values[:,k],self.cells)))
        coverage = []
        for g,mask in groups(self.data,np.full(len(self.folds),'validation')):
            coverage.append(dict(g,n=int(sum(mask)),qa_low=int(sum(mask)) < self.qa_min,
                                 **{f'fold{f}':int(sum(mask & (self.folds == f))) for f in np.unique(self.folds)}))
            if g['level'] != 'foil_delta': continue
            for j,alpha,l1 in complete:
                for k,model in enumerate(('full_svd','refit_svd')):
                    details.append(dict(g,target=TARGETS[j],alpha=alpha,l1=l1,model=model,
                                        **scores(self.values[j,alpha,l1][mask,k],self.cells[mask])))
        write_tsv(self.out/'coverage.tsv',coverage)
        if summary: write_tsv(self.out/'pooled.tsv',summary)
        if details: write_tsv(self.out/'foil_delta.tsv',details)
        np.savez_compressed(self.out/'residuals.npz',rungroup=self.data['rungroup'],entry=self.data['entry'],
                            fold=self.folds,models=np.array(['full_svd','refit_svd']),**archive)
        text = '# Reduced-basis seed study: pooled training validation\n\n'
        text += f'{len(complete)}/{len(expected)} target/penalty cases completed every fold. '
        text += 'Incomplete cases are excluded from pooled scores and plots.\n\n'
        text += ('Each balanced training event supplies exactly one out-of-fold residual per complete case. '
                 'Every fold independently selects terms, learns scaling, and refits scaled X directly by '
                 'unpenalized SVD on the other folds. No protected events are used.\n\n'
                 'Read plots/plateau.png first: pooled median absolute residual, mean foil/delta MSE, '
                 'P90 absolute residual, and the range of fitting-matrix condition numbers. '
                 'Faint points show separate folds; horizontal spans show their term-count range. '
                 'The condition-number point is the median of fold condition numbers, not a condition '
                 'number computed from pooled residuals. The constant is counted as a term but handled '
                 'separately from the centered/scaled nonconstant columns used for conditioning.\n\n'
                 'Median and P90 pool individual events. cell_mse averages MSE equally over populated '
                 'physical foil/delta cells, retaining the earlier scoring convention; mse in the tables '
                 'is the ordinary event-weighted value. Residual units are mrad for angles and cm for ytar. '
                 'Quantiles are computed from pooled residuals, never by averaging fold quantiles.\n\n'
                 'plots/foil_delta.png checks mean-squared residuals throughout acceptance. '
                 'foil_delta.tsv also supplies median and P90 by cell. coverage.tsv gives pooled and '
                 'individual-fold counts down to setting/hole, with qa_low as a reporting flag only. '
                 'Pooling does not create extra events in sparse holes.\n\n'
                 'Use modest, consistent changes across nearby term counts to nominate several seeds. '
                 'These four penalty settings can suggest a stable range, not establish a finely resolved '
                 'plateau. There is no automatic winning model or requirement to beat full-basis SVD. '
                 'Basis membership may differ between folds; terms.tsv and coefficients.npz preserve each '
                 'actual seed. Do not treat the average term count as a single fitted basis. Beam search '
                 'may reintroduce omitted terms and refit without penalties.\n\n'
                 'fold_scores.tsv preserves individual-fold metrics; pooled.tsv preserves exact pooled '
                 'metrics. residuals.npz stores complete-case residual pairs in [full_svd, refit_svd] order '
                 'with event IDs and folds. seed.dat, folds.tsv, code/, checkpoints.tsv and manifest.json '
                 'make the study reproducible. No final replay matrix is exported.\n')
        (self.out/'RESULTS.md').write_text(text)
        if summary: self.plots(summary,details)
        return len(complete)

    def plots(self, summary, details):
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm
        (self.out/'plots').mkdir()
        fig,axes = plt.subplots(4,3,figsize=(14,12),constrained_layout=True)
        keys = ('median','cell_mse','p90','condition')
        labels = ('Median |residual|','Mean foil/delta MSE','P90 |residual|','Condition of scaled X')
        for j,target in enumerate(TARGETS):
            rows = [r for r in summary if r['target'] == target and r['model'] == 'refit_svd']
            base = next((r for r in summary if r['target'] == target and r['model'] == 'full_svd'),None)
            if base is None: continue
            for c,r in enumerate(rows):
                fold = [v for v in self.rows if (v['target'],v['alpha'],v['l1'],v['model']) ==
                        (target,r['alpha'],r['l1'],'refit_svd')]
                color = f'C{c}'
                for i,key in enumerate(keys):
                    ax = axes[i,j]
                    ax.scatter([v['terms'] for v in fold],[v[key] for v in fold],color=color,alpha=.3,s=18)
                    ax.errorbar(r['terms'],r[key],xerr=[[r['terms']-r['terms_min']],[r['terms_max']-r['terms']]],
                                fmt='o',color=color,label=f"α={r['alpha']:g}, L1={r['l1']:g}")
                    if key == 'condition': ax.vlines(r['terms'],r['condition_min'],r['condition_max'],color=color)
            for i,key in enumerate(keys):
                ax = axes[i,j]
                ax.axhline(base[key],color='black',ls='--',lw=1,label='Full basis reference')
                ax.set(xlabel='Selected terms (constant included)',ylabel=labels[i],title=LABELS[j])
                if key == 'condition': ax.set_yscale('log')
            axes[0,j].legend(fontsize=7)
        fig.suptitle('Reduced-basis seed candidates: pooled validation; faint points show individual folds')
        fig.savefig(self.out/'plots/plateau.png',dpi=140);plt.close(fig)
        cases = sorted({(r['alpha'],r['l1']) for r in summary})
        foils,deltas = np.unique(self.data['ztarT']),np.unique(self.data['ndel'])
        lookup = {(r['target'],r['alpha'],r['l1'],r['model'],r['zfoil'],r['ndel']):r for r in details}
        fig,axes = plt.subplots(len(cases),3,figsize=(13,3*len(cases)),squeeze=False,constrained_layout=True)
        norm = TwoSlopeNorm(vmin=.5,vcenter=1,vmax=1.5)
        for i,(alpha,l1) in enumerate(cases):
            for j,target in enumerate(TARGETS):
                ax = axes[i,j]; values = np.full((len(foils),len(deltas)),np.nan)
                for zi,z in enumerate(foils):
                    for di,d in enumerate(deltas):
                        r = lookup.get((target,alpha,l1,'refit_svd',z,d))
                        b = lookup.get((target,alpha,l1,'full_svd',z,d))
                        if r and b and b['mse'] > 0:
                            values[zi,di] = r['mse']/b['mse']
                            ax.text(di,zi,f"{values[zi,di]:.2f}\nN={r['n']}",ha='center',va='center',fontsize=7)
                ax.imshow(values,norm=norm,cmap='coolwarm',aspect='auto')
                ax.set(title=f'{target}: α={alpha:g}, L1={l1:g}',xticks=range(len(deltas)),xticklabels=deltas,
                       yticks=range(len(foils)),yticklabels=foils,xlabel='Delta slice',ylabel='Foil (cm)')
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='coolwarm'),ax=axes.ravel().tolist(),
                     label='Pooled refit MSE / full-basis MSE (reference, not a pass/fail test)',extend='both',shrink=.6)
        fig.savefig(self.out/'plots/foil_delta.png',dpi=140);plt.close(fig)


def run_study(out, a, x, y, folds, cells, cfg, tol, rcond, record, terms, qa_min):
    from elastic_convergence import study
    pool = Pool(out,a,y,folds,cells,terms,qa_min)
    baseline = {}
    for fold in sorted(set(folds)):
        pool.fold = int(fold)
        print(f'Validation fold {fold+1}/{len(set(folds))}',flush=True)
        baseline[str(fold)] = study(x,y,folds,cells,dict(cfg,fold=int(fold)),tol,rcond,record,pool)
    complete = pool.finish(cfg)
    return dict(baselines=baseline,converged_cases=len(pool.rows)//2,
                cases=len(set(folds))*len(TARGETS)*len(cfg['alphas'])*len(cfg['l1']),pooled_cases=complete)
