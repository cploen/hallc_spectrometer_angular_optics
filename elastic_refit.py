"""Unpenalized selected-term refits and validation diagnostics for one saved fold."""
import numpy as np
from core_sample import write_tsv
from elastic_net import TARGETS, SCALES, macro_mse
from elastic_diagnostics import groups, residuals

MODELS = ('full_svd', 'enet', 'refit_svd')


def subset_svd(x, y, selected, rcond):
    """Solve the selected columns of scaled X directly; y is centered on fitting data."""
    beta = np.zeros(x.shape[1])
    if not np.any(selected): return beta, 0, np.array([])
    values, _, rank, singular = np.linalg.lstsq(x[:, selected], y, rcond=rcond)
    beta[selected] = values
    if not np.isfinite(beta).all(): raise ValueError('Nonfinite selected-term SVD coefficients')
    return beta, int(rank), singular


def native(beta, state, j):
    mean, scale, ym, ys, constant = state
    b = beta*ys[j]/scale
    b[constant] = 0.
    return np.r_[ym[j]-mean @ b, b]


class Comparison:
    def __init__(self, out, data, truth, offset, terms, qa_min, spec):
        self.out, self.data, self.truth, self.offset = out, data, truth, offset
        self.terms, self.qa_min, self.spec = terms, qa_min, spec
        self.summary, self.details, self.term_rows = [], [], []
        self.cases, self.coefficients = {}, {}
        (out/'plots').mkdir()

    def __call__(self, j, alpha, l1, enet, fit):
        selected = enet != 0
        refit, rank, singular = subset_svd(fit['x'], fit['y'][:, j], selected, fit['rcond'])
        # Refit coefficients may differ; selected terms never come from validation.
        beta = dict(full_svd=fit['svd'][:, j], enet=enet, refit_svd=refit)
        mask, state = fit['mask'], fit['state']
        a = {k: v[mask] for k, v in self.data.items()}
        y = self.truth[mask, j]
        predictions = {m: (fit['validation'] @ b)*state[3][j]+state[2][j] for m, b in beta.items()}
        cells = np.array([f'{z}:{d}' for z,d in zip(a['ztarT'], a['ndel'])])
        errors = {m: float(macro_mse(p-y, cells)) for m,p in predictions.items()}
        scores = dict(target=TARGETS[j], alpha=alpha, l1=l1, terms=int(sum(selected))+1,
                      rank=rank, slopes=int(sum(selected)), svd_cutoff=fit['rcond'],
                      condition=float(singular[0]/singular[-1]) if len(singular) and singular[-1] > 0 else '',
                      full_mse=errors['full_svd'], enet_mse=errors['enet'], refit_mse=errors['refit_svd'],
                      enet_ratio=errors['enet']/errors['full_svd'] if errors['full_svd'] else '',
                      refit_ratio=errors['refit_svd']/errors['full_svd'] if errors['full_svd'] else '')
        self.summary.append(scores)
        for g, h in groups(a, np.full(len(y), 'validation')):
            for m,p in predictions.items():
                v=(p[h]-y[h])*SCALES[j]
                self.details.append(dict(g, target=TARGETS[j], alpha=alpha, l1=l1, model=m,
                                         n=int(sum(h)), qa_low=int(sum(h)) < self.qa_min,
                                         bias=float(v.mean()), rms=float(np.sqrt(np.mean(v*v))),
                                         p95=float(np.percentile(np.abs(v),95))))
        coeff = {m: native(b, state, j) for m,b in beta.items()}
        for i,t in enumerate(self.terms):
            self.term_rows.append(dict(target=TARGETS[j], alpha=alpha, l1=l1, term=''.join(map(str,t)),
                                      selected=bool(i == 0 or selected[i-1]),
                                      **{m: float(b[i]) for m,b in coeff.items()}))
        key = (alpha,l1)
        bundle = self.cases.setdefault(key, {})
        bundle[j] = predictions
        # Transport-unit coefficients include the unpenalized constant. Fixed seed terms are separate.
        for m,b in coeff.items(): self.coefficients[f'{TARGETS[j]}_a{alpha:g}_l{l1:g}_{m}'] = b
        write_tsv(self.out/'comparison.tsv',self.summary)
        write_tsv(self.out/'validation.tsv',self.details)
        write_tsv(self.out/'terms.tsv',self.term_rows)
        np.savez_compressed(self.out/'coefficients.npz', terms=self.terms, **self.coefficients)
        print(f'  Selected-term SVD: {int(sum(selected))+1} terms including constant; '
              f'validation MSE/full SVD={scores["refit_ratio"]:.4g}',flush=True)
        self.validation = a
        self.validation_mask = mask

    def finish(self):
        if not self.summary:
            return 'No cases converged; no term sets were refitted. Inspect checkpoints.tsv.\n'
        # ztar needs paired ytar/yptar predictions. Report complete three-target cases only.
        for (alpha,l1), bundle in self.cases.items():
            if set(bundle) != {0,1,2}: continue
            a = self.validation
            truth, offset = self.truth[self.validation_mask], self.offset[self.validation_mask]
            preds = {m: np.column_stack([bundle[j][m] for j in range(3)]) for m in MODELS}
            rr = {m: residuals(a, truth, p, offset)[:,3] for m,p in preds.items()}
            for g,h in groups(a,np.full(len(truth),'validation')):
                for m,r in rr.items():
                    v=r[h]
                    self.details.append(dict(g,target='ztar',alpha=alpha,l1=l1,model=m,n=int(sum(h)),
                                             qa_low=int(sum(h)) < self.qa_min,bias=float(v.mean()),
                                             rms=float(np.sqrt(np.mean(v*v))),p95=float(np.percentile(np.abs(v),95))))
        write_tsv(self.out/'validation.tsv',self.details)
        self.plots()
        text = '## Unpenalized selected-term SVD comparison\n\n'
        text += '| Target | Alpha | L1 | Terms including constant | Penalized MSE / full SVD | Refit MSE / full SVD |\n|---|---:|---:|---:|---:|---:|\n'
        for r in self.summary:
            text += f"| {r['target']} | {r['alpha']:g} | {r['l1']:g} | {r['terms']} | {r['enet_ratio']:.4g} | {r['refit_ratio']:.4g} |\n"
        text += '\nBelow 1 means lower validation MSE than full-basis SVD. Only converged elastic-net cases are refitted. '
        text += 'Term selection, centering, scaling and both SVD solves use only the fitting portion of the fold. '
        text += 'The constant stays free of penalties; fixed xtar terms and delta coefficients are retained in seed.dat. '
        text += 'A score near or below 1 is a candidate for the three-fold study, not a final matrix choice.\n\n'
        text += 'comparison.tsv records terms, rank and scores; validation.tsv records N, bias, RMS and P95 '
        text += 'by physical foil/delta, setting and setting/hole. ztar is included only where all three target fits '
        text += 'converged for the same penalty settings. terms.tsv and coefficients.npz preserve the fitting-fold '
        text += 'selection and native coefficients. These are development fits, not full-training replay matrices.\n\n'
        text += 'See plots/comparison.png and each case foil/delta and hole map. Error ratios use common scales; '
        text += 'black outlines flag sparse validation holes, and blanks mean no validation events. '
        text += 'Counts are not training allocation counts. Sparse populations may be absent from this fold. '
        text += 'Complete all folds before choosing the basis; keep the protected pools closed.\n'
        return text

    def plots(self):
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm
        # Geometry uses the same profile as the TFit adapter; there is no campaign foil-count assumption.
        spec = self.spec
        fig, axes = plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
        for j,ax in enumerate(axes):
            rows=[r for r in self.summary if r['target']==TARGETS[j]]
            for i,r in enumerate(rows):
                ax.plot([i,i],[r['enet_ratio'],r['refit_ratio']],color='gray')
            ax.scatter(range(len(rows)),[r['enet_ratio'] for r in rows],label='Penalized elastic net',marker='x')
            ax.scatter(range(len(rows)),[r['refit_ratio'] for r in rows],label='Selected-term SVD',marker='o')
            ax.axhline(1,color='black',ls='--',label='Full-basis SVD')
            ax.set(xticks=range(len(rows)),xticklabels=[f"α={r['alpha']:g}\nL1={r['l1']:g}\n{r['terms']} terms" for r in rows],
                   title=TARGETS[j],ylabel='Validation MSE / full-basis SVD')
            ax.tick_params(axis='x',labelsize=8);ax.legend(fontsize=8)
        fig.savefig(self.out/'plots/comparison.png',dpi=140,bbox_inches='tight');plt.close(fig)
        norm=TwoSlopeNorm(vmin=.5,vcenter=1,vmax=1.5)
        foils, deltas=np.unique(self.validation['ztarT']),np.unique(self.validation['ndel'])
        for alpha,l1 in self.cases:
            case=f'a{alpha:g}_l{l1:g}'
            rows=[r for r in self.details if r['alpha']==alpha and r['l1']==l1]
            available=[t for t in (*TARGETS,'ztar') if any(r['target']==t for r in rows)]
            fig,axes=plt.subplots(2,len(available),figsize=(4*len(available),7),squeeze=False,constrained_layout=True)
            for j,t in enumerate(available):
                subset=[r for r in rows if r['level']=='foil_delta' and r['target']==t]
                lookup={(r['model'],float(r['zfoil']),int(r['ndel'])):r for r in subset}
                for i,m in enumerate(('enet','refit_svd')):
                    values=np.full((len(foils),len(deltas)),np.nan)
                    for zi,z in enumerate(foils):
                        for di,d in enumerate(deltas):
                            r=lookup.get((m,z,d));base=lookup.get(('full_svd',z,d))
                            if r and base and base['rms']>0:
                                values[zi,di]=r['rms']/base['rms']
                                axes[i,j].text(di,zi,f"{values[zi,di]:.2f}\nN={r['n']}",ha='center',va='center',fontsize=7)
                    im=axes[i,j].imshow(values,norm=norm,cmap='coolwarm',aspect='auto')
                    axes[i,j].set(title=f'{t}: {m}',xticks=range(len(deltas)),xticklabels=deltas,
                                  yticks=range(len(foils)),yticklabels=foils,xlabel='Delta slice',ylabel='Foil (cm)')
            fig.colorbar(im,ax=axes.ravel().tolist(),label='Validation RMS / full-basis SVD',extend='both',shrink=.6)
            fig.suptitle(f'α={alpha:g}, L1={l1:g}: saved validation fold')
            fig.savefig(self.out/f'plots/{case}_foil_delta.png',dpi=140,bbox_inches='tight');plt.close(fig)
            # Pool settings only for the map; exact setting/hole rows remain in the table.
            for target in (t for t in available if t != 'ztar'):
                r=[v for v in rows if v['level']=='hole' and v['target']==target]
                pooled={}
                for v in r:
                    key=(v['model'],v['zfoil'],v['ndel'],v['xscol'],v['yscol'])
                    n,s=pooled.get(key,(0,0.))
                    pooled[key]=(n+v['n'],s+v['n']*v['rms']**2)
                fig,axes=plt.subplots(len(foils),len(deltas),figsize=(3.4*len(deltas),3.1*len(foils)),squeeze=False,constrained_layout=True)
                for zi,z in enumerate(foils):
                    for di,d in enumerate(deltas):
                        ax=axes[zi,di]
                        for key,(n,s) in pooled.items():
                            m,zf,nd,xh,yh=key
                            if m!='refit_svd' or zf!=z or nd!=d:continue
                            bn,bs=pooled[('full_svd',zf,nd,xh,yh)]
                            if bs<=0:continue
                            ratio=np.sqrt((s/n)/(bs/bn))
                            xx,yy=spec.ys(yh),spec.xs(xh)
                            ax.scatter(xx,yy,c=[ratio],norm=norm,cmap='coolwarm',s=220,edgecolors='black' if n<self.qa_min else 'none')
                            ax.text(xx,yy,f'{ratio:.2f}',ha='center',va='center',fontsize=6,
                                    color='white' if abs(ratio-1)>.3 else 'black')
                        ax.set(xlim=(spec.ys(0)-1,spec.ys(spec.ny-1)+1),ylim=(spec.xs(0)-1,spec.xs(spec.nx-1)+1),
                               aspect='equal',xlabel='Y sieve (cm)',ylabel='X sieve (cm)',title=f'Foil {z:g}, delta {d}')
                fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='coolwarm'),ax=axes.ravel().tolist(),label='Refit / full SVD RMS',extend='both',shrink=.5)
                fig.suptitle(f'{target}, α={alpha:g}, L1={l1:g}: validation hole RMS ratios\nBlack outline: N<{self.qa_min}; blank: no validation estimate; N is in validation.tsv')
                fig.savefig(self.out/f'plots/{case}_{target}_holes.png',dpi=140,bbox_inches='tight');plt.close(fig)
