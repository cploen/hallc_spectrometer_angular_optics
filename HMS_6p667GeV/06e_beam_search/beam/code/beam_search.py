"""Deterministic add/remove beam search with direct-X SVD and pooled cell scoring."""
import time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from threadpoolctl import threadpool_limits
from elastic_net import TARGETS, SCALES, scaling, transform
from elastic_refit import native

DEFAULTS = dict(beam=4, threads=8, steps=20, patience=3, gain=.001, slack=.005)


def validate(cfg):
    if set(cfg) != set(DEFAULTS): raise ValueError('Unknown/missing beam configuration keys')
    for k in ('beam','threads','steps','patience'):
        if type(cfg[k]) is not int or cfg[k] < 1: raise ValueError(f'Invalid {k}')
    for k in ('gain','slack'):
        if type(cfg[k]) not in (int,float) or not np.isfinite(cfg[k]) or not 0 <= cfg[k] < 1:
            raise ValueError(f'Invalid {k}')


def solve(r, y, selected, rcond):
    """R is an orthogonal compression of X, not the normal-equation matrix."""
    if not selected: return np.zeros(0),0,np.array([])
    beta,_,rank,s = np.linalg.lstsq(r[:,selected],y,rcond=rcond)
    if not np.isfinite(beta).all(): raise ValueError('Nonfinite candidate coefficients')
    return beta,int(rank),s


class Folds:
    def __init__(self,x,y,folds,cells,rcond,threads=8):
        self.p,self.n,self.rcond = x.shape[1],len(x),rcond
        self.x,self.y,self.folds = x,y,folds
        self.parts = []
        _,inverse,counts = np.unique(cells,return_inverse=True,return_counts=True)
        # Across all folds, these weights give exactly the pooled mean cell MSE.
        weights = 1/np.sqrt(len(counts)*counts[inverse])
        with threadpool_limits(limits=threads):
            for f in sorted(set(folds)):
                tr,va = folds != f,folds == f
                state = scaling(x[tr],y[tr])
                z,t = transform(x[tr],y[tr],state)
                zv,tv = transform(x[va],y[va],state)
                fit = np.linalg.qr(np.column_stack((z,t)),mode='r')
                score = np.linalg.qr(np.column_stack((zv,tv))*weights[va,None],mode='r')
                self.parts.append(dict(fold=int(f),mask=va,state=state,fit=fit,score=score))
                print(f'Prepared fold {f+1}: {sum(tr):,} fit, {sum(va):,} validation events',flush=True)

    def evaluate(self,j,selected):
        score,conditions,ranks = 0.,[],[]
        for f in self.parts:
            beta,rank,s = solve(f['fit'][:,:self.p],f['fit'][:,self.p+j],selected,self.rcond)
            r = f['score'][:,selected] @ beta-f['score'][:,self.p+j]
            score += float(r @ r)*(f['state'][3][j]*SCALES[j])**2
            conditions.append(float(s[0]/s[-1]) if len(s) and s[-1] > 0 else 1. if not selected else np.inf)
            ranks.append(rank)
        return dict(mask=selected,terms=len(selected)+1,score=score,condition=max(conditions),
                    rank_min=min(ranks),valid=bool(all(r == len(selected) for r in ranks) and np.isfinite(score)))

    def predict(self,j,selected):
        predicted = np.empty(self.n)
        rows = []
        for f in self.parts:
            beta,rank,s = solve(f['fit'][:,:self.p],f['fit'][:,self.p+j],selected,self.rcond)
            full = np.zeros(self.p); full[list(selected)] = beta
            z,_ = transform(self.x[f['mask']],self.y[f['mask']],f['state'])
            predicted[f['mask']] = z @ full*f['state'][3][j]+f['state'][2][j]
            rows.append(dict(fold=f['fold'],terms=len(selected)+1,rank=rank,
                             condition=float(s[0]/s[-1]) if len(s) and s[-1] > 0 else 1. if not selected else np.inf))
        return predicted,rows


def neighbours(mask,p):
    s = set(mask)
    return [tuple(sorted(s ^ {i})) for i in range(p)]


def search(evaluate,seeds,p,cfg,record):
    """Score every one-term add/remove neighbour; keep a bounded next frontier."""
    validate(cfg)
    cache,history = {},[]
    clock = time.monotonic()
    order = lambda r:(r['score'],r['terms'],r['condition'],r['mask'])
    with threadpool_limits(limits=1),ThreadPoolExecutor(max_workers=cfg['threads']) as workers:
        def assess(masks,step):
            todo = sorted(set(masks)-set(cache))
            for r in workers.map(evaluate,todo):
                cache[r['mask']] = r
                record(dict(r,step=step))
            return [cache[m] for m in todo if cache[m]['valid']]
        initial = assess(seeds,0)
        if not initial: raise ValueError('No full-rank seed basis in every fold')
        start = min(initial,key=order)
        best = start
        frontier = sorted(initial,key=order)[:cfg['beam']]
        anchor,stale = best['score'],0
        reason = 'step limit'
        for step in range(cfg['steps']+1):
            if step:
                candidates = assess([m for r in frontier for m in neighbours(r['mask'],p)],step)
                if not candidates:
                    reason = 'no unvisited candidates';break
                frontier = sorted(candidates,key=order)[:cfg['beam']]
                best = min([best,*frontier],key=order)
                if best['score'] < anchor*(1-cfg['gain']): anchor,stale = best['score'],0
                else: stale += 1
            history.append(dict(step=step,score=best['score'],terms=best['terms'],condition=best['condition'],
                                evaluated=len(cache),seconds=time.monotonic()-clock))
            print(f'  Step {step}: {len(cache):,} candidates, best MSE={best["score"]:.6g}, '
                  f'{best["terms"]} terms, condition={best["condition"]:.3g}',flush=True)
            if stale >= cfg['patience']:
                reason = 'no material gain';break
    # A small, explicit residual allowance keeps negligible gains from requiring extra terms.
    eligible = [r for r in cache.values() if r['valid'] and r['score'] <= best['score']*(1+cfg['slack'])]
    chosen = min(eligible,key=lambda r:(r['terms'],r['score'],r['condition'],r['mask']))
    return dict(start=start,best=best,chosen=chosen,history=history,stop=reason,evaluated=len(cache))


def fit_all(x,y,masks,rcond):
    """Final coefficient fit uses all training events and the selected fixed bases."""
    state = scaling(x,y)
    z,t = transform(x,y,state)
    coeff,info = [],[]
    for j,selected in enumerate(masks):
        beta,rank,s = solve(z,t[:,j],selected,rcond)
        if rank != len(selected): raise ValueError(f'Final {TARGETS[j]} basis is rank deficient')
        full = np.zeros(x.shape[1]);full[list(selected)] = beta
        coeff.append(native(full,state,j))
        info.append(dict(target=TARGETS[j],terms=len(selected)+1,rank=rank,
                         condition=float(s[0]/s[-1]) if len(s) else 1.))
    return np.column_stack(coeff),info
