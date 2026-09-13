"""Scaled elastic-net angular fit; no campaign geometry or event selection here."""
import hashlib
import warnings
import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet

TARGETS = ('xptar', 'ytar', 'yptar')
SCALES = np.array([1000., 100., 1000.])
DEFAULTS = dict(seed=667, folds=3, alphas=[.1, .03, .01, .003, .001, .0003, .0001, .00003, .00001],
                l1=[.1, .5, .9], tol=1e-6, max_iter=20000, qa_min=10, rcond=1e-12)


def validate(cfg):
    if set(cfg) != set(DEFAULTS):
        raise ValueError('Unknown/missing elastic configuration keys')
    for k in ('seed', 'folds', 'max_iter', 'qa_min'):
        if type(cfg[k]) is not int or cfg[k] < (2 if k == 'folds' else 0 if k == 'seed' else 1):
            raise ValueError(f'Invalid {k}')
    for k in ('alphas', 'l1'):
        a = cfg[k]
        if not isinstance(a, list) or not a or any(type(v) not in (float, int) or not np.isfinite(v) or v <= 0 for v in a):
            raise ValueError(f'Invalid {k}')
        if len(set(a)) != len(a) or (k == 'l1' and max(a) > 1):
            raise ValueError(f'Invalid {k}')
    for k in ('tol', 'rcond'):
        if type(cfg[k]) not in (float, int) or not 0 < cfg[k] < 1:
            raise ValueError(f'Invalid {k}')


def folds_for(ids, strata, nfold, seed):
    """Round-robin within each labeled leaf; hash order and offset are stable."""
    out = np.empty(len(ids), dtype=int)
    grouped = {}
    for i, s in enumerate(strata): grouped.setdefault(s, []).append(i)
    for s, idx in sorted(grouped.items()):
        key = lambda i: hashlib.sha256(f'{seed}:{ids[i]}'.encode()).digest()
        idx.sort(key=key)
        offset = int.from_bytes(hashlib.sha256(f'{seed}:{s}'.encode()).digest()[:8], 'big') % nfold
        out[idx] = (np.arange(len(idx)) + offset) % nfold
    if set(out) != set(range(nfold)):
        raise ValueError('Too few training events for all validation folds')
    return out


def scaling(x, y):
    mean, scale = x.mean(axis=0), x.std(axis=0)
    # An exactly constant column conveys no slope information; preserve an audit.
    constant = scale <= np.finfo(float).eps * np.maximum(np.abs(mean), np.finfo(float).tiny)
    scale[constant] = 1.
    ym, ys = y.mean(axis=0), y.std(axis=0)
    if np.any(ys <= 0):
        raise ValueError('Training target has zero variance')
    return mean, scale, ym, ys, constant


def transform(x, y, state):
    xm, xs, ym, ys, constant = state
    z = (x-xm)/xs
    z[:, constant] = 0.
    return np.asfortranarray(z), (y-ym)/ys


def raw_coefficients(beta, state):
    xm, xs, ym, ys, constant = state
    slopes = beta * ys / xs[:, None]
    slopes[constant] = 0.
    return np.vstack((ym-xm @ slopes, slopes))


def macro_mse(residual, cells):
    return np.mean([np.mean(residual[cells == c]**2, axis=0) for c in np.unique(cells)], axis=0)


def fit_path(x, y, cfg):
    """Yield independent-target fits, warming only within one l1 path."""
    gram = np.ascontiguousarray(x.T @ x)
    for l1 in cfg['l1']:
        models = [ElasticNet(l1_ratio=l1, fit_intercept=False, precompute=gram,
                             warm_start=True, max_iter=cfg['max_iter'], tol=cfg['tol']) for _ in TARGETS]
        for alpha in sorted(cfg['alphas'], reverse=True):
            beta, ok, gaps = [], [], []
            for j, model in enumerate(models):
                model.set_params(alpha=alpha)
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always', ConvergenceWarning)
                    model.fit(x, y[:, j])
                beta.append(model.coef_.copy())
                good = not any(issubclass(w.category, ConvergenceWarning) for w in caught)
                ok.append(good and np.isfinite(model.coef_).all() and np.isfinite(model.dual_gap_))
                gaps.append(float(model.dual_gap_))
            yield l1, alpha, np.array(beta).T, ok, gaps


def tune(x, y, folds, cells, cfg):
    records, supports = [], {}
    baseline = []
    for f in range(cfg['folds']):
        tr, va = folds != f, folds == f
        state = scaling(x[tr], y[tr])
        z, t = transform(x[tr], y[tr], state)
        design = np.column_stack((np.ones(sum(va)), x[va]))
        ls, _, rank, singular = np.linalg.lstsq(z, t, rcond=cfg['rcond'])
        residual = design @ raw_coefficients(ls, state)-y[va]
        baseline.append(macro_mse(residual, cells[va]))
        print(f'Elastic validation fold {f+1}/{cfg["folds"]}: {sum(tr)} fit, {sum(va)} validation', flush=True)
        for l1, alpha, beta, ok, gaps in fit_path(z, t, cfg):
            errors = macro_mse(design @ raw_coefficients(beta, state)-y[va], cells[va])
            for j, name in enumerate(TARGETS):
                records.append(dict(target=name, fold=f, l1=l1, alpha=alpha, mse=float(errors[j]),
                                    active=int(np.count_nonzero(beta[:, j])), converged=ok[j], dual_gap=gaps[j]))
                supports[(j, f, l1, alpha)] = beta[:, j] != 0
    choices, summary, frequency = [], [], []
    for j, name in enumerate(TARGETS):
        candidates = []
        for l1 in cfg['l1']:
            for alpha in cfg['alphas']:
                r = [v for v in records if v['target'] == name and v['l1'] == l1 and v['alpha'] == alpha]
                mse = np.array([v['mse'] for v in r])
                row = dict(target=name, l1=l1, alpha=alpha, mse=float(mse.mean()),
                           se=float(mse.std(ddof=1)/np.sqrt(len(mse))), active=float(np.mean([v['active'] for v in r])),
                           converged=all(v['converged'] for v in r))
                summary.append(row)
                if row['converged']: candidates.append(row)
        if not candidates:
            raise ValueError(f'No converged candidate in every fold for {name}; increase max_iter or regularization')
        best = min(candidates, key=lambda r: r['mse'])
        # One-SE rule: prefer fewer terms among statistically similar CV scores.
        chosen = min([r for r in candidates if r['mse'] <= best['mse']+best['se']],
                     key=lambda r: (r['active'], -r['alpha'], -r['l1']))
        frequency.append(np.mean([supports[j, f, chosen['l1'], chosen['alpha']] for f in range(cfg['folds'])], axis=0))
        choices.append(dict(chosen, best_mse=best['mse'], svd_cv_mse=float(np.mean(baseline, axis=0)[j]),
                            alpha_edge=chosen['alpha'] in (min(cfg['alphas']), max(cfg['alphas']))))
    state = scaling(x, y)
    z, t = transform(x, y, state)
    beta = np.zeros((x.shape[1], 3))
    for l1, alpha, b, ok, gaps in fit_path(z, t, cfg):
        for j, c in enumerate(choices):
            if (l1, alpha) == (c['l1'], c['alpha']):
                if not ok[j]: raise ValueError(f'Final {TARGETS[j]} fit did not converge')
                beta[:, j] = b[:, j]
                c['final_active'] = int(np.count_nonzero(b[:, j]))
                c['final_dual_gap'] = gaps[j]
    ls, _, rank, singular = np.linalg.lstsq(z, t, rcond=cfg['rcond'])
    info = dict(rank=int(rank), columns=x.shape[1], rcond=cfg['rcond'],
                constant_columns=np.flatnonzero(state[4]).tolist(), singular_values=singular.tolist())
    return dict(enet=raw_coefficients(beta, state), svd=raw_coefficients(ls, state),
                state=state, beta=beta, frequency=np.array(frequency).T, choices=choices,
                cv=records, summary=summary, svd_info=info)
