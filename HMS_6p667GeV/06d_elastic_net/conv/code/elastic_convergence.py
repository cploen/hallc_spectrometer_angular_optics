"""Bounded continuation study on one saved training-validation fold."""
import importlib.metadata
import json
import shutil
import time
import warnings
from pathlib import Path
import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet
from core_sample import PROJECT, digest, read_tsv, write_tsv
from elastic_net import TARGETS, SCALES, scaling, transform, macro_mse

DEFAULTS = dict(fold=0, alphas=[.0001, .00003], l1=[.1, .5], budgets=[20000, 100000, 300000])


def validate(cfg, nfold):
    if set(cfg) != set(DEFAULTS): raise ValueError('Unknown/missing convergence configuration keys')
    if type(cfg['fold']) is not int or not 0 <= cfg['fold'] < nfold: raise ValueError('Invalid fold')
    for k in ('alphas', 'l1'):
        values = cfg[k]
        if not isinstance(values, list) or not values or any(type(v) not in (int, float) or not np.isfinite(v) or v <= 0 for v in values):
            raise ValueError(f'Invalid {k}')
        if len(set(values)) != len(values) or (k == 'l1' and max(values) > 1): raise ValueError(f'Invalid {k}')
    b = cfg['budgets']
    if not isinstance(b, list) or not b or any(type(v) is not int or v <= 0 for v in b) or b != sorted(set(b)):
        raise ValueError('Budgets must be strictly increasing positive integers')


def saved_folds(a, rows, nfold):
    lookup = {(r['rungroup'], int(r['entry'])): int(r['fold']) for r in rows}
    keys = list(zip(a['rungroup'], map(int, a['entry'])))
    if len(lookup) != len(rows) or len(set(keys)) != len(keys) or set(keys) != set(lookup):
        raise ValueError('Saved fold IDs differ from current training IDs')
    folds = np.array([lookup[k] for k in keys])
    if set(folds) != set(range(nfold)): raise ValueError('Invalid saved fold numbers')
    return folds


def checkpoints(x, y, gram, alpha, l1, budgets, tol):
    """Each case starts at zero; later budgets continue that same fixed objective."""
    model = ElasticNet(alpha=alpha, l1_ratio=l1, fit_intercept=False, precompute=gram,
                       warm_start=True, tol=tol, selection='cyclic')
    used, elapsed = 0, 0.
    xy = x.T @ y / len(y)
    for budget in budgets:
        model.set_params(max_iter=budget-used)
        start = time.monotonic()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always', ConvergenceWarning)
            model.fit(x, y)
        seconds = time.monotonic()-start
        elapsed += seconds
        used += int(model.n_iter_)
        beta = model.coef_.copy()
        if not np.isfinite(beta).all() or not np.isfinite(model.dual_gap_):
            raise ValueError('Nonfinite solver result')
        residual = x @ beta-y
        objective = .5*np.mean(residual**2)+alpha*l1*np.abs(beta).sum()+.5*alpha*(1-l1)*(beta @ beta)
        gradient = gram @ beta/len(y)-xy+alpha*(1-l1)*beta
        kkt = np.where(beta != 0, np.abs(gradient+alpha*l1*np.sign(beta)), np.maximum(np.abs(gradient)-alpha*l1, 0))
        warned = any(issubclass(w.category, ConvergenceWarning) for w in caught)
        # Keep the original sklearn criterion; the recorded gap check is explicit too.
        limit = tol*float(y @ y)/len(y)
        converged = not warned and float(model.dual_gap_) <= limit
        yield beta, dict(budget=budget, iterations=used, seconds=elapsed, converged=converged,
                         dual_gap=float(model.dual_gap_), gap_limit=limit, gap_ratio=float(model.dual_gap_)/limit,
                         objective=float(objective), kkt_max=float(kkt.max()), active=int(np.count_nonzero(beta)))
        if converged: break


def study(x, y, folds, cells, cfg, tol, rcond, record):
    tr, va = folds != cfg['fold'], folds == cfg['fold']
    state = scaling(x[tr], y[tr])
    z, t = transform(x[tr], y[tr], state)
    zv, _ = transform(x[va], y[va], state)
    gram = np.ascontiguousarray(z.T @ z)
    ls, _, rank, singular = np.linalg.lstsq(z, t, rcond=rcond)
    baseline = macro_mse(zv @ ls*state[3]+state[2]-y[va], cells[va])
    for j, target in enumerate(TARGETS):
        for l1 in cfg['l1']:
            for alpha in sorted(cfg['alphas'], reverse=True):
                print(f'Study {target}: L1={l1:g}, alpha={alpha:g}', flush=True)
                previous = None
                for beta, info in checkpoints(z, t[:, j], gram, alpha, l1, cfg['budgets'], tol):
                    prediction = (zv @ beta)*state[3][j]+state[2][j]
                    mse = float(macro_mse(prediction-y[va, j], cells[va]))
                    change = '' if previous is None else float(np.sqrt(np.mean((prediction-previous)**2))*SCALES[j])
                    row = dict(target=target, fold=cfg['fold'], alpha=alpha, l1=l1, **info,
                               validation_mse=mse, svd_mse=float(baseline[j]),
                               mse_ratio=mse/baseline[j] if baseline[j] > 0 else '',
                               validation_rms=float(np.sqrt(mse)*SCALES[j]),
                               prediction_change_rms=change)
                    record(row)
                    previous = prediction
                    print(f'  {info["iterations"]:,} iterations; gap/tolerance={info["gap_ratio"]:.3g}; '
                          f'MSE/SVD={row["mse_ratio"]:.4g}; converged={info["converged"]}', flush=True)
    return dict(training=int(sum(tr)), validation=int(sum(va)), rank=int(rank),
                singular_values=singular.tolist(), svd_mse=baseline.tolist())


def plot(out, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for j, target in enumerate(TARGETS):
        cases = sorted({(r['l1'], r['alpha']) for r in rows if r['target'] == target})
        for l1, alpha in cases:
            r = [v for v in rows if v['target'] == target and v['l1'] == l1 and v['alpha'] == alpha]
            label = f'L1 {l1:g}, α {alpha:g}'
            line, = axes[0,j].plot([v['iterations'] for v in r], [max(v['gap_ratio'],1e-12) for v in r], label=label)
            for ax, key in ((axes[0,j], 'gap_ratio'), (axes[1,j], 'mse_ratio')):
                if key == 'mse_ratio': ax.plot([v['iterations'] for v in r], [v[key] for v in r], color=line.get_color())
                for v in r:
                    ax.scatter(v['iterations'], max(v[key],1e-12), color=line.get_color(),
                               marker='o' if v['converged'] else 'x', s=40)
        for ax in axes[:,j]: ax.set_xscale('log'); ax.axhline(1,color='black',ls='--',lw=1)
        axes[0,j].set(yscale='log', title=target, ylabel='Solver gap / tolerance limit')
        axes[1,j].set(xlabel='Actual cumulative iterations', ylabel='Validation MSE / scaled SVD')
        axes[0,j].legend(fontsize=8)
    fig.suptitle('One saved validation fold — circles: converged; crosses: unfinished\nTolerance is unchanged; unfinished scores cannot select a final model')
    fig.savefig(out/'convergence.png',dpi=140,bbox_inches='tight'); plt.close(fig)


def run(campaign, tag, name, source='enet'):
    from fit_elastic import checked, load_sample, matrix_rows, problem
    parent = campaign/'06d_elastic_net'/source
    out = campaign/'06d_elastic_net'/name
    if out.exists(): raise FileExistsError(f'Output exists: {out}')
    base = json.loads((parent/'manifest.json').read_text())
    if base.get('schema') != 'elastic_v1' or base['sample'] != tag: raise ValueError('Source fit/tag mismatch')
    hashes = {}
    for rel, sha in base['outputs'].items(): checked(parent/rel, sha, hashes)
    for path, sha in base['inputs'].items(): checked(Path(path), sha, hashes)
    a, inputs = load_sample(campaign, tag, 'fit')
    folds = saved_folds(a, read_tsv(parent/'tsv/folds.tsv'), base['config']['folds'])
    cfg = dict(DEFAULTS)
    policy = campaign/'config/elastic_conv.json'
    if policy.exists(): cfg.update(json.loads(policy.read_text()))
    validate(cfg, base['config']['folds'])
    x, y, _ = problem(a, matrix_rows(parent/'seed.dat'))
    cells = np.array([f'{z}:{d}' for z,d in zip(a['ztarT'],a['ndel'])])
    out.mkdir(parents=True)
    (out/'code').mkdir()
    for file in ('elastic_convergence.py','elastic_net.py','fit_elastic.py','core_sample.py',
                 'preallocated_svd.py','spectrometer_config.py','spectrometer_profiles.def','run_elastic.sh'):
        shutil.copy2(PROJECT/file,out/'code'/file)
    rows = []
    manifest = dict(schema='elastic_convergence_v1', status='running', sample=tag, source=source,
                    source_manifest=digest(parent/'manifest.json'), config=cfg, tol=base['config']['tol'],
                    rcond=base['config']['rcond'], inputs=inputs,
                    versions={m: importlib.metadata.version(m) for m in ('numpy','scipy','scikit-learn','uproot','matplotlib')})
    def save_manifest(): (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    def record(row):
        rows.append(row)
        write_tsv(out/'checkpoints.tsv',rows)
    save_manifest()
    try:
        manifest['baseline'] = study(x[:,1:], y, folds, cells, cfg, manifest['tol'],manifest['rcond'],record)
        plot(out, rows)
        final = {}
        for r in rows: final[r['target'],r['l1'],r['alpha']] = r
        text = '# Elastic-net convergence study\n\n'
        text += f"Saved fold {cfg['fold']} of {source}; tolerance {manifest['tol']:g}, unchanged. Each case starts at zero and continues the same objective between checkpoints.\n\n"
        text += '| Target | L1 | Alpha | Iterations | Converged | Gap / limit | Validation MSE / SVD |\n|---|---:|---:|---:|---|---:|---:|\n'
        for r in final.values():
            text += f"| {r['target']} | {r['l1']:g} | {r['alpha']:g} | {r['iterations']} | {r['converged']} | {r['gap_ratio']:.3g} | {r['mse_ratio']:.4g} |\n"
        text += '\nSee convergence.png and checkpoints.tsv. Crosses remain unfinished even if their validation scores look good. '
        text += 'One-fold scores guide the next experiment; they do not select a final matrix. '
        text += 'No protected events were fitted or evaluated, and no matrix was exported. The SVD reference solves scaled X directly.\n\n'
        text += 'If convergence is reached and scores stabilize, extend the useful settings to all saved folds. '
        text += 'If not, inspect objective, gap, KKT residual, prediction changes and runtime before raising the budget again or changing solver. '
        text += 'Do not loosen tolerance to obtain a passing status. Starts differ from the original descending-alpha warm path, so this is a controlled continuation study rather than an exact replay of its 20,000-iteration fit.\n'
        (out/'RESULTS.md').write_text(text)
        manifest['status'] = 'complete'
        manifest['converged_cases'] = sum(r['converged'] for r in final.values())
        manifest['cases'] = len(final)
        manifest['outputs'] = {str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and p.name != 'manifest.json'}
        save_manifest()
    except BaseException as exc:
        manifest['status'] = 'incomplete'; manifest['error'] = str(exc); save_manifest()
        raise
    print(f'Convergence study complete: {out}; {manifest["converged_cases"]}/{manifest["cases"]} cases converged. No protected evaluation.',flush=True)
