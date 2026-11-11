#!/usr/bin/env python3
"""Conditional x-target angular refit on frozen selections, with geometric xtarT."""
import argparse
import json
from pathlib import Path
import re
import shutil
import numpy as np
from threadpoolctl import threadpool_limits
from core_sample import PROJECT, digest, read_tsv, write_tsv
from fit_elastic import load_sample, matrix_rows, design, checked
from elastic_net import DEFAULTS, TARGETS, SCALES, tune, validate
from elastic_convergence import saved_folds
from beam_search import Folds, search, fit_all, DEFAULTS as BEAM_DEFAULTS


def conditional_problem(a, rows, max_order):
    terms = [t for c, t in rows if t[4] > 0 and sum(t) <= max_order]
    if not terms:
        raise ValueError('No x-target terms within the requested order')
    coordinate = dict(a, xtar=a['xtarT'])
    fixed = [(c, t) for c, t in rows if t not in terms]
    offset = design(coordinate, [t for c, t in fixed]) @ np.array([c[:3] for c, t in fixed])
    truth = np.column_stack((a['xptarT'], a['ytarT']/100., a['yptarT']))
    return terms, design(coordinate, terms), truth-offset, offset


def predict_matrix(a, rows, coordinate):
    data = dict(a, xtar=a[coordinate])
    return design(data, [t for c, t in rows]) @ np.array([c[:3] for c, t in rows])


def export_matrix(source, target, terms, coefficients):
    """Preserve frozen rows byte-for-byte and delta tokens on every free row."""
    lookup = {''.join(map(str, t)): coefficients[i] for i, t in enumerate(terms)}
    lines = []
    for line in source.read_text().splitlines(keepends=True):
        fields = line.split()
        if len(fields) == 5 and fields[-1] in lookup:
            values = lookup[fields[-1]]
            line = ' ' + ' '.join(f'{v:.16e}' for v in values) + ' ' + fields[3] + ' ' + fields[4] + '\n'
        lines.append(line)
    target.write_text(''.join(lines))
    old = {t: c for c, t in matrix_rows(source)}
    new = {t: c for c, t in matrix_rows(target)}
    if old.keys() != new.keys():
        raise ValueError('Export changed the matrix row set')
    for t, c in old.items():
        if new[t][3] != c[3] or (t not in terms and not np.array_equal(new[t], c)):
            raise ValueError('Export changed a frozen coefficient')


def inputs(campaign, tag, baseline, folds_path, max_order):
    a, hashes = load_sample(campaign, tag, 'fit', xtar_coordinate='xtarT')
    hashes[str(baseline)] = digest(baseline)
    hashes[str(folds_path)] = digest(folds_path)
    rows = matrix_rows(baseline)
    for c, t in rows:
        if not any(t) and np.any(c):
            raise ValueError('This zero-offset experiment requires zero matrix constants')
    table = read_tsv(folds_path)
    nfold = len({int(r['fold']) for r in table})
    folds = saved_folds(a, table, nfold)
    terms, x, y, offset = conditional_problem(a, rows, max_order)
    return a, rows, terms, x, y, offset, folds, hashes


def run(args):
    a, rows, terms, x, y, offset, folds, hashes = inputs(
        args.campaign, args.tag, args.baseline, args.folds_from, args.max_order)
    cfg = dict(DEFAULTS)
    policy = args.campaign/'config/elastic.json'
    if policy.exists():
        cfg.update(json.loads(policy.read_text()))
        hashes[str(policy)] = digest(policy)
    cfg['folds'] = len(set(folds))
    validate(cfg)
    beam = dict(BEAM_DEFAULTS, beam=args.beam, steps=args.steps, threads=args.threads)
    print(f'{len(x):,} frozen training events; {len(terms)} free x-target terms through order {args.max_order}', flush=True)
    print('Coordinate=xtarT; constants fixed at zero; non-x-target, higher-order and delta coefficients frozen.', flush=True)
    print('Baseline SHA256:', digest(args.baseline), flush=True)
    if args.check:
        print('Check only; no fit, protected evaluation, or output writes.', flush=True)
        return
    out = args.campaign/'06g_xtar_refit'/args.name
    if out.exists():
        raise FileExistsError('Immutable output exists: ' + str(out))
    for folder in ('matrices', 'tsv', 'plots', 'code'):
        (out/folder).mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.baseline, out/'seed.dat')
    write_tsv(out/'tsv/folds.tsv', [dict(rungroup=r, entry=int(e), fold=int(f))
              for r, e, f in zip(a['rungroup'], a['entry'], folds)])
    codes = [''.join(map(str, t)) for t in terms]
    manifest = dict(schema='xtar_refit_v1', status='running', campaign=args.campaign.name,
        sample=args.tag, coordinate='xtarT', fit_intercept=False, max_order=args.max_order,
        coefficient_mode='full free-block coefficients; penalty toward zero',
        baseline=str(args.baseline), baseline_sha256=digest(args.baseline),
        free_terms=codes, frozen='all non-x-target, higher-order x-target, delta and zero-order coefficients',
        training=len(x), folds=len(set(folds)), enet_config=cfg, beam_config=beam,
        inputs=hashes, holdout_status='not evaluated',
        score='training-only conditional CV; adaptive term selection; frozen other terms previously fitted on training events')
    def save():
        (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    save()
    try:
        for name in ('fit_xtar.py', 'run_xtar_refit.sh', 'fit_elastic.py', 'elastic_net.py',
                     'beam_search.py', 'elastic_refit.py', 'elastic_convergence.py',
                     'spectrometer_config.py', 'spectrometer_profiles.def', 'core_sample.py',
                     'preallocated_svd.py', 'mixed_sample.py', 'elastic_diagnostics.py'):
            shutil.copyfile(PROJECT/name, out/'code'/name)
        cells = np.array([f'{z}:{d}' for z, d in zip(a['ztarT'], a['ndel'])])
        with threadpool_limits(limits=args.threads):
            result = tune(x, y, folds, cells, cfg, fit_intercept=False)
        write_tsv(out/'tsv/enet_cv.tsv', result['cv'])
        write_tsv(out/'tsv/enet_summary.tsv', result['summary'])
        from elastic_diagnostics import training_plots
        training_plots(out, result)
        manifest['enet_choices'] = result['choices']
        singular = result['svd_info']['singular_values']
        manifest['full_basis'] = dict(result['svd_info'], condition=singular[0]/singular[-1] if singular and singular[-1] else None)
        save()
        prepared = Folds(x, y, folds, cells, cfg['rcond'], args.threads, fit_intercept=False)
        results, candidates = [], []
        for j, target in enumerate(TARGETS):
            choice = result['choices'][j]
            seeds = {(), tuple(map(int,np.flatnonzero(result['beta'][:, j])))}
            for f in range(cfg['folds']):
                support = result['fold_supports'][j, f, choice['l1'], choice['alpha']]
                seeds.add(tuple(map(int,np.flatnonzero(support))))
            print('Search x-target block:', target, flush=True)
            def record(row):
                candidates.append(dict(row, target=target, mask=' '.join(codes[i] for i in row['mask'])))
            selected = search(lambda mask: prepared.evaluate(j, mask), sorted(seeds), len(terms), beam, record)
            results.append(selected)
            (out/'search.json').write_text(json.dumps(results, indent=2)+'\n')
        write_tsv(out/'tsv/candidates.tsv', candidates)
        write_tsv(out/'tsv/search.tsv', [dict(target=TARGETS[j], **v) for j, r in enumerate(results) for v in r['history']])
        models = dict(enet=result['enet'][1:], full=result['svd'][1:])
        info = []
        oof = {}
        for model, key in [('start', 'start'), ('beam', 'chosen')]:
            masks = [r[key]['mask'] for r in results]
            with threadpool_limits(limits=args.threads):
                fitted, stats = fit_all(x, y, masks, cfg['rcond'], fit_intercept=False)
            if np.any(fitted[0]):
                raise ValueError('Conditional fit generated an intercept')
            models[model] = fitted[1:]
            info.extend(dict(model=model, **r) for r in stats)
            oof[model] = np.column_stack([prepared.predict(j, m)[0] for j, m in enumerate(masks)])+offset
        for model, coefficients in models.items():
            destination = out/f'matrices/{model}.dat'
            export_matrix(out/'seed.dat', destination, terms, coefficients)
            pred = predict_matrix(a, matrix_rows(destination), 'xtarT')
            if not np.allclose(pred, x@coefficients+offset, rtol=1e-10, atol=1e-11):
                raise ValueError('Export prediction mismatch: ' + model)
        write_tsv(out/'tsv/conditioning.tsv', info)
        write_tsv(out/'tsv/terms.tsv', [dict(target=t, term=code,
            **{m:float(c[i,j]) for m,c in models.items()}) for j,t in enumerate(TARGETS) for i,code in enumerate(codes)])
        baseline_geom = predict_matrix(a, rows, 'xtarT')
        baseline_replay = predict_matrix(a, rows, 'xtar')
        diagnostics(out, a, dict(baseline_geometry=baseline_geom, baseline_replay_coordinate=baseline_replay,
                                start=oof['start'], beam=oof['beam']), 'development_cv')
        for path, expected in hashes.items():
            checked(Path(path), expected, {})
        manifest.update(status='complete', selected_terms={t:r['chosen']['terms'] for t,r in zip(TARGETS,results)})
        manifest['outputs'] = {str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and p.name!='manifest.json'}
        save()
        text = '# Conditional x-target refit\n\n'
        text += 'Training coordinate: xtarT. All constants and the other coefficient blocks are frozen.\n\n'
        text += '| Target | Selected x-target terms | Full-block columns | Final scaled condition |\n|---|---:|---:|---:|\n'
        for j,t in enumerate(TARGETS):
            stat = next(r for r in info if r['model']=='beam' and r['target']==t)
            text += f"| {t} | {results[j]['chosen']['terms']} | {len(terms)} | {stat['condition']:.6g} |\n"
        text += ('\nElastic-net and beam selection use only the frozen training sample and saved folds. '
                 'Development CV is conditional on the already-fitted frozen coefficients and adaptive basis search; '
                 'it is not an independent generalization estimate. Free-block coefficients are fitted in full, '
                 'without an intercept; excluded terms in this block are zero. Higher-order terms retain simulation values. '
                 'The full-block reference may be rank-truncated; inspect manifest full_basis.\n\n'
                 'baseline_geometry evaluates the original matrix at xtarT; baseline_replay_coordinate uses saved xtar. '
                 'Their difference isolates the coordinate substitution. Compare candidate against baseline_geometry '
                 'to assess coefficient changes. No afterburner closure claim or replay installation is made.\n')
        (out/'RESULTS.md').write_text(text)
        manifest['outputs']['RESULTS.md'] = digest(out/'RESULTS.md')
        save()
    except BaseException as exc:
        manifest.update(status='incomplete', error=str(exc)); save(); raise
    print('Completed:', out, flush=True)


def diagnostics(out, a, predictions, label):
    from elastic_diagnostics import metrics
    from spectrometer_config import from_campaign
    truth = np.column_stack((a['xptarT'], a['ytarT']/100., a['yptarT']))
    rr = {m:(p-truth)*SCALES for m,p in predictions.items()}
    pool = np.full(len(truth), label)
    details = metrics(a, truth, predictions, pool, 10, np.zeros_like(truth))
    write_tsv(out/f'tsv/{label}_metrics.tsv', details)
    np.savez_compressed(out/f'{label}_residuals.npz', rungroup=a['rungroup'], entry=a['entry'],
                        xtar=a['xtar'], xtarT=a['xtarT'], **rr)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    foils = np.unique(a['ztarT']); deltas = np.unique(a['ndel'])
    models = list(predictions)
    fig, axes = plt.subplots(len(models), 3, figsize=(18, 3.6*len(models)), squeeze=False, constrained_layout=True)
    for j,t in enumerate(TARGETS):
        values = [r['rms'] for r in details if r['level']=='foil_delta' and r['target']==t]
        norm = Normalize(0, max(values) if values else 1.)
        for i,m in enumerate(models):
            ax = axes[i,j]; grid = np.full((len(foils),len(deltas)), np.nan)
            for zi,z in enumerate(foils):
                for di,d in enumerate(deltas):
                    mask = (a['ztarT']==z)&(a['ndel']==d)
                    if mask.any():
                        value = float(np.sqrt(np.mean(rr[m][mask,j]**2))); grid[zi,di] = value
                        ax.text(di,zi,f'{value:.3f}\nN={sum(mask)}',ha='center',va='center',fontsize=7,
                                color='white' if value < .5*norm.vmax else 'black')
            image = ax.imshow(grid, cmap='viridis',norm=norm,aspect='auto')
            labels = [f"{a['delta_low'][np.flatnonzero(a['ndel']==d)[0]]:g} to {a['delta_high'][np.flatnonzero(a['ndel']==d)[0]]:g}" for d in deltas]
            ax.set(title=f'{m}: {t}',xticks=range(len(deltas)),xticklabels=labels,
                   yticks=range(len(foils)),yticklabels=foils,xlabel='Delta (%)',ylabel='Foil (cm)')
            ax.tick_params(axis='x',labelrotation=45)
        fig.colorbar(image,ax=axes[:,j].tolist(),label='RMS ('+('cm' if t=='ytar' else 'mrad')+')',shrink=.6)
    fig.suptitle(label+': RMS includes bias; shared scale per target; N<10 is sparse, not excluded')
    fig.savefig(out/f'plots/{label}_foil_delta.png',dpi=140); plt.close(fig)
    spec = from_campaign(str(out))
    for lo,hi in [(-8,10),(10,12),(12,14),(-15,-8)]:
        mask = (a['delta']>=lo)&(a['delta']<hi)
        if not mask.any():
            continue
        fig,axes = plt.subplots(1,len(models),figsize=(5*len(models),7),squeeze=False,constrained_layout=True)
        for i,m in enumerate(models):
            ax=axes[0,i]; p=predictions[m]
            xs=a['xtarT']+168*p[:,0]; ys=100*p[:,1]+168*p[:,2]
            ax.hist2d(ys[mask],xs[mask],bins=160,cmin=1)
            ax.set(title=f'{m}\nN={sum(mask):,}',xlabel='Y sieve (cm)',ylabel='X sieve (cm)',aspect='equal')
        fig.suptitle(f'{label}: delta [{lo:g}, {hi:g})%; geometric xtarT projection; not a new replay')
        fig.savefig(out/f'plots/{label}_sieve_{lo:g}_{hi:g}.png',dpi=120); plt.close(fig)


def evaluate(args):
    out = args.campaign/'06g_xtar_refit'/args.name
    manifest = json.loads((out/'manifest.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('Only a completed frozen fit can be evaluated')
    for path, expected in manifest['inputs'].items():
        checked(Path(path), expected, {})
    for path, expected in manifest['outputs'].items():
        checked(out/path, expected, {})
    dest=out/'evaluation'
    if dest.exists():
        raise FileExistsError('Immutable evaluation exists')
    a, hashes = load_sample(args.campaign, args.tag, 'holdout', xtar_coordinate='xtarT')
    training={(r['rungroup'],int(r['entry'])) for r in read_tsv(out/'tsv/folds.tsv')}
    if training.intersection(zip(a['rungroup'],map(int,a['entry']))):
        raise ValueError('Training/holdout overlap')
    for folder in ('plots','tsv'):
        (dest/folder).mkdir(parents=True,exist_ok=True)
    rows=matrix_rows(out/'seed.dat')
    predictions=dict(baseline_geometry=predict_matrix(a,rows,'xtarT'),
                     baseline_replay_coordinate=predict_matrix(a,rows,'xtar'))
    for model in ('full','beam'):
        predictions[model]=predict_matrix(a,matrix_rows(out/f'matrices/{model}.dat'),'xtarT')
    diagnostics(dest,a,predictions,'reserved_holdout')
    counts = dict(all_reserved=len(a['entry']), blocked=int(sum(a['excluded'])))
    for pool, quality in [('core',2),('noncore',1),('unsupported',0)]:
        mask=(a['quality']==quality)&~a['excluded']
        counts[pool]=int(sum(mask))
        if pool=='unsupported' or not mask.any():
            continue
        sub=dest/pool
        for folder in ('plots','tsv'):(sub/folder).mkdir(parents=True,exist_ok=True)
        diagnostics(sub,{k:v[mask] for k,v in a.items()},
                    {k:v[mask] for k,v in predictions.items()},'reserved_'+pool)
    (dest/'manifest.json').write_text(json.dumps(dict(schema='xtar_evaluation_v1',events=len(a['entry']),
        source_manifest_sha256=digest(out/'manifest.json'),inputs=hashes,
        counts=counts,
        status='complete',coordinate='xtarT',interpretation='Reserved event IDs, disjoint from fit; previously inspected holdouts are development evidence, not a fresh independent test',
        outputs={str(p.relative_to(dest)):digest(p) for p in dest.rglob('*') if p.is_file()}),indent=2)+'\n')
    print('Evaluation complete:',dest,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path); p.add_argument('tag'); p.add_argument('name')
    p.add_argument('--baseline',type=Path); p.add_argument('--folds-from',type=Path)
    p.add_argument('--max-order',type=int,default=4); p.add_argument('--beam',type=int,default=30)
    p.add_argument('--steps',type=int,default=40); p.add_argument('--threads',type=int,default=8)
    mode=p.add_mutually_exclusive_group();mode.add_argument('--check',action='store_true');mode.add_argument('--evaluate',action='store_true')
    args=p.parse_args()
    for name in (args.tag,args.name):
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',name):p.error('Invalid sample/result name')
    args.campaign=args.campaign.resolve()
    if args.evaluate:evaluate(args)
    else:
        if args.baseline is None or args.folds_from is None:p.error('--baseline and --folds-from are required')
        args.baseline=args.baseline.resolve();args.folds_from=args.folds_from.resolve()
        run(args)


if __name__=='__main__':main()
