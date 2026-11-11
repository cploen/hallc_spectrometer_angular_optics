"""Four-variable diagnostics on fixed reserved IDs and original delta coordinates.

This module consumes predictions from an audited reconstruction adapter. It
neither reconstructs events nor changes event selection. All prediction
bundles must retain exactly the same ordered (rungroup, entry) IDs.
"""
from pathlib import Path
import csv
import numpy as np

EDGES = np.array([-15., -12., -10., -8., -5., 0., 5., 8., 10., 12., 14., 15.])
TARGETS = ('xptar', 'yptar', 'ytar', 'xtar')
# Input slopes rad and positions cm; reporting angles mrad, positions mm.
SCALES = np.array([1000., 1000., 10., 10.])
UNITS = ('mrad', 'mrad', 'mm', 'mm')
MODELS = ('oldfit_reference', 'completed_refit', 'refit_plus_xtar_predictor')


def ordered_ids(a):
    ids = list(zip(map(str, a['rungroup']), map(int, a['entry'])))
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate event IDs in diagnostic bundle')
    return ids


def metrics(a, truth, predictions, *, physical_foils, supported, qa_min=10):
    """Return rectangular foil/delta metrics and explicit population counts.

    a supplies frozen delta (%), physical ztarT (cm), quality and excluded.
    supported must be determined from TRAINING coverage, not candidate values.
    predictions maps three model names to {rungroup, entry, values}; each values
    array is (n,4), ordered as TARGETS. truth has the same order and native units.
    Empty cells have n=0 and blank values; no nonfinite or missing predictions
    are silently dropped. Delta bins are half-open, including [14,15).
    """
    ids = ordered_ids(a)
    n = len(ids)
    truth = np.asarray(truth, dtype=float)
    supported = np.asarray(supported)
    if truth.shape != (n, 4) or not np.isfinite(truth).all():
        raise ValueError('Expected finite four-variable diagnostic truth')
    if supported.shape != (n,) or supported.dtype.kind != 'b':
        raise ValueError('Expected explicit training-support mask')
    if set(predictions) != set(MODELS):
        raise ValueError('Expected the three fixed comparison models')
    residuals = {}
    for model, bundle in predictions.items():
        if ordered_ids(bundle) != ids:
            raise ValueError('Candidate changed reserved IDs/order')
        p = np.asarray(bundle['values'], dtype=float)
        if p.shape != truth.shape or not np.isfinite(p).all():
            raise ValueError('Missing/nonfinite model predictions; do not silently drop events')
        residuals[model] = (p-truth)*SCALES
    delta = np.asarray(a['delta'], dtype=float)
    foil = np.asarray(a['ztarT'], dtype=float)
    quality = np.asarray(a['quality'])
    blocked = np.asarray(a['excluded'])
    if (any(v.shape != (n,) for v in (delta,foil,quality,blocked)) or
            blocked.dtype.kind != 'b' or not np.isfinite(delta).all() or
            not np.isfinite(foil).all() or not np.isin(quality, [0,1,2]).all()):
        raise ValueError('Invalid original coordinate/population metadata')
    foils = np.asarray(sorted(physical_foils), dtype=float)
    if not len(foils) or len(np.unique(foils)) != len(foils) or not np.isin(foil, foils).all():
        raise ValueError('Physical-foil list must cover all accepted events')
    bins = np.searchsorted(EDGES, delta, side='right')-1
    in_range = (bins >= 0)&(bins < len(EDGES)-1)
    populations = dict(core=(~blocked)&(quality==2),
                       supported_noncore=(~blocked)&(quality!=2)&supported,
                       unsupported_noncore=(~blocked)&(quality!=2)&(~supported),
                       blocked=blocked)
    counts = dict(total=n, in_extended_delta=int(in_range.sum()),
                  outside_extended_delta=int((~in_range).sum()),
                  unsupported_core=int(((quality==2)&(~supported)&(~blocked)).sum()))
    for name, mask in populations.items():
        counts[name] = int(mask.sum())
        counts[name+'_outside_extended_delta'] = int((mask&(~in_range)).sum())
    rows = []
    for population in ('core','supported_noncore'):
        for z in foils:
            for b, (lo,hi) in enumerate(zip(EDGES[:-1], EDGES[1:])):
                mask = populations[population]&(foil==z)&(bins==b)&in_range
                count = int(mask.sum())
                values = {}
                for model, residual in residuals.items():
                    for j, target in enumerate(TARGETS):
                        v = residual[mask,j]
                        values[model,target] = (float(v.mean()), float(np.sqrt(np.mean(v*v)))) if count else ('','')
                for j,target in enumerate(TARGETS):
                    for model in MODELS:
                        bias,rms = values[model,target]
                        row = dict(population=population,foil_cm=float(z),delta_low=float(lo),delta_high=float(hi),
                                   target=target,unit=UNITS[j],model=model,n=count,sparse=bool(0<count<qa_min),
                                   bias=bias,rms_including_bias=rms)
                        for reference, suffix in [(MODELS[0],'oldfit'),(MODELS[1],'refit')]:
                            rb,rr = values[reference,target]
                            row['bias_change_vs_'+suffix] = bias-rb if count else ''
                            row['rms_change_vs_'+suffix] = rms-rr if count else ''
                        rows.append(row)
    return rows, counts


def write_report(out, rows, counts, *, labels, assumptions):
    """Write readable wide plots and CSV. Caller supplies factual model labels.

    Negative RMS changes mean improvement. Shared absolute plot ranges prevent
    a candidate from appearing improved merely by independent color rescaling.
    The label for oldfit must retain a provenance qualification until resolved.
    """
    import json
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if set(labels) != set(MODELS) or not assumptions:
        raise ValueError('Model labels and offline reconstruction assumptions required')
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    with (out/'metrics.csv').open('w', newline='') as stream:
        writer=csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    (out/'population_counts.json').write_text(json.dumps(counts,indent=2)+'\n')
    (out/'assumptions.txt').write_text(assumptions+'\nNegative RMS changes mean improvement. Empty cells are uncovered. Asterisk denotes a sparse cell.\n')
    for population in ('core','supported_noncore'):
        for target,unit in zip(TARGETS,UNITS):
            selected=[r for r in rows if r['population']==population and r['target']==target]
            foils=sorted({r['foil_cm'] for r in selected})
            lookup={(r['model'],r['foil_cm'],r['delta_low']):r for r in selected}
            panels=[]
            for statistic in ('bias','rms_including_bias'):
                for model in MODELS:
                    panels.append((model,statistic,labels[model]+' — '+statistic.replace('_',' ')))
            for reference,suffix in [(MODELS[0],'oldfit'),(MODELS[1],'refit')]:
                for stat in ('bias','rms'):
                    panels.append((MODELS[2],stat+'_change_vs_'+suffix,
                        labels[MODELS[2]]+' minus '+labels[reference]+' — '+stat+' change'))
            arrays=[]
            for model,stat,title in panels:
                arrays.append(np.array([[lookup[model,z,float(lo)][stat]
                    if lookup[model,z,float(lo)][stat]!='' else np.nan for lo in EDGES[:-1]] for z in foils]))
            category=lambda s:next(p for p in ('rms_including_bias','bias_change','rms_change','bias') if s.startswith(p))
            limits={}
            for prefix in ('bias','rms_including_bias','bias_change','rms_change'):
                values=np.concatenate([a[np.isfinite(a)] for a,(_,s,_) in zip(arrays,panels) if category(s)==prefix])
                limits[prefix]=max(float(np.max(np.abs(values))) if len(values) else 0.,1e-12)
            fig,axes=plt.subplots(len(panels),1,figsize=(30,4.6*len(panels)),constrained_layout=True)
            for ax,a,(model,stat,title) in zip(axes,arrays,panels):
                prefix=category(stat)
                limit=limits[prefix]
                positive=stat=='rms_including_bias'
                cmap=plt.get_cmap('viridis' if positive else 'RdBu_r').copy();cmap.set_bad('white')
                im=ax.imshow(np.ma.masked_invalid(a),origin='lower',aspect='auto',cmap=cmap,
                             vmin=0 if positive else -limit,vmax=limit)
                for i,z in enumerate(foils):
                    for j,lo in enumerate(EDGES[:-1]):
                        r=lookup[model,z,float(lo)]
                        if r['n']:
                            ax.text(j,i,f"{a[i,j]:+.3g}"+ ('*' if r['sparse'] else '')+f"\nn={r['n']}",
                                    ha='center',va='center',fontsize=11,
                                    bbox=dict(facecolor='white',alpha=.82,edgecolor='none'))
                ax.set_xticks(range(len(EDGES)-1),[f'[{lo:g}, {hi:g})' for lo,hi in zip(EDGES[:-1],EDGES[1:])])
                ax.set_yticks(range(len(foils)),[f'{z:g}' for z in foils])
                ax.set(title=title+' ('+unit+')',ylabel='Physical foil (cm)',xlabel='Original delta (%)')
                fig.colorbar(im,ax=ax,label=unit)
            fig.suptitle(target+' — '+population+' — identical reserved IDs; fixed original delta\n'
                         'Offline proof of principle. Negative RMS change = improvement. * Sparse cell. Blank = uncovered.',fontsize=18)
            fig.savefig(out/f'{population}_{target}.png',dpi=130)
            plt.close(fig)
