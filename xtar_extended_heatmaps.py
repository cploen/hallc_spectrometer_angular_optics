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


def rms_percent_change(candidate, reference):
    """Prior comparison convention: negative is smaller RMS, positive larger."""
    if reference > 0:
        return 100*(candidate/reference-1)
    return 0. if candidate == 0 else float('nan')


def write_report(out, rows, counts, *, labels, assumptions):
    """Prior four-column layout and ±20% RMS-change color convention.

    Blue ONLY means smaller RMS and red ONLY larger RMS. Absolute RMS and bias
    use neutral grays. Signed bias is text; gray shading reflects its magnitude.
    Every occupied cell has a count. Existing artifacts are never overwritten.
    """
    import json
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    if set(labels) != set(MODELS) or not assumptions:
        raise ValueError('Model labels and offline reconstruction assumptions required')
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    with (out/'metrics.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    (out/'population_counts.json').write_text(json.dumps(counts,indent=2)+'\n')
    guide = '''START WITH core_RMS_change.png or supported_noncore_RMS_change.png.
Each column is one variable: xptar, yptar, ytar, xtar.
Rows compare completed refit vs oldfit, predictor vs oldfit, predictor vs refit.
BLUE = smaller RMS (better); RED = larger RMS (worse); WHITE = no change.
Colors saturate at +/-20%, matching the previous comparison plot convention.
Every colored cell prints the actual percentage, signed absolute RMS change
in the column's units, and event count. Saturation never caps the printed value.
Absolute-RMS and bias plots use neutral gray, never the improvement colors.
Gray intensity represents absolute RMS or magnitude of bias, respectively.
N is the number of identical reserved events in that foil/delta cell.
* marks N<10. Uncovered cells are blank. Original delta bins never change.
Negative bias describes a sign; it does not by itself mean improvement.
'''
    (out/'PLOT_GUIDE.txt').write_text(guide)
    (out/'assumptions.txt').write_text(assumptions+'\n'+guide)
    short = dict(zip(MODELS,['Preserved oldfit reference','Completed 2018 refit','Refit + separate xtar predictor']))
    pairs = [(MODELS[1],MODELS[0]),(MODELS[2],MODELS[0]),(MODELS[2],MODELS[1])]
    for population in ('core','supported_noncore'):
        local = [r for r in rows if r['population']==population]
        foils = sorted({r['foil_cm'] for r in local})
        lookup = {(r['model'],r['target'],r['foil_cm'],r['delta_low']):r for r in local}
        for kind in ('RMS_change','absolute_RMS','bias'):
            cases = pairs if kind=='RMS_change' else [(m,None) for m in MODELS]
            fig,grid = plt.subplots(3,4,figsize=(48,15),squeeze=False,constrained_layout=True)
            for j,(target,unit) in enumerate(zip(TARGETS,UNITS)):
                if kind=='RMS_change':
                    norm=Normalize(-20,20)
                    cmap=plt.get_cmap('RdBu_r').copy()
                else:
                    field='rms_including_bias' if kind=='absolute_RMS' else 'bias'
                    limit=max([abs(r[field]) for r in local if r['target']==target and r['n']] or [1.])
                    norm=Normalize(0,max(limit,1e-12))
                    cmap=plt.get_cmap('Greys').copy()
                cmap.set_bad('white')
                for i,(model,base) in enumerate(cases):
                    ax=grid[i,j]
                    values=np.full((len(foils),len(EDGES)-1),np.nan)
                    annotations=[]
                    for zi,z in enumerate(foils):
                        for di,lo in enumerate(EDGES[:-1]):
                            r=lookup[model,target,z,float(lo)]
                            if not r['n']:continue
                            count=f"N={r['n']:,}"+('*' if r['sparse'] else '')
                            if kind=='RMS_change':
                                b=lookup[base,target,z,float(lo)]
                                delta=r['rms_including_bias']-b['rms_including_bias']
                                percent=rms_percent_change(r['rms_including_bias'],b['rms_including_bias'])
                                # Zero reference with nonzero candidate is unambiguously worse,
                                # but its percentage is undefined. Mark it rather than divide by zero.
                                v=percent if np.isfinite(percent) else 20.
                                label=(f'{percent:+.1f}%' if np.isfinite(percent) else 'reference RMS=0')+f'\nΔ={delta:+.3g}\n'+count
                            elif kind=='absolute_RMS':
                                v=r['rms_including_bias'];label=f'{v:.3g}\n'+count
                            else:
                                v=abs(r['bias']);label=f"{r['bias']:+.3g}\n"+count
                            values[zi,di]=v
                            annotations.append((di,zi,label,v))
                    im=ax.imshow(np.ma.masked_invalid(values),norm=norm,cmap=cmap,aspect='auto',origin='upper')
                    for di,zi,label,v in annotations:
                        rgba=cmap(norm(v));luminance=np.dot(rgba[:3],[.2126,.7152,.0722])
                        ax.text(di,zi,label,ha='center',va='center',fontsize=10,
                                color='black' if luminance>.5 else 'white')
                    subtitle=short[model]+('\nminus '+short[base] if base else '')
                    ax.set_title(f'{target} ({unit})\n{subtitle}',fontsize=13)
                    ax.set_xticks(range(len(EDGES)-1),[f'[{lo:g},{hi:g})' for lo,hi in zip(EDGES[:-1],EDGES[1:])],rotation=35,ha='right')
                    ax.set_yticks(range(len(foils)),[f'{z:g}' for z in foils])
                    ax.tick_params(labelsize=10)
                    ax.set(xlabel='Original delta (%)',ylabel='Physical foil (cm)')
                    ax.set_xticks(np.arange(-.5,len(EDGES)-1,1),minor=True)
                    ax.set_yticks(np.arange(-.5,len(foils),1),minor=True)
                    ax.grid(which='minor',color='0.85',linewidth=.4)
                    ax.tick_params(which='minor',bottom=False,left=False)
                if kind!='RMS_change':
                    bar_label=('RMS including bias' if kind=='absolute_RMS' else '|bias|')+' ('+unit+')'
                    fig.colorbar(im,ax=grid[:,j].tolist(),label=bar_label,shrink=.7)
            if kind=='RMS_change':
                fig.colorbar(im,ax=grid.ravel().tolist(),label='RMS change (%): BLUE better / RED worse',extend='both',shrink=.75)
                title='RMS CHANGE: BLUE = BETTER (smaller)   |   RED = WORSE (larger)   |   WHITE = unchanged\nColors saturate at ±20%; actual percentages and signed ΔRMS are printed in every occupied cell.'
            else:
                title=('ABSOLUTE RMS INCLUDING BIAS — gray intensity = RMS; numbers in column units' if kind=='absolute_RMS'
                       else 'SIGNED BIAS — numbers give direction; neutral gray intensity = |bias|, not improvement')
            fig.suptitle(population+' — '+title+'\nIdentical reserved events; fixed original delta. * N<10; blank = uncovered. Offline comparison; oldfit historical provenance remains qualified.',fontsize=18)
            path=out/f'{population}_{kind}.png'
            fig.savefig(path,dpi=130)
            plt.close(fig)
            print('Wrote',path,flush=True)


def replot(source, output):
    """Read saved cell metrics only. No event loading, selection or fitting."""
    import hashlib
    import json
    source,output=Path(source),Path(output)
    csv_path=source/'metrics.csv'
    before=hashlib.sha256(csv_path.read_bytes()).hexdigest()
    with csv_path.open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    numeric={'foil_cm','delta_low','delta_high','bias','rms_including_bias',
             'bias_change_vs_oldfit','bias_change_vs_refit','rms_change_vs_oldfit','rms_change_vs_refit'}
    for row in rows:
        row['n']=int(row['n']);row['sparse']=row['sparse'].lower()=='true'
        for key in numeric:
            if row[key]!='':row[key]=float(row[key])
    counts=json.loads((source/'population_counts.json').read_text())
    parent_manifest=source.parent/'manifest.json'
    if parent_manifest.exists():
        labels=json.loads(parent_manifest.read_text())['model_labels']
    else:
        labels={m:m for m in MODELS}
    assumptions=(source/'assumptions.txt').read_text()
    write_report(output,rows,counts,labels=labels,assumptions=assumptions)
    after=hashlib.sha256(csv_path.read_bytes()).hexdigest()
    if before!=after:raise RuntimeError('Original metrics changed during replot')
    (output/'replot_provenance.json').write_text(json.dumps(dict(source=str(source),
        original_metrics_sha256=before,source_metrics_unchanged=True,
        plotting_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        method='plot saved cell metrics only; no fitting or event reselection'),indent=2)+'\n')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description='Replot saved xtar comparison metrics using the prior blue-better/red-worse scheme. No fit.')
    parser.add_argument('--replot-from',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    replot(args.replot_from,args.out)
