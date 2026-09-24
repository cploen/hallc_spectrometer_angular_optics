#!/usr/bin/env python3
"""Frozen-sample expansion study: equal15 -> all cores -> supported shoulders.

Reuses the historical elastic-net-selected and beam-search term sets. Does not
repeat their searches, alter density labels, or install a replay matrix.
"""
import argparse
import importlib.metadata
import json
import re
import shutil
import time
from pathlib import Path
import numpy as np
import uproot
from threadpoolctl import threadpool_limits
from core_sample import PROJECT, digest, ranks, read_tsv, write_tsv
from reallocate_core import load_inputs, checked
from fit_elastic import matrix_rows, basis, problem, design, write_matrix
from elastic_net import TARGETS, SCALES
from huber_optics import FixedBasis, robust_scale

LADDER = ('balanced', 'half_surplus', 'all_core', 'half_shoulders', 'all_supported')


def load_frozen(campaign):
    inp = load_inputs(campaign, 'min10')
    hashes = {}
    def record(path, expected=None):
        if expected is not None: checked(path, expected)
        hashes[str(path.resolve())] = digest(path)
    for path in (inp['source']/'manifest.json', inp['table'], inp['metadata'], inp['maskpath']):
        record(path)
    pieces = []
    for row in inp['settings']:
        name = row['rungroup']
        path = inp['source']/f'root/CoreSample_{name}.root'
        record(path, inp['manifest']['outputs'][str(path.relative_to(inp['source']))])
        with uproot.open(path) as root: a = root['CoreSample'].arrays(library='np')
        n = len(a['entry'])
        if a['entry'].dtype.kind not in 'iu' or len(np.unique(a['entry'])) != n:
            raise ValueError('Duplicate or noninteger frozen IDs')
        if not np.isin(a['sample'], [0,1,2,3]).all() or not np.isin(a['quality'], [0,1,2]).all():
            raise ValueError('Invalid frozen sample or quality')
        if np.any(a['run'] != int(row['optics_id'])): raise ValueError('Wrong run identity')
        meta = inp['meta'][name]
        spec = __import__('spectrometer_config').from_campaign(campaign)
        if np.any(a['foil'] < 0) or np.any(a['foil'] >= len(meta['foils'])):
            raise ValueError('Invalid foil index')
        if np.any(a['ndel'] < 0) or np.any(a['ndel'] >= len(meta['edges'])-1):
            raise ValueError('Invalid delta index')
        if np.any(a['xscol'] < 0) or np.any(a['xscol'] >= spec.nx) or np.any(a['yscol'] < 0) or np.any(a['yscol'] >= spec.ny):
            raise ValueError('Invalid hole index')
        lo, hi = np.array(meta['edges'])[a['ndel']], np.array(meta['edges'])[a['ndel']+1]
        if not np.array_equal(a['zfoil'], np.array(meta['foils'])[a['foil']]):
            raise ValueError('Physical foil mismatch')
        if not (np.array_equal(lo,a['delta_low']) and np.array_equal(hi,a['delta_high']) and
                np.all((a['delta'] >= lo)&(a['delta'] < hi))):
            raise ValueError('Delta metadata mismatch')
        a['rungroup'] = np.full(n, name)
        a['angle'] = np.full(n, np.deg2rad(meta['angle_deg']))
        a['ztarT'] = a['zfoil'].copy()
        a['xsT'] = np.array([spec.xs(int(v)) for v in a['xscol']])
        a['ysT'] = np.array([spec.ys(int(v)) for v in a['yscol']])
        # Same HMS targetTruth equations as spectrometer_config.h and TFit export.
        c, s = np.cos(a['angle']), np.sin(a['angle'])
        deg = min(abs(meta['angle_deg']), 40.)
        ymis = .1*(.52-.012*deg+.002*deg**2)
        beam = -a['xbpm_tar']
        a['xptarT'] = a['xsT']/(spec.sieve_distance-a['zfoil']*c)
        a['yptarT'] = (a['ysT']-a['zfoil']*s-beam*c+ymis)/(spec.sieve_distance-a['zfoil']*c)
        a['ytarT'] = a['zfoil']*(s-a['yptarT']*c)+beam*(c+a['yptarT']*s)-ymis
        a['excluded'] = np.array([(int(x),int(y)) in inp['blocked'] for x,y in zip(a['xscol'],a['yscol'])])
        fields = ('entry','rungroup','sample','quality','excluded','zfoil','foil','ndel','xscol','yscol',
                  'delta','angle','xfp','xpfp','yfp','ypfp','xtar','xptarT','yptarT','ytarT','ztarT','xsT','ysT')
        for k in fields:
            if k != 'rungroup' and not np.isfinite(a[k]).all(): raise ValueError(f'Nonfinite {name}/{k}')
        pieces.append({k:a[k] for k in fields})
    a = {k:np.concatenate([p[k] for p in pieces]) for k in pieces[0]}
    source = campaign/'06e_beam_search/beam10'
    manifest = json.loads((source/'manifest.json').read_text())
    record(source/'manifest.json')
    for rel in ('seed.dat','matrices/svd.dat','matrices/start.dat','matrices/beam.dat','tsv/folds.tsv'):
        record(source/rel, manifest['outputs'][rel])
    # Independently join the saved equal15 ID table and archived beam training IDs.
    idpath = campaign/'05c_core_sample/equal15/tsv/selected_ids.tsv'
    record(idpath)
    ids = [(r['rungroup'],int(r['entry'])) for r in read_tsv(idpath)]
    beam_ids = [(r['rungroup'],int(r['entry'])) for r in read_tsv(source/'tsv/folds.tsv')]
    if len(set(ids)) != len(ids) or set(ids) != set(beam_ids) or len(ids) != len(beam_ids):
        raise ValueError('equal15 and frozen beam training IDs differ')
    wanted = set(ids)
    a['balanced'] = np.array([(r,int(e)) in wanted for r,e in zip(a['rungroup'],a['entry'])])
    if sum(a['balanced']) != len(ids): raise ValueError('Training IDs missing from frozen trees')
    valid = (a['sample'] != 2)&(a['quality']==2)&~a['excluded']
    if np.any(a['balanced']&~valid): raise ValueError('Balanced IDs include held-out or invalid events')
    _, counts = np.unique(a['zfoil'][a['balanced']], return_counts=True)
    if len(set(counts)) != 1: raise ValueError('Baseline no longer foil balanced')
    return a, hashes, source


def sample_ladder(a):
    available = (a['sample'] != 2)&~a['excluded']
    core = available&(a['quality']==2)
    shoulder = available&(a['quality']==1)
    half = np.zeros(len(core), dtype=bool)
    for rg in np.unique(a['rungroup']):
        h = a['rungroup']==rg
        half[h] = ranks(667, rg, a['entry'][h], 'huber_expansion') < np.uint64(2**63)
    b = a['balanced']
    return dict(balanced=b, half_surplus=b|(core&half), all_core=core,
                half_shoulders=core|(shoulder&half), all_supported=core|shoulder)


def coefficients(path, terms):
    rows = matrix_rows(path)
    values = {t:c for c,t in rows}
    if set(terms)-values.keys(): raise ValueError('Matrix missing fitted basis rows')
    return np.array([values[t][:3] for t in terms])


def summarize(a, residuals, out):
    rows, groups = [], []
    pools = dict(core=a['quality']==2, shoulder=a['quality']==1, supported=a['quality']>=1)
    for pool,h in pools.items():
        groups.append((pool,'overall','all',h))
        for foil in np.unique(a['zfoil'][h]):
            groups.append((pool,'foil',f'{foil:g}',h&(a['zfoil']==foil)))
            for d in np.unique(a['ndel'][h]):
                groups.append((pool,'foil_delta',f'{foil:g}:{d}',h&(a['zfoil']==foil)&(a['ndel']==d)))
    for model,res in residuals.items():
        for pool,kind,group,h in groups:
            if not h.any(): continue
            e = res[h]
            for j,target in enumerate(TARGETS):
                v = e[:,j]
                rows.append(dict(model=model,pool=pool,group_kind=kind,group=group,target=target,n=len(v),
                    bias=float(v.mean()),spread=float(v.std()),rms=float(np.sqrt(np.mean(v*v))),
                    width68=float((np.quantile(v,.84)-np.quantile(v,.16))/2),
                    p95=float(np.quantile(np.abs(v),.95))))
    write_tsv(out/'tsv/metrics.tsv',rows)
    # Record sparse holes rather than letting high-rate holes hide their behavior.
    labels = np.array([f'{r}:{z:g}:{d}:{x}:{y}' for r,z,d,x,y in zip(a['rungroup'],a['zfoil'],a['ndel'],a['xscol'],a['yscol'])])
    hole_rows = []
    for pool,h in pools.items():
        names,inv = np.unique(labels[h],return_inverse=True)
        n = np.bincount(inv)
        for model,res in residuals.items():
            e = res[h]
            for j,target in enumerate(TARGETS):
                mean = np.bincount(inv,weights=e[:,j])/n
                mse = np.bincount(inv,weights=e[:,j]**2)/n
                for i,name in enumerate(names):
                    hole_rows.append(dict(model=model,pool=pool,hole=name,target=target,n=int(n[i]),
                                          qa_low=bool(n[i]<10),bias=float(mean[i]),rms=float(np.sqrt(mse[i]))))
    write_tsv(out/'tsv/holes.tsv',hole_rows)
    return rows, labels


def paired_intervals(a, residuals, labels, out):
    """Conditional uncertainty: paired, stratified run/foil/delta/hole bootstrap.

    Resamples holes within physical foil/delta, never refits models. Does not
    account for training uncertainty, label error or previous holdout inspection.
    """
    rng = np.random.default_rng(667)
    baseline = residuals['balanced_squared_full']
    records = []
    for pool in ('core','shoulder','supported'):
        h = (a['quality']==(2 if pool=='core' else 1)) if pool!='supported' else a['quality']>=1
        names,inv = np.unique(labels[h],return_inverse=True)
        count = np.bincount(inv)
        cell = np.array([':'.join(s.split(':')[1:3]) for s in names])
        strata = [np.flatnonzero(cell==c) for c in np.unique(cell)]
        draw = np.concatenate([rng.choice(s,(500,len(s)),replace=True) for s in strata],axis=1)
        denom = count[draw].sum(axis=1)
        base_sums = np.array([np.bincount(inv,weights=baseline[h,j]**2) for j in range(3)]).T
        base_rms = np.sqrt(base_sums[draw].sum(axis=1)/denom[:,None])
        for model,res in residuals.items():
            sums = np.array([np.bincount(inv,weights=res[h,j]**2) for j in range(3)]).T
            boot = 100*(1-np.sqrt(sums[draw].sum(axis=1)/denom[:,None])/base_rms)
            point = 100*(1-np.sqrt(np.mean(res[h]**2,axis=0)/np.mean(baseline[h]**2,axis=0)))
            for j,target in enumerate(TARGETS):
                lo,hi = np.quantile(boot[:,j],[.025,.975])
                records.append(dict(model=model,pool=pool,target=target,rms_gain_percent=float(point[j]),
                                    low95=float(lo),high95=float(hi),clusters=len(names),replicates=500))
    write_tsv(out/'tsv/paired_intervals.tsv',records)
    return records


def run(campaign, name, threads=4, check=False):
    a, hashes, saved = load_frozen(campaign)
    ladder = sample_ladder(a)
    evaluation = (a['sample']==2)&~a['excluded']&(a['quality']>=1)
    counts = {k:int(h.sum()) for k,h in ladder.items()}
    print('Training sample sizes:',counts,flush=True)
    print('Reserved evaluation:',int(sum(evaluation)),flush=True)
    if any(np.any(h&evaluation) for h in ladder.values()): raise ValueError('Evaluation leakage')
    if check: return counts
    out = campaign/'06f_huber'/name
    if out.exists(): raise FileExistsError(f'Immutable output exists: {out}')
    out.mkdir(parents=True)
    for sub in ('tsv','plots','matrices','code'): (out/sub).mkdir()
    rows = matrix_rows(saved/'seed.dat'); terms = basis(rows)
    shutil.copy2(saved/'seed.dat',out/'seed.dat')
    for f in ('fit_huber.py','huber_optics.py','huber_plots.py','fit_elastic.py','elastic_net.py',
              'core_sample.py','reallocate_core.py','core_balance.py','spectrometer_config.py',
              'spectrometer_config.h','spectrometer_profiles.def','preallocated_svd.py'):
        shutil.copy2(PROJECT/f,out/'code'/f)
    manifest = dict(schema='huber_expansion_v1',status='running',counts=counts,inputs=hashes,
        loss='pseudo_huber',delta=1.5,threads=threads,seed=667,versions={m:importlib.metadata.version(m)
        for m in ('numpy','scipy','uproot','matplotlib','scikit-learn')},
        policy='Fixed threshold 1.5 balanced-training residual MAD scales; no validation tuning; no population weights',
        basis_policy='Full polynomial for main ladder; reuse frozen EN-selected and beam10 supports for controls')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    x,y,offset = problem(a,rows)
    ae = {k:v[evaluation] for k,v in a.items()}
    xe,ye = x[evaluation],y[evaluation]
    residuals, fits, audit = {}, {}, []
    coverage = []
    for pool,sel in ladder.items():
        for z in np.unique(a['zfoil']):
            for d in np.unique(a['ndel']):
                h = sel&(a['zfoil']==z)&(a['ndel']==d)
                holes = set(zip(a['rungroup'][h],a['xscol'][h],a['yscol'][h]))
                coverage.append(dict(sample=pool,zfoil=z,ndel=int(d),n=int(sum(h)),holes=len(holes)))
    write_tsv(out/'tsv/coverage.tsv',coverage)
    # Save every assigned ID so expansion and reserved membership are reviewable.
    write_tsv(out/'tsv/membership.tsv', [dict(rungroup=r,entry=int(e),quality=int(q),reserved=bool(v),
        **{k:bool(h[i]) for k,h in ladder.items()}) for i,(r,e,q,v) in enumerate(zip(a['rungroup'],a['entry'],a['quality'],evaluation))])
    sigma = None
    with threadpool_limits(limits=threads):
        for pool,selected in ladder.items():
            clock=time.monotonic()
            print(f'{pool}: factor full design ({sum(selected):,} x {x.shape[1]})',flush=True)
            factor = FixedBasis(x[selected])
            squared = np.column_stack([factor.fit(y[selected,j])[0] for j in range(3)])
            if pool=='balanced':
                sigma=np.array([robust_scale(x[selected]@squared[:,j]-y[selected,j]) for j in range(3)])
                # Adapter equivalence is checked before considering expanded samples.
                old=coefficients(saved/'matrices/svd.dat',terms)
                difference=np.max(np.abs((xe@squared-xe@old)*SCALES),axis=0)
                if np.max(difference)>1e-5: raise ValueError(f'Balanced fit does not reproduce archived SVD: {difference}')
                manifest['balanced_prediction_max_abs_difference']=difference.tolist()
                manifest['residual_scales_physical']=(sigma*SCALES).tolist()
            for loss,coef in [('squared',squared),('huber',None)]:
                label=f'{pool}_{loss}_full'
                if loss=='huber':
                    col=[]
                    for j,target in enumerate(TARGETS):
                        b,info=factor.fit(y[selected,j],sigma[j])
                        col.append(b);audit.append(dict(model=label,target=target,**info))
                        print(f'  {target}: {info["iterations"]} iterations, gradient {info["gradient"]:.2g}',flush=True)
                    coef=np.column_stack(col)
                fits[label]=coef
                residuals[label]=(xe@coef-ye)*SCALES
            del factor
            print(f'  Completed in {time.monotonic()-clock:.1f} s',flush=True)
        # Compact bases are fixed by previous development work; do not choose anew
        # on these reserved residuals. All supports include the fitted constant.
        for kind,file in [('enet_basis','start.dat'),('beam_basis','beam.dat')]:
            archived=coefficients(saved/'matrices'/file,terms)
            residuals[f'balanced_squared_{kind}']=(xe@archived-ye)*SCALES
            fits[f'balanced_squared_{kind}']=archived
            for pool in ('all_core','all_supported'):
                selected=ladder[pool]
                robust=np.zeros_like(archived); square=np.zeros_like(archived)
                for j,target in enumerate(TARGETS):
                    support=np.flatnonzero((archived[:,j]!=0)|(np.arange(len(terms))==0))
                    print(f'{pool}/{kind}/{target}: {len(support)} terms',flush=True)
                    factor=FixedBasis(x[selected][:,support])
                    square[support,j]=factor.fit(y[selected,j])[0]
                    b,info=factor.fit(y[selected,j],sigma[j])
                    robust[support,j]=b
                    audit.append(dict(model=f'{pool}_huber_{kind}',target=target,**info))
                    del factor
                for loss,coef in [('squared',square),('huber',robust)]:
                    label=f'{pool}_{loss}_{kind}'
                    fits[label]=coef;residuals[label]=(xe@coef-ye)*SCALES
    for label,coef in fits.items():
        write_matrix(out/f'matrices/{label}.dat',rows,terms,coef)
        restored=matrix_rows(out/f'matrices/{label}.dat')
        pred=design(ae,[t for c,t in restored])@np.array([c[:3] for c,t in restored])
        if not np.allclose(pred,xe@coef+offset[evaluation],rtol=1e-10,atol=1e-10):
            raise ValueError('Exported matrix does not reproduce predictions')
    write_tsv(out/'tsv/convergence.tsv',audit)
    stats,labels=summarize(ae,residuals,out)
    intervals=paired_intervals(ae,residuals,labels,out)
    np.savez_compressed(out/'residuals.npz',**residuals,quality=ae['quality'],zfoil=ae['zfoil'],
                        ndel=ae['ndel'],entry=ae['entry'],rungroup=ae['rungroup'])
    from huber_plots import plots,report
    plots(out,stats,counts)
    (out/'RESULTS.md').write_text(report(stats,intervals,counts,manifest))
    manifest['status']='complete'
    manifest['reserved_counts']={p:int(sum(ae['quality']==q)) for p,q in [('core',2),('shoulder',1)]}
    manifest['excluded_counts']=dict(blocked=int(sum(a['excluded'])),unsupported=int(sum((a['quality']==0)&~a['excluded'])))
    manifest['outputs']={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and p.name!='manifest.json'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Complete: {out}',flush=True)
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign',type=Path)
    parser.add_argument('name',nargs='?',default='expansion')
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*',args.name) or args.threads<1: parser.error('Invalid name/threads')
    run(args.campaign.resolve(),args.name,args.threads,args.check)
