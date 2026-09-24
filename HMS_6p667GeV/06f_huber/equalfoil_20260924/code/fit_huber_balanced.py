#!/usr/bin/env python3
"""Four fixed full-basis fits with equal total base weight per physical foil."""
import argparse
import importlib.metadata
import json
import re
import shutil
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from core_sample import PROJECT, digest, write_tsv
from fit_elastic import matrix_rows, basis, problem, design, write_matrix
from fit_huber import load_frozen, sample_ladder, summarize, paired_intervals, coefficients
from elastic_net import TARGETS, SCALES
from huber_optics import FixedBasis

REFERENCES=('balanced_squared_full','all_core_squared_full','all_core_huber_full',
            'all_supported_squared_full','all_supported_huber_full')
NEW_MODELS=tuple(f'{pool}_{loss}_equalfoil' for pool in ('all_core','all_supported') for loss in ('squared','huber'))


def foil_weights(foils):
    foils=np.asarray(foils)
    if foils.ndim!=1 or not len(foils) or not np.isfinite(foils).all():
        raise ValueError('Finite physical foil labels are required')
    values,inverse,counts=np.unique(foils,return_inverse=True,return_counts=True)
    return len(foils)/(len(values)*counts[inverse])


def weight_audit(a,base,residual,sigma,delta,model,target,robust):
    """IRLS weights are a contribution diagnostic, not full coefficient influence."""
    factor=1/np.hypot(1.,residual/(delta*sigma)) if robust else np.ones(len(base))
    combined=base*factor
    records=[]
    masks=[('all','all',np.ones(len(base),bool))]
    for z in np.unique(a['zfoil']): masks.append(('foil',f'{z:g}',a['zfoil']==z))
    for q,name in ((2,'core'),(1,'shoulder')): masks.append(('quality',name,a['quality']==q))
    for z in np.unique(a['zfoil']):
        for q,name in ((2,'core'),(1,'shoulder')):
            masks.append(('foil_quality',f'{z:g}:{name}',(a['zfoil']==z)&(a['quality']==q)))
    for kind,group,h in masks:
        if not h.any(): continue
        b,e=base[h],combined[h]
        records.append(dict(model=model,target=target,group_kind=kind,group=group,n=int(sum(h)),
            raw_share=float(h.mean()),base_sum=float(b.sum()),base_share=float(b.sum()/base.sum()),
            irls_sum=float(e.sum()),irls_share=float(e.sum()/combined.sum()),
            mean_robust_factor=float(e.sum()/b.sum()),factor_below_half=float(np.mean(factor[h]<.5)),
            effective_n=float(e.sum()**2/(e@e))))
    return records


def run(campaign,name,source,threads):
    parent=campaign/'06f_huber'/source
    previous=json.loads((parent/'manifest.json').read_text())
    if previous.get('status')!='complete': raise ValueError('Requires completed unweighted study')
    a,inputs,archived=load_frozen(campaign)
    for path,sha in previous['inputs'].items():
        if inputs.get(path)!=sha: raise ValueError(f'Source input identity changed: {path}')
    inputs[str(parent/'manifest.json')]=digest(parent/'manifest.json')
    for rel in ('residuals.npz','seed.dat',*[f'matrices/{m}.dat' for m in REFERENCES]):
        path=parent/rel
        if digest(path)!=previous['outputs'][rel]: raise ValueError(f'Source output changed: {rel}')
        inputs[str(path)]=digest(path)
    ladder=sample_ladder(a)
    evaluation=(a['sample']==2)&~a['excluded']&(a['quality']>=1)
    train={k:ladder[k] for k in ('all_core','all_supported')}
    for k,h in train.items():
        if np.any(h&evaluation) or h.sum()!=previous['counts'][k]: raise ValueError('Membership changed')
    ae={k:v[evaluation] for k,v in a.items()}
    with np.load(parent/'residuals.npz') as saved:
        for k in ('entry','rungroup','quality','ndel','zfoil'):
            if not np.array_equal(ae[k],saved[k]): raise ValueError('Reserved membership mismatch')
        reference_residuals={m:saved[m].copy() for m in REFERENCES}
    sigma=np.array(previous['residual_scales_physical'])/SCALES
    delta=previous['delta']
    rows=matrix_rows(parent/'seed.dat');terms=basis(rows)
    out=campaign/'06f_huber'/name
    if out.exists(): raise FileExistsError(f'Immutable output exists: {out}')
    out.mkdir(parents=True)
    for sub in ('tsv','plots','matrices','code'): (out/sub).mkdir()
    for file in ('fit_huber_balanced.py','balanced_huber_results.py','huber_optics.py','fit_huber.py',
                 'fit_elastic.py','elastic_net.py','core_sample.py','reallocate_core.py','core_balance.py',
                 'spectrometer_config.py','spectrometer_config.h','spectrometer_profiles.def','preallocated_svd.py'):
        shutil.copy2(PROJECT/file,out/'code'/file)
    (out/'code/tests').mkdir()
    for file in ('test_huber.py','test_huber_weights.py'):
        shutil.copy2(PROJECT/'tests'/file,out/'code/tests'/file)
    shutil.copy2(parent/'seed.dat',out/'seed.dat')
    report=dict(schema='huber_equalfoil_v1',status='running',source=source,inputs=inputs,
        counts={k:int(h.sum()) for k,h in train.items()},delta=delta,residual_scales_physical=previous['residual_scales_physical'],
        weight_policy='w_i = N/(K*N_foil); equal total base weight per physical foil, normalized mean one',
        residual_policy='Base weights multiply loss outside pseudo-Huber; threshold fixed in physical residual units',
        sample_policy='Exact previous all-core/all-supported training IDs; reserved events unchanged',
        basis_terms=len(terms),threads=threads,
        versions={m:importlib.metadata.version(m) for m in ('numpy','scipy','matplotlib','uproot')})
    (out/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    x,y,offset=problem(a,rows);xe,ye=x[evaluation],y[evaluation]
    for m in REFERENCES:
        coef=coefficients(parent/f'matrices/{m}.dat',terms)
        if not np.allclose((xe@coef-ye)*SCALES,reference_residuals[m],atol=1e-7,rtol=1e-8):
            raise ValueError('Reference residuals cannot be reproduced')
        shutil.copy2(parent/f'matrices/{m}.dat',out/f'matrices/{m}.dat')
    fits={};audits=[];weights=[];weight_archive={};coverage=[]
    # Fit all four prespecified candidates before computing any new reserved scores.
    with threadpool_limits(limits=threads):
        for pool,sel in train.items():
            at={k:v[sel] for k,v in a.items()}
            w=foil_weights(at['zfoil'])
            for z in np.unique(at['zfoil']):
                if not np.isclose(w[at['zfoil']==z].sum()/w.sum(),.2):
                    raise ValueError('Expected equal weights for the five physical foils')
                for d in np.unique(at['ndel']):
                    h=(at['zfoil']==z)&(at['ndel']==d)
                    coverage.append(dict(pool=pool,zfoil=z,ndel=int(d),n=int(sum(h)),
                        core=int(sum(h&(at['quality']==2))),shoulder=int(sum(h&(at['quality']==1))),
                        base_sum=float(w[h].sum()),base_share=float(w[h].sum()/w.sum())))
            weight_archive[pool+'_rungroup']=at['rungroup'];weight_archive[pool+'_entry']=at['entry']
            weight_archive[pool+'_base_weight']=w
            print(f'{pool}: {sum(sel):,} events; each physical foil has 20% base weight',flush=True)
            xt,yt=x[sel],y[sel]
            factor=FixedBasis(xt,weights=w)
            for loss in ('squared','huber'):
                label=f'{pool}_{loss}_equalfoil';coef=[]
                for j,target in enumerate(TARGETS):
                    b,info=factor.fit(yt[:,j],sigma[j] if loss=='huber' else None,delta=delta)
                    coef.append(b);audits.append(dict(model=label,target=target,**info))
                    weights+=weight_audit(at,w,xt@b-yt[:,j],sigma[j],delta,label,target,loss=='huber')
                    print(f'  {loss}/{target}: rank={info["rank"]}, iterations={info["iterations"]}',flush=True)
                fits[label]=np.column_stack(coef)
            del factor,xt,yt
    residuals=dict(reference_residuals)
    for label,coef in fits.items():
        residuals[label]=(xe@coef-ye)*SCALES
        write_matrix(out/f'matrices/{label}.dat',rows,terms,coef)
        restored=matrix_rows(out/f'matrices/{label}.dat')
        pred=design(ae,[t for c,t in restored])@np.array([c[:3] for c,t in restored])
        if not np.allclose(pred,xe@coef+offset[evaluation],rtol=1e-10,atol=1e-10):
            raise ValueError('Exported matrix does not reproduce predictions')
    write_tsv(out/'tsv/convergence.tsv',audits,fields=sorted(set().union(*(r.keys() for r in audits))))
    write_tsv(out/'tsv/weight_contributions.tsv',weights)
    write_tsv(out/'tsv/coverage.tsv',coverage)
    np.savez_compressed(out/'training_weights.npz',**weight_archive)
    print('All four fits complete; computing common reserved-sample diagnostics',flush=True)
    stats,labels=summarize(ae,residuals,out)
    intervals=paired_intervals(ae,residuals,labels,out)
    np.savez_compressed(out/'residuals.npz',**residuals,quality=ae['quality'],zfoil=ae['zfoil'],
                        ndel=ae['ndel'],entry=ae['entry'],rungroup=ae['rungroup'])
    from balanced_huber_results import products
    products(out,stats,intervals,weights,report)
    report['status']='complete'
    report['reserved_counts']={p:int(sum(ae['quality']==q)) for p,q in [('core',2),('shoulder',1)]}
    report['outputs']={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and p.name!='manifest.json'}
    (out/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'Complete: {out}',flush=True)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path);p.add_argument('name',nargs='?',default='equalfoil_20260924')
    p.add_argument('--source',default='expansion_20260924');p.add_argument('--threads',type=int,default=4)
    args=p.parse_args()
    if args.threads<1 or any(not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*',v) for v in (args.name,args.source)):
        p.error('Invalid name or threads')
    run(args.campaign.resolve(),args.name,args.source,args.threads)
