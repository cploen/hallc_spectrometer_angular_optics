#!/usr/bin/env python3
"""Search angular bases from pooled elastic-net seeds; protected events stay closed."""
import argparse
import csv
import importlib.metadata
import json
import re
import shutil
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from core_sample import PROJECT,digest,read_tsv,write_tsv
from fit_elastic import load_sample,checked,matrix_rows,basis,problem,write_matrix,design
from elastic_convergence import saved_folds
from elastic_net import TARGETS
from beam_search import DEFAULTS,validate,Folds,search,fit_all

SEEDS = dict(xptar=dict(alpha=.00003,l1=.5),ytar=dict(alpha=.00003,l1=.5),
             yptar=dict(alpha=.0001,l1=.5))


def inputs(campaign,tag,source,selectors):
    parent = campaign/'06d_elastic_net'/source
    manifest = json.loads((parent/'manifest.json').read_text())
    if (manifest.get('schema'),manifest.get('status'),manifest.get('sample')) != ('elastic_pool_v1','complete',tag):
        raise ValueError('Beam search requires a completed pooled study for the same sample tag')
    hashes = {str(parent/'manifest.json'):digest(parent/'manifest.json')}
    for rel,sha in manifest['outputs'].items(): checked(parent/rel,sha,hashes)
    a,current = load_sample(campaign,tag,'fit')
    for path,sha in manifest['inputs'].items():
        if path in current:
            if current[path] != sha: raise ValueError(f'Pooled training input changed: {path}')
        else: checked(Path(path),sha,hashes)
    hashes.update(current)
    table = read_tsv(parent/'folds.tsv')
    nfold = len({int(r['fold']) for r in table})
    if nfold < 2: raise ValueError('At least two saved training folds are required')
    folds = saved_folds(a,table,nfold)
    rows = matrix_rows(parent/'seed.dat'); terms = basis(rows)
    tr = read_tsv(parent/'terms.tsv'); pooled = read_tsv(parent/'pooled.tsv')
    codes = [''.join(map(str,t)) for t in terms]
    seeds,selected_rows = [],[]
    with np.load(parent/'coefficients.npz') as coeff:
        if not np.array_equal(coeff['terms'],terms): raise ValueError('Seed coefficient basis mismatch')
        matrices = {f:np.zeros((len(terms),3)) for f in range(nfold)}
        for j,target in enumerate(TARGETS):
            if set(selectors[target]) != {'alpha','l1'}: raise ValueError(f'Invalid seed selector for {target}')
            alpha,l1 = selectors[target]['alpha'],selectors[target]['l1']
            match = lambda r:r['target']==target and float(r['alpha'])==alpha and float(r['l1'])==l1
            if not any(match(r) and r['model']=='refit_svd' for r in pooled):
                raise ValueError(f'{target}: seed setting did not complete every fold')
            masks = []
            for f in range(nfold):
                selected = [r for r in tr if match(r) and int(r['fold'])==f]
                if len(selected) != len(codes) or {r['term'] for r in selected} != set(codes):
                    raise ValueError(f'{target}/fold {f}: missing or duplicate seed terms')
                lookup = {r['term']:r['selected'] for r in selected}
                if set(lookup.values())-{'True','False'} or lookup['00000'] != 'True':
                    raise ValueError('Invalid seed selection or missing constant')
                mask = tuple(i-1 for i,c in enumerate(codes) if i and lookup[c]=='True')
                masks.append(mask);selected_rows.extend(selected)
                key = f'{target}_a{alpha:g}_l{l1:g}_f{f}'
                values = coeff[key]
                if values.shape != (len(terms),) or not np.isfinite(values).all():
                    raise ValueError('Invalid seed coefficients')
                if any(values[i] != 0 for i,c in enumerate(codes) if lookup[c]=='False'):
                    raise ValueError('Seed coefficients disagree with selected terms')
                matrices[f][:,j] = values
            seeds.append(sorted(set(masks)))
    return a,rows,folds,seeds,matrices,selected_rows,manifest,hashes


def run(campaign,tag,name,source,cfg,selectors,check=False):
    from beam_diagnostics import report
    validate(cfg)
    if set(selectors) != set(TARGETS): raise ValueError('Expected one seed selector per target')
    a,seed,folds,seeds,seed_matrices,seed_terms,parent,hashes = inputs(campaign,tag,source,selectors)
    if check:
        print(f'OK: {len(folds):,} exact training cores, {len(set(folds))} folds; '+
              ', '.join(f'{t}: {len(s)} distinct seeds' for t,s in zip(TARGETS,seeds)),flush=True)
        print(f'Beam {cfg["beam"]}; threads {cfg["threads"]}; at most {cfg["steps"]} rounds per target. No outputs written.')
        return
    out = campaign/'06e_beam_search'/name
    if out.exists(): raise FileExistsError(f'Immutable output exists: {out}')
    terms = basis(seed)
    x,y,offset = problem(a,seed)
    cells = np.array([f'{z}:{d}' for z,d in zip(a['ztarT'],a['ndel'])])
    for folder in ('code','matrices','seeds','tsv','plots'): (out/folder).mkdir(parents=True)
    write_tsv(out/'seeds/terms.tsv',seed_terms)
    for f,c in seed_matrices.items(): write_matrix(out/f'seeds/fold{f}.dat',seed,terms,c)
    shutil.copy2(campaign/'06d_elastic_net'/source/'seed.dat',out/'seed.dat')
    shutil.copy2(campaign/'06d_elastic_net'/source/'folds.tsv',out/'tsv/folds.tsv')
    for filename in ('fit_beam.py','beam_search.py','beam_diagnostics.py','fit_elastic.py','elastic_net.py',
                     'elastic_refit.py','elastic_convergence.py','elastic_diagnostics.py','preallocated_svd.py',
                     'core_sample.py','spectrometer_config.py','spectrometer_profiles.def','run_beam.sh'):
        shutil.copy2(PROJECT/filename,out/'code'/filename)
    manifest = dict(schema='beam_v1',status='running',sample=tag,source=source,config=cfg,seeds=selectors,
                    inputs=hashes,rcond=parent['rcond'],training=len(folds),qa_min=10,
                    score='pooled mean physical-foil/delta MSE; adaptive development score',
                    versions={m:importlib.metadata.version(m) for m in ('numpy','scipy','threadpoolctl','uproot','matplotlib')})
    def save(): (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    try:
        prepared = Folds(x[:,1:],y,folds,cells,parent['rcond'],cfg['threads'])
        results = []
        candidate_fields = ['target','step','terms','score','condition','rank_min','valid','mask']
        with (out/'tsv/candidates.tsv').open('w') as fp:
            writer = csv.DictWriter(fp,fieldnames=candidate_fields,delimiter='\t');writer.writeheader()
            for j,target in enumerate(TARGETS):
                def record(r):
                    writer.writerow(dict(r,target=target,mask=' '.join(''.join(map(str,terms[i+1])) for i in r['mask'])))
                print(f'Search {target}',flush=True)
                result = search(lambda mask:prepared.evaluate(j,mask),seeds[j],x.shape[1]-1,cfg,record)
                fp.flush()
                results.append(result)
                write_tsv(out/'tsv/search.tsv',[dict(target=TARGETS[k],**r) for k,result in enumerate(results) for r in result['history']])
                # Preserve a compact, usable search result even if a later target fails.
                (out/'search.json').write_text(json.dumps(results,indent=2)+'\n')
        masks = dict(svd=[tuple(range(x.shape[1]-1))]*3,
                     start=[r['start']['mask'] for r in results],beam=[r['chosen']['mask'] for r in results])
        predictions,coefficients,fit_info,fold_info = {},{},{},[]
        with threadpool_limits(limits=cfg['threads']):
            for model,selected in masks.items():
                pred = []
                for j,mask in enumerate(selected):
                    values,stats = prepared.predict(j,mask);pred.append(values)
                    fold_info.extend(dict(r,model=model,target=TARGETS[j]) for r in stats)
                predictions[model] = np.column_stack(pred)
                coefficients[model],fit_info[model] = fit_all(x[:,1:],y,selected,parent['rcond'])
                write_matrix(out/f'matrices/{model}.dat',seed,terms,coefficients[model])
                reread = matrix_rows(out/f'matrices/{model}.dat')
                got = design(a,[t for c,t in reread]) @ np.array([c[:3] for c,t in reread])
                if not np.allclose(got,x @ coefficients[model]+offset,rtol=1e-10,atol=1e-11):
                    raise ValueError('Matrix export does not reproduce predictions')
        write_tsv(out/'tsv/conditioning.tsv',fold_info)
        write_tsv(out/'tsv/final_fit.tsv',[dict(model=m,**r) for m,info in fit_info.items() for r in info])
        write_tsv(out/'tsv/terms.tsv',[dict(target=TARGETS[j],term=''.join(map(str,t)),
                  **{m:bool(i == 0 or i-1 in selected[j]) for m,selected in masks.items()})
                  for j in range(3) for i,t in enumerate(terms)])
        np.savez_compressed(out/'model.npz',terms=terms,**coefficients)
        report(out,a,y,predictions,offset,folds,results,cells,manifest['qa_min'])
        manifest.update(status='complete',stops={t:r['stop'] for t,r in zip(TARGETS,results)},
                        selected_terms={t:r['chosen']['terms'] for t,r in zip(TARGETS,results)})
        manifest['outputs'] = {str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and p.name!='manifest.json'}
        save()
    except BaseException as exc:
        manifest.update(status='incomplete',error=str(exc));save();raise
    print(f'Beam search complete: {out}. Candidate matrices exported; no protected evaluation or replay installation.',flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path);p.add_argument('tag',nargs='?',default='equal15')
    p.add_argument('name',nargs='?',default='beam');p.add_argument('--source',default='pooled')
    for k in ('beam','threads','steps'): p.add_argument('--'+k,type=int)
    p.add_argument('--check',action='store_true')
    args = p.parse_args()
    for v in (args.tag,args.name,args.source):
        if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*',v): p.error('Invalid sample/output/source name')
    campaign = (args.campaign if args.campaign.is_absolute() else PROJECT/args.campaign).resolve()
    try:
        cfg = dict(DEFAULTS);selectors = dict(SEEDS)
        policy = campaign/'config/beam.json'
        if policy.exists():
            custom = json.loads(policy.read_text());selectors = custom.pop('seeds',selectors);cfg.update(custom)
        cfg.update({k:getattr(args,k) for k in ('beam','threads','steps') if getattr(args,k) is not None})
        run(campaign,args.tag,args.name,args.source,cfg,selectors,args.check)
    except (OSError,ValueError,KeyError,RuntimeError) as exc: p.exit(1,f'ERROR: {exc}\n')


if __name__ == '__main__': main()
