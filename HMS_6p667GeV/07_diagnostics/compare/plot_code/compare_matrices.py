#!/usr/bin/env python3
"""Frozen five-matrix comparison on protected and surplus events; never fit coefficients."""
import argparse
import csv
import json
import re
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from core_sample import PROJECT,digest,read_tsv,write_tsv
from fit_elastic import load_sample,checked,design
from elastic_diagnostics import residuals
from elastic_net import SCALES

MODELS = ('old','gmm','core','enet','beam')


def read_matrix(path, ignored=()):
    """Accept transport tables with or without separators; comments are never coefficients."""
    rows = []; remaining=list(ignored)
    for n,line in enumerate(path.read_text().splitlines(),1):
        if line.strip() in remaining:
            remaining.remove(line.strip()); continue
        line = line.split('!',1)[0].strip()
        if not line or line.startswith('#') or line.startswith('---'): continue
        fields = line.split()
        if len(fields)!=5 or not re.fullmatch(r'\d{5}',fields[4]):
            raise ValueError(f'Malformed matrix row: {path}:{n}')
        values = np.array([float(v.replace('D','E').replace('d','e')) for v in fields[:4]])
        if not np.isfinite(values).all(): raise ValueError(f'Nonfinite coefficient: {path}:{n}')
        rows.append((values,tuple(map(int,fields[4]))))
    if remaining: raise ValueError(f'Documented excluded matrix row was not found: {path}')
    if not rows or len({t for c,t in rows})!=len(rows): raise ValueError(f'Empty/duplicate matrix terms: {path}')
    return rows


def predict(a,matrices,threads=8,chunk=20000):
    """Use identical saved focal-plane/xtar inputs for every complete frozen matrix."""
    terms = sorted({t for rows in matrices.values() for c,t in rows})
    lookup = {t:i for i,t in enumerate(terms)}
    coeff = np.zeros((len(terms),3*len(MODELS)))
    for j,m in enumerate(MODELS):
        for c,t in matrices[m]: coeff[lookup[t],3*j:3*j+3]=c[:3]
    result = np.empty((len(a['entry']),3*len(MODELS)))
    with threadpool_limits(limits=threads):
        for start in range(0,len(result),chunk):
            part = {k:v[start:start+chunk] for k,v in a.items()}
            result[start:start+chunk] = design(part,terms) @ coeff
    if not np.isfinite(result).all(): raise ValueError('Nonfinite matrix predictions')
    return {m:result[:,3*j:3*j+3] for j,m in enumerate(MODELS)}


def historical_overlap(a,selected,table):
    known,seen = np.zeros(len(a['entry']),bool),np.zeros(len(a['entry']),bool)
    with np.load(selected,allow_pickle=False) as saved:
        for row in table:
            mask = a['rungroup']==row['rungroup'];key=f"{int(row['optics_id'])}_entry"
            if key not in saved: continue
            ids = saved[key]
            if ids.dtype.kind not in 'iu':
                if not np.isfinite(ids).all() or np.any(ids!=np.floor(ids)) or np.any(np.abs(ids)>2**53):
                    raise ValueError('Historical event IDs cannot be represented exactly')
                ids = ids.astype(np.int64)
            if len(np.unique(ids)) != len(ids): raise ValueError('Duplicate historical fit IDs')
            known[mask]=True;seen[mask]=np.isin(a['entry'][mask],ids)
    return known,seen


def offset_rows(rungroups,policy):
    """External corrections are explicit physical-unit additions, never fitted to evaluation data."""
    values=policy.get('offsets',{});notes=policy.get('offset_sources',{})
    if set(values)!=set(MODELS):raise ValueError('Specify external offsets for all five matrices')
    rows=[]
    for m in MODELS:
        value=values[m]
        if value is None:raise ValueError(f'{m}: external offsets are unresolved; supply historical replay settings')
        if not isinstance(notes.get(m),str) or not notes[m].strip():raise ValueError(f'{m}: document the offset source')
        for rg in sorted(set(rungroups)):
            v=value.get(rg) if isinstance(value,dict) else value
            if not isinstance(v,list) or len(v)!=3 or any(type(n) not in (int,float) or not np.isfinite(n) for n in v):
                raise ValueError(f'{m}/{rg}: require [xptar mrad, ytar cm, yptar mrad] external offsets')
            rows.append(dict(model=m,rungroup=rg,xptar_mrad=v[0],ytar_cm=v[1],yptar_mrad=v[2],source=notes[m]))
    return rows


def apply_offsets(predictions,a,rows):
    for r in rows:
        mask=a['rungroup']==r['rungroup']
        predictions[r['model']][mask]+=np.array([r['xptar_mrad'],r['ytar_cm'],r['yptar_mrad']])/SCALES


def prepare(campaign,tag,source,policy):
    src=campaign/'06e_beam_search'/source
    beam=json.loads((src/'manifest.json').read_text())
    if (beam.get('schema'),beam.get('status'),beam.get('sample'))!=('beam_v1','complete',tag):
        raise ValueError('Requires a completed beam result for the same sample')
    hashes={str(src/'manifest.json'):digest(src/'manifest.json')}
    required=['seed.dat','matrices/svd.dat','matrices/start.dat','matrices/beam.dat','tsv/folds.tsv','tsv/final_fit.tsv']
    for rel in required: checked(src/rel,beam['outputs'][rel],hashes)
    # Large candidate logs are not required for frozen matrix evaluation.
    paths=dict(old=src/'seed.dat',gmm=campaign/policy['gmm'],core=src/'matrices/svd.dat',
               enet=src/'matrices/start.dat',beam=src/'matrices/beam.dat')
    if policy.get('old'): paths['old']=campaign/policy['old']
    exclusions=policy.get('exclude_rows',{})
    if set(exclusions)-set(MODELS):raise ValueError('Unknown matrix row exclusion')
    for m, item in exclusions.items():
        if not policy.get(m+'_sha256') or not item.get('reason') or not item.get('rows'):
            raise ValueError('Row exclusions require a pinned matrix checksum, exact rows and reason')
    matrices={m:read_matrix(p,exclusions.get(m,{}).get('rows',[])) for m,p in paths.items()}
    for m,p in paths.items():
        sha=digest(p)
        if policy.get(m+'_sha256') and sha!=policy[m+'_sha256']:
            raise ValueError(f'{m} matrix differs from documented historical checksum')
        hashes[str(p)]=sha
    train={(r['rungroup'],int(r['entry'])) for r in read_tsv(src/'tsv/folds.tsv')}
    corrections=offset_rows([r for r,e in train],policy)
    samples=[];pools=[]
    for sample in ('holdout','surplus'):
        a,h=load_sample(campaign,tag,sample);hashes.update(h)
        if train & set(zip(a['rungroup'],map(int,a['entry']))): raise ValueError('Evaluation/training overlap')
        samples.append(a)
        names=np.full(len(a['entry']),'surplus_core',dtype='U32') if sample=='surplus' else np.where(
            a['quality']==2,'protected_core',np.where(a['quality']==1,'protected_noncore','protected_unsupported'))
        names=np.where(a['excluded'],'blocked',names);pools.append(names)
    a={k:np.concatenate([s[k] for s in samples]) for k in samples[0]};pool=np.concatenate(pools)
    corrections=offset_rows(a['rungroup'],policy)
    ids=list(zip(a['rungroup'],map(int,a['entry'])))
    if len(set(ids))!=len(ids):raise ValueError('Duplicate events across evaluation pools')
    allocation=campaign/'05c_core_sample'/tag
    suffix=f'/05c_core_sample/{tag}/manifest.json'
    expected={h for p,h in beam['inputs'].items() if p.endswith(suffix)}
    if expected!={digest(allocation/'manifest.json')}: raise ValueError('Allocation differs from beam training allocation')
    tables=list(allocation.glob('rungroups_*_inputs.tsv'))
    if len(tables)!=1:raise ValueError('Expected one frozen campaign table')
    selected=campaign/policy['gmm_ids']
    provenance=campaign/policy['gmm_provenance']
    history=json.loads(provenance.read_text())
    for p in (selected,provenance):hashes[str(p)]=digest(p)
    known,seen=historical_overlap(a,selected,read_tsv(tables[0]))
    return dict(source=src,beam=beam,paths=paths,matrices=matrices,data=a,pool=pool,
                known=known,seen=seen,history=history,hashes=hashes,offsets=corrections)


def run(campaign,tag,name,source,policy,threads=8,check=False):
    from comparison_plots import report
    if threads<1:raise ValueError('threads must be positive')
    out=campaign/'07_diagnostics'/name
    if not check and out.exists():raise FileExistsError(f'Immutable comparison exists: {out}')
    job=prepare(campaign,tag,source,policy)
    counts={p:int(sum(job['pool']==p)) for p in np.unique(job['pool'])}
    print('Evaluation pools:',counts,flush=True)
    if check:
        print('Five frozen matrices and exact memberships verified. No residuals evaluated; no outputs written.');return
    for folder in ('plots','tsv','matrices','code'):(out/folder).mkdir(parents=True)
    manifest=dict(schema='matrix_comparison_v1',status='running',sample=tag,source=source,policy=policy,
                  inputs=job['hashes'],counts=counts,threads=threads,
                  historical_overlap='reconstructed GMM membership; old replay training membership unknown')
    def save():(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    try:
        for m,p in job['paths'].items():shutil.copy2(p,out/f'matrices/{m}.dat')
        for file in ('compare_matrices.py','comparison_plots.py','fit_elastic.py','elastic_net.py',
                     'elastic_diagnostics.py','core_sample.py','preallocated_svd.py','spectrometer_config.py',
                     'spectrometer_profiles.def','run_compare.sh'):shutil.copy2(PROJECT/file,out/'code'/file)
        shutil.copy2(job['source']/'tsv/final_fit.tsv',out/'tsv/training_conditioning.tsv')
        write_tsv(out/'tsv/offsets.tsv',job['offsets'])
        a=job['data'];predictions=predict(a,job['matrices'],threads)
        apply_offsets(predictions,a,job['offsets'])
        truth=np.column_stack((a['xptarT'],a['ytarT']/100.,a['yptarT']))
        rr={m:residuals(a,truth,p,np.zeros_like(truth)) for m,p in predictions.items()}
        report(out,a,job['pool'],rr,job['known'],job['seen'],job['history'])
        manifest['status']='complete'
        manifest['outputs']={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and p.name!='manifest.json'}
        save()
    except BaseException as exc:
        manifest.update(status='incomplete',error=str(exc));save();raise
    print(f'Frozen comparison complete: {out}. Coefficients were not refitted.',flush=True)


def replot(campaign,tag,name):
    """Redraw a completed comparison from saved residuals, without prediction or statistics."""
    from comparison_plots import plots,TARGETS
    out=campaign/'07_diagnostics'/name
    manifest_path=out/'manifest.json'
    manifest=json.loads(manifest_path.read_text())
    if (manifest.get('schema'),manifest.get('status'),manifest.get('sample'))!=('matrix_comparison_v1','complete',tag):
        raise ValueError('Replot requires a completed comparison for the same sample')
    hashes={}
    for rel in ('residuals.npz','tsv/residuals.tsv'):
        if not (out/rel).exists():raise FileNotFoundError(f'Replot needs saved {out/rel}; run on the machine holding the original comparison')
        checked(out/rel,manifest['outputs'][rel],hashes)
    print('Saved residuals and tables verified; redrawing only. No matrix evaluation or grouped statistics.',flush=True)
    with np.load(out/'residuals.npz',allow_pickle=False) as saved:
        if tuple(saved['targets'])!=TARGETS:raise ValueError('Unexpected saved residual target order')
        pool=saved['pool'];rr={m:saved[m] for m in MODELS}
    if pool.ndim!=1 or any(v.shape!=(len(pool),len(TARGETS)) or not np.isfinite(v).all() for v in rr.values()):
        raise ValueError('Invalid saved residual arrays')
    with (out/'tsv/residuals.tsv').open(newline='') as f:
        details=[r for r in csv.DictReader(f,delimiter='\t') if r['level']=='foil_delta']
    for r in details:
        for k in ('zfoil','rms'):r[k]=float(r[k])
        for k in ('ndel','n'):r[k]=int(r[k])
    # Generate everything before replacing existing figures; retain the original numerical code snapshot.
    with tempfile.TemporaryDirectory(prefix='.replot-',dir=out) as tmp:
        stage=Path(tmp);plots(stage,pool,rr,details)
        (stage/'plot_code').mkdir()
        for file in ('compare_matrices.py','comparison_plots.py'):shutil.copy2(PROJECT/file,stage/'plot_code'/file)
        report_path=out/'MATRIX_COMPARISON.md'
        if report_path.exists():
            text=report_path.read_text()
            text=re.sub(r'(?ms)^plots/\*_center\.png.*?\n\n',
                        'Read [plots/README.md](plots/README.md) for the revised figures. Plotting changes do not change the saved residuals or numerical tables.\n\n',text)
            (stage/report_path.name).write_text(text)
        updated={str(p.relative_to(stage)):digest(p) for p in stage.rglob('*') if p.is_file()}
        for rel in updated:
            (out/rel).parent.mkdir(exist_ok=True)
            (stage/rel).replace(out/rel)
    manifest['outputs'].update(updated)
    manifest['replot']=dict(time=datetime.now(timezone.utc).isoformat(),inputs=hashes,
                            code='plot_code',note='Presentation only; residuals, tables, matrices and offsets unchanged')
    pending=out/'manifest.json.tmp';pending.write_text(json.dumps(manifest,indent=2)+'\n');pending.replace(manifest_path)
    print(f'Plots refreshed: {out}/plots. Numerical results unchanged.',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path);p.add_argument('tag',nargs='?',default='equal15')
    p.add_argument('name',nargs='?',default='compare');p.add_argument('--source',default='beam10')
    p.add_argument('--threads',type=int,default=8);p.add_argument('--check',action='store_true')
    p.add_argument('--replot',action='store_true',help='refresh plots from an existing comparison; do not evaluate matrices again')
    a=p.parse_args()
    if a.replot and a.check:p.error('--replot and --check are separate modes')
    for v in (a.tag,a.name,a.source):
        if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*',v):p.error('Invalid tag/name/source')
    campaign=(a.campaign if a.campaign.is_absolute() else PROJECT/a.campaign).resolve()
    try:
        if a.replot:
            replot(campaign,a.tag,a.name);return
        policy=json.loads((campaign/'config/comparison.json').read_text())
        allowed={'old','gmm','old_sha256','gmm_sha256','gmm_ids','gmm_provenance','offsets','offset_sources','exclude_rows'}
        if set(policy)-allowed or not {'gmm','gmm_ids','gmm_provenance'}<=set(policy):raise ValueError('Invalid comparison policy')
        run(campaign,a.tag,a.name,a.source,policy,a.threads,a.check)
    except (OSError,ValueError,KeyError,RuntimeError) as exc:p.exit(1,f'ERROR: {exc}\n')


if __name__=='__main__':main()
