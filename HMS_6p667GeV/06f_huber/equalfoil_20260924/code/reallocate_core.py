#!/usr/bin/env python3
"""Preview or reallocate a frozen core tag; never estimate density or run SVD."""
import argparse
import importlib.metadata
import sys
import json
import re
import shutil
import tempfile
from pathlib import Path
import numpy as np
import uproot
from core_sample import PROJECT, digest, read_tsv, write_tsv, groups, ranks
from core_balance import LEAF, allocate_counts, validate
from spectrometer_config import from_campaign, run_metadata


def checked(path, expected):
    if not path.is_file():
        raise FileNotFoundError(path)
    if digest(path) != expected:
        raise ValueError(f'Frozen input hash mismatch: {path}')
    return path


def load_inputs(campaign, parent):
    source = campaign / '05c_core_sample' / parent
    manifest = json.loads((source/'manifest.json').read_text())
    outputs = manifest['outputs']
    for rel in ('config.json', 'tsv/allocation.tsv', 'tsv/counts.tsv'):
        checked(source/rel, outputs[rel])
    names = [k for k in outputs if k.startswith('rungroups_') and k.endswith('_inputs.tsv')]
    if len(names) != 1:
        raise ValueError('Parent must name exactly one campaign table')
    # Locally mirrored tags may omit the saved copy. Only a hash-identical copy is allowed.
    table = source/names[0]
    if not table.exists():
        table = campaign/'config'/names[0]
    checked(table, outputs[names[0]])
    settings = read_tsv(table)
    if len({r['rungroup'] for r in settings}) != len(settings) or len({r['optics_id'] for r in settings}) != len(settings):
        raise ValueError('Duplicate rungroup or optics identity')
    metadata = source/'metadata/optics.dat'
    if metadata.exists():
        checked(metadata, outputs['metadata/optics.dat'])
    else:
        metadata = campaign/'06b_svd_fit/diagnostics/balance/optics.dat'
        prov = json.loads((metadata.parent/'provenance.json').read_text())
        expected, = [v for k,v in prov['input_sha256'].items() if k.endswith('DATfiles/list_of_optics_run.dat')]
        checked(metadata, expected)
    meta = {r['rungroup']: run_metadata(int(r['optics_id']), metadata) for r in settings}
    edges = list(meta.values())[0]['edges']
    if any(m['edges'] != edges for m in meta.values()):
        raise ValueError('Incompatible delta boundaries: frozen settings cannot be rebinned')
    spec = from_campaign(campaign)
    if spec.name != 'HMS':
        raise ValueError('Frozen balance adapter currently verified for HMS only')
    required = [-8., -3., 0., 3., 8.] if campaign.name == 'HMS_6p667GeV' else sorted({f for m in meta.values() for f in m['foils']})
    maskpath = campaign/'config/sieve_mask.json'
    mask = json.loads(maskpath.read_text())
    blocked = {tuple(p) for p in mask['blocked']}
    if not all(len(p)==2 and all(type(i) is int for i in p) and 0<=p[0]<spec.nx and 0<=p[1]<spec.ny for p in blocked):
        raise ValueError('Invalid campaign geometry mask')
    if campaign.name=='HMS_6p667GeV' and not {(2,3),(5,5)} <= blocked:
        raise ValueError('Required confirmed HMS blocked positions missing')
    cfg = dict(manifest['config'])
    cfg.update(balance=dict(foil='equal',delta=1.5,hole=1.5), qa_min=10)
    policy = campaign/'config/core_balance.json'
    if policy.exists():
        changes = json.loads(policy.read_text())
        if changes.keys()-{'balance','qa_min','fit_max'}:
            raise ValueError('Frozen reallocation only changes balance, qa_min, fit_max; seed/holdout stay frozen')
        cfg.update(changes)
    validate(cfg)
    alloc = {}
    for r in read_tsv(source/'tsv/allocation.tsv'):
        key = (r['rungroup'],float(r['zfoil']),float(r['delta_low']),float(r['delta_high']),int(r['xscol']),int(r['yscol']))
        if key in alloc or int(r['available']) != int(r['fit'])+int(r['surplus']):
            raise ValueError('Invalid original allocation summary')
        alloc[key] = int(r['available'])
    rows = []
    for r in read_tsv(source/'tsv/counts.tsv'):
        name = r['rungroup']
        m = meta[name]
        foil, nd = int(r['foil']), int(r['ndel'])
        if not 0 <= foil < len(m['foils']) or not 0 <= nd < len(edges)-1:
            raise ValueError('Invalid local foil or delta label')
        key = (name,m['foils'][foil],edges[nd],edges[nd+1],int(r['xscol']),int(r['yscol']))
        original = int(r['fit'])+int(r['surplus'])
        if original != alloc.pop(key, 0):
            raise ValueError(f'Allocation/quality counts or physical metadata disagree: {key}')
        if not (0<=key[-2]<spec.nx and 0<=key[-1]<spec.ny):
            raise ValueError(f'Out-of-geometry label: {key}')
        excluded = key[-2:] in blocked
        row = dict(zip(LEAF,key))
        row.update(foil=foil, ndel=nd, original_core=original, available=0 if excluded else original,
                   excluded_core=(original+int(r.get('excluded_core',0))) if excluded else 0, blocked=excluded,
                   total=int(r['total']), holdout=int(r['holdout']), holdout_core=int(r['holdout_core']),
                   unsupported=int(r['unsupported']), noncore_development=int(r['total'])-int(r['holdout'])-original-int(r.get('excluded_core',0)),
                   geometry_excluded_development=int(r['total'])-int(r['holdout']) if excluded else 0,
                   pool=cfg.get('pools',{}).get(name,name))
        rows.append(row)
    if alloc:
        raise ValueError('Original allocation leaves absent from quality counts')
    inputs = dict(source=source, manifest=manifest, settings=settings, table=table, metadata=metadata,
                  meta=meta, edges=edges, required=required, maskpath=maskpath, blocked=blocked, cfg=cfg, rows=rows)
    return inputs


def event_allocate(all_data, inputs, result):
    """Apply already computed quotas once and verify against frozen summaries."""
    chosen, audit = [], []
    quotas = {tuple(r[k] for k in LEAF):r for r in result['leaves']}
    actual_counts = {}
    for name, data in sorted(all_data.items()):
        n = len(data['entry'])
        if len(np.unique(data['entry'])) != n:
            raise ValueError(f'Duplicate event identity in {name}')
        if not np.issubdtype(data['entry'].dtype,np.integer):
            raise ValueError('entry must retain integer identity')
        if not np.isin(data['quality'],[0,1,2]).all() or not np.isin(data['sample'],[0,1,2,3]).all():
            raise ValueError('Invalid frozen event categories')
        m = inputs['meta'][name]
        for (foil,nd), idx in groups(data,['foil','ndel']).items():
            if not 0<=foil<len(m['foils']) or not 0<=nd<len(m['edges'])-1:
                raise ValueError('Invalid event foil/delta indices')
            for field, value in [('zfoil',m['foils'][foil]),('delta_low',m['edges'][nd]),('delta_high',m['edges'][nd+1])]:
                if not np.all(data[field][idx]==value):
                    raise ValueError(f'Event metadata mismatch: {name} {field}')
        protected = data['sample']==2
        oldsample = data['sample'].copy()
        excluded = np.array([(int(x),int(y)) in inputs['blocked'] for x,y in zip(data['xscol'],data['yscol'])])
        if np.any(np.isin(oldsample,[1,3]) & (data['quality']!=2)):
            raise ValueError('Frozen fit/surplus contains noncore events')
        if np.any((data['quality']==2)&~protected&~np.isin(oldsample,[1,3])&~excluded):
            raise ValueError('Frozen accepted core missing from old fit/surplus')
        eligible = (data['quality']==2)&~protected&~excluded
        data['sample'][eligible]=3
        data['sample'][excluded&~protected]=0
        data['balance_excluded']=excluded.astype(np.int32)
        # Validate all observed populations, including zero-core and excluded labels.
        for key, idx in groups(data, list(LEAF[1:])).items():
            full = (name,*key)
            if full not in quotas:
                raise ValueError(f'Unexpected event leaf: {full}')
            row = quotas[full]
            av = idx[eligible[idx]]
            if len(idx)!=row['total'] or int(protected[idx].sum())!=row['holdout'] or len(av)!=row['available']:
                raise ValueError(f'Frozen counts/event membership disagreement: {full}')
            if 'holdout_core' in row:
                if int(np.sum(protected[idx] & (data['quality'][idx]==2)))!=row['holdout_core'] or int(np.sum(data['quality'][idx]==0))!=row['unsupported']:
                    raise ValueError(f'Frozen quality categories disagree: {full}')
            actual_counts[full]=len(av)
            order=np.lexsort((data['entry'][av],ranks(inputs['cfg']['seed'],name,data['entry'][av],'fit')))
            selected=av[order[:row['quota']]]
            data['sample'][selected]=1
            for i in selected:
                chosen.append(dict(rungroup=name,entry=int(data['entry'][i]),**{k:row[k] for k in LEAF[1:]},foil=row['foil'],ndel=row['ndel']))
            row['actual_selected']=len(selected)
        for i in np.flatnonzero(excluded):
            audit.append(dict(rungroup=name,entry=int(data['entry'][i]),reason='geometry_blocked',
                              protected=bool(protected[i]),quality=int(data['quality'][i]),old_sample=int(oldsample[i]),sample=int(data['sample'][i])))
        assert np.array_equal(data['sample']==2,protected)
        assert int((data['sample']==1).sum())==sum(r['quota'] for r in result['leaves'] if r['rungroup']==name)
    if set(actual_counts)!=set(quotas):
        raise ValueError('Missing frozen event populations')
    chosen.sort(key=lambda r:(r['rungroup'],r['entry']))
    assert len(chosen)==result['training_total']
    return chosen,audit


def run(campaign,parent,tag,preview=False):
    inp=load_inputs(campaign,parent)
    result=allocate_counts(inp['rows'],inp['cfg'],inp['required'],zip(inp['edges'],inp['edges'][1:]))
    source=inp['source']
    required=[source/rel for rel in inp['manifest']['outputs'] if rel.startswith(('root/','models/','code/')) or rel.endswith('_excluded.tsv')]
    missing=[str(p) for p in required if not p.is_file()]
    origin_tables=[Path(k) for k in inp['manifest'].get('inputs',{}) if Path(k).name==inp['table'].name and '/config/' in k]
    origin_missing=[str(origin_tables[0].parent.parent/'05c_core_sample'/parent/p.relative_to(source)) for p in required if not p.is_file()] if origin_tables else []
    replay_missing=[r['rootfile'] for r in inp['settings'] if not Path(r['rootfile']).is_file()]
    if not preview and missing and result['training_total']:
        raise FileNotFoundError('Frozen reallocation requires these exact parent files (no density rerun):\n'+'\n'.join(missing))
    target=campaign/('05d_core_balance' if preview else '05c_core_sample')/tag
    if target.exists():
        raise FileExistsError(f'Immutable output exists: {target}; choose a new tag')
    target.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.balance_',dir=target.parent))
    try:
        for sub in ('tsv','code','metadata','parent'):
            (stage/sub).mkdir()
        for name in ('core_balance.py','reallocate_core.py','balance_diagnostics.py','core_sample.py','spectrometer_config.py','spectrometer_profiles.def','build_core_fit.py','preallocated_svd.py','preallocated_sample.h','fit_opt_matrix_gmm.C','make_fit_ntuple_from_gmm.C','spectrometer_config.h','spectrometer_root.h'):
            shutil.copy2(PROJECT/name,stage/'code'/name)
        (stage/'code/tests').mkdir()
        for name in ('test_core_balance.py','test_frozen_reallocation.py','test_preallocated_svd.py','test_core_sample.py'):
            shutil.copy2(PROJECT/'tests'/name,stage/'code/tests'/name)
        shutil.copy2(PROJECT/'CORE_SAMPLE.md',stage/'code/CORE_SAMPLE.md')
        shutil.copy2(source/'manifest.json',stage/'parent/manifest.json')
        shutil.copy2(source/'config.json',stage/'parent/config.json')
        shutil.copy2(inp['table'],stage/inp['table'].name)
        shutil.copy2(inp['metadata'],stage/'metadata/optics.dat')
        shutil.copy2(inp['maskpath'],stage/'metadata/sieve_mask.json')
        for p in (source/'tsv').glob('*excluded.tsv'):
            checked(p,inp['manifest']['outputs'][str(p.relative_to(source))])
            shutil.copy2(p,stage/'tsv'/p.name)
        (stage/'config.json').write_text(json.dumps(inp['cfg'],indent=2)+'\n')
        (stage/'INPUT_AVAILABILITY.md').write_text('# Frozen input availability\n\nReallocation needs the missing frozen files below. Replay files are needed only for TFit export.\n\n## Frozen files missing locally\n\n'+''.join(f'- `{p}`\n' for p in missing)+'\n## Corresponding originating frozen paths (ifarm, inferred from parent input manifest)\n\n'+''.join(f'- `{p}`\n' for p in origin_missing)+'\n## Replay files missing locally\n\n'+''.join(f'- `{p}`\n' for p in replay_missing))
        data=None
        if not preview and result['training_total']:
            for p in required:
                checked(p,inp['manifest']['outputs'][str(p.relative_to(source))])
            data={}
            for r in inp['settings']:
                with uproot.open(source/f"root/CoreSample_{r['rungroup']}.root") as f:
                    data[r['rungroup']]=f['CoreSample'].arrays(library='np')
            ids,audit=event_allocate(data,inp,result)
            write_tsv(stage/'tsv/selected_ids.tsv',ids)
            write_tsv(stage/'tsv/geometry_exclusions.tsv',audit)
            (stage/'root').mkdir()
            for name,a in data.items():
                with uproot.recreate(stage/f'root/CoreSample_{name}.root') as f:
                    f.mktree('CoreSample',a)
            shutil.copytree(source/'models',stage/'models')
            # Preserve the established frozen-tag table interface for subsequent reallocations.
            allocation=[]; category_counts=[]
            for r in result['leaves']:
                allocation.append(dict(pool=r['pool'],**{k:r[k] for k in LEAF},available=r['available'],fit=r['quota'],surplus=r['surplus']))
            write_tsv(stage/'tsv/allocation.tsv',allocation)
            for name,a in sorted(data.items()):
                for (nf,nd,x,y),idx in groups(a,['foil','ndel','xscol','yscol']).items():
                    category_counts.append(dict(rungroup=name,foil=nf,ndel=nd,xscol=x,yscol=y,total=len(idx),
                        fit=int(np.sum(a['sample'][idx]==1)),holdout=int(np.sum(a['sample'][idx]==2)),
                        holdout_core=int(np.sum((a['sample'][idx]==2)&(a['quality'][idx]==2))),
                        surplus=int(np.sum(a['sample'][idx]==3)),unsupported=int(np.sum(a['quality'][idx]==0)),
                        excluded_core=int(np.sum((a['balance_excluded'][idx]!=0)&(a['sample'][idx]!=2)&(a['quality'][idx]==2)))))
            write_tsv(stage/'tsv/counts.tsv',category_counts)
        for row in result['leaves']:
            row.setdefault('actual_selected','')
        for level in ('foils','deltas'):
            for row in result[level]:
                row['actual_selected']=row['quota'] if data is not None else ''
        result['mode']='counts_preview' if preview else ('event_allocation' if data is not None else 'zero_budget')
        result['config']=inp['cfg']
        (stage/'allocation.json').write_text(json.dumps(result,indent=2)+'\n')
        for level in ('leaves','deltas','foils'):
            write_tsv(stage/'tsv'/f'{level}.tsv',result[level])
        from balance_diagnostics import diagnostics
        diagnostics(stage,campaign,inp,result,data)
        manifest=dict(python=sys.version,versions={m:importlib.metadata.version(m) for m in ('numpy','uproot','matplotlib','scipy','scikit-image')},schema='core_balance_v1',config=inp['cfg'],mode=result['mode'],
                      status=result['status'],totals=dict(fit=result['training_total']),
                      parent_manifest=digest(source/'manifest.json'),
                      inputs={str(p):digest(p) for p in (inp['table'],inp['metadata'],inp['maskpath'],source/'tsv/allocation.tsv',source/'tsv/counts.tsv')},
                      grouping={r['rungroup']:inp['cfg'].get('pools',{}).get(r['rungroup'],r['rungroup']) for r in inp['settings']},
                      metadata_origin=str(inp['metadata']),missing_frozen=missing,missing_replay=replay_missing,
                      inactive_legacy_controls=result['inactive_legacy_controls'])
        manifest['outputs']={str(p.relative_to(stage)):digest(p) for p in sorted(stage.rglob('*')) if p.is_file()}
        (stage/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        stage.rename(target)
        print(f"{result['mode']}: {target}\nBudget: {result['training_total']} = {len(inp['required'])} x {result['foil_budget']}")
        if not result['training_total']:
            print(f"STOP: zero campaign budget; limiting foils {result['zero_capacity_foils']}; constraint {result['budget_constraint']}; see foil/delta/hole constraints. No export or solver use.")
        return target
    except BaseException:
        shutil.rmtree(stage)
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path)
    p.add_argument('parent')
    p.add_argument('tag')
    p.add_argument('--preview',action='store_true')
    args=p.parse_args()
    if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',v) for v in (args.parent,args.tag)):
        p.error('Invalid tag')
    try:
        run(args.campaign.resolve(),args.parent,args.tag,args.preview)
    except (ValueError,OSError,KeyError) as e:
        p.exit(1,f'ERROR: {e}\n')
