#!/usr/bin/env python3
"""Fixed central-hole cores, configured GMM matrix, existing HMS afterburner."""
import argparse,csv,json,subprocess,shutil
from pathlib import Path
import numpy as np
import uproot
from core_sample import digest
from spectrometer_config import from_campaign,run_metadata
from compare_matrices import read_matrix
from sieve_afterburner import reconstruct

PROJECT=Path(__file__).resolve().parent

def output_tag(policy):
    offsets=policy['offsets']['gmm']
    status='zero' if all(v==0 for v in offsets) else '_'.join(f'{v:g}'.replace('-','m').replace('.','p') for v in offsets)
    return Path(policy['gmm']).stem+'__external_offsets_'+status+'__fitted_constants_included'

def name_outputs(out,policy):
    tag=output_tag(policy)
    for path in list(out.iterdir()):
        if path.is_file() and '__external_offsets_' not in path.stem:
            path.rename(path.with_name(path.stem+'__'+tag+path.suffix))
    return tag


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign');p.add_argument('--rungroup',required=True)
    args=p.parse_args();campaign=Path(args.campaign).resolve()
    spec=from_campaign(campaign)
    if spec.name!='HMS':raise ValueError('Existing afterburner supports HMS only')
    cfg=json.loads((campaign/'config/geometry_study.json').read_text())
    policy=json.loads((campaign/'config/comparison.json').read_text())
    tables=list((campaign/'config').glob('rungroups_*_inputs.tsv'))
    if len(tables)!=1:raise ValueError('Expected one rungroup table')
    row=next(r for r in csv.DictReader(tables[0].open(),delimiter='\t') if r['rungroup']==args.rungroup)
    angle=float(row.get('angle_deg') or row['hms_angle_deg']);foil=cfg.get('foil',0.)
    if foil not in run_metadata(int(row['optics_id']),PROJECT/'DATfiles/list_of_optics_run.dat')['foils']:
        raise ValueError('Foil absent from group')
    core=campaign/cfg['core_file'].format(rungroup=args.rungroup)
    with uproot.open(core) as f:a=f['CoreSample'].arrays(library='np')
    mask=(a['core_keep']==1)&(a['xscol']==spec.nx//2)&(a['yscol']==spec.ny//2)&(np.abs(a['zfoil']-foil)<1e-8)
    mask&=(a['delta']>spec.delta_min)&(a['delta']<spec.delta_max)&(a['delta']>=a['delta_low'])&(a['delta']<a['delta_high'])
    mask&=(a['sumnpe']>2)&(a['etracknorm']>.65)&(np.abs(a['reactz']-foil)<cfg.get('foil_width',2.))
    a={k:v[mask] for k,v in a.items()};a['ry']=a['reacty'];a['rz']=a['reactz']
    dense=np.zeros(len(a['entry']),bool)
    for ndel in np.unique(a['ndel']):
        ids=np.flatnonzero(a['ndel']==ndel);scores=np.sort(a['core_score'][ids])[::-1]
        dense[ids]=a['core_score'][ids]>=scores[(len(ids)-1)//2]
    out=campaign/cfg['output']/f'foil_{foil:g}cm'/'gmm_afterburner'/args.rungroup
    out.mkdir(parents=True,exist_ok=True)
    manifest=dict(rungroup=args.rungroup,core=str(core),source_replay=row['rootfile'],foil=foil,angle=angle,
                  n=len(dense),n_dense=int(dense.sum()),selection='Frozen core_keep and core_score; original PID/delta/vertex',matrices={})
    targets={}
    for model in ('old','gmm'):
        path=campaign/policy[model]
        actual=digest(path)
        if actual!=policy[model+'_sha256']:raise ValueError(f'{model} matrix checksum mismatch')
        excluded=policy.get('exclude_rows',{}).get(model,{}).get('rows',[])
        rows=read_matrix(path,excluded)
        offsets=policy['offsets'][model]
        xy,it,active,target,xtar=reconstruct(a,rows,angle,offsets,return_target=True)
        targets[model]=target
        manifest['matrices'][model]=dict(path=str(path),sha256=actual,excluded_rows=excluded,offsets=offsets,
                                        unconverged=int(active.sum()),max_iterations=int(it.max()),
                                        fitted_zero_order_coefficients=[c.tolist() for c,t in rows if t==(0,0,0,0,0)])
    old=targets['old'];new=targets['gmm']
    manifest['closure']=dict(xptar_max_abs_mrad=float(np.max(np.abs(1000*(old[:,0]-a['xptar'])))),
                            ytar_max_abs_cm=float(np.max(np.abs(100*old[:,1]-a['ytar']))))
    print('Original-matrix closure:',manifest['closure'],flush=True)
    if manifest['closure']['xptar_max_abs_mrad']>1e-5 or manifest['closure']['ytar_max_abs_cm']>1e-6:
        raise ValueError('Original matrix does not reproduce saved values; inspect closure')
    prediction=out/'predictions.tsv'
    with prediction.open('w') as f:
        for entry,target in zip(a['entry'],new):f.write(f'{int(entry)}\t{target[0]:.17g}\t{100*target[1]:.17g}\n')
    np.savez_compressed(out/'event_comparison.npz',entry=a['entry'],dense=dense,reacty=a['reacty'],
                        xptar_old=a['xptar'],xptar_gmm=new[:,0],ytar_old=a['ytar'],ytar_gmm=100*new[:,1])
    values=[str(core),row['rootfile'],str(out),spec.name,int(row['optics_id']),angle,foil,
            f'{spec.nx//2},{spec.ny//2}',cfg.get('min_events',30),cfg.get('fp_fraction',.5),cfg.get('foil_width',2.),True,str(prediction)]
    expression=str(PROJECT/'diagnostics/validation/geometry/core_geometry.C')+'('+','.join(json.dumps(v) for v in values)+')'
    with (out/'terminal.txt').open('w') as log:
        subprocess.run(['root','-l','-b','-q',expression],stdout=log,stderr=subprocess.STDOUT,check=True)
    manifest['filename_tag']=output_tag(policy)
    manifest['offset_convention']='External additions [xptar mrad, ytar cm, yptar mrad]; active fitted 00000 coefficients retained'
    shutil.copy2(campaign/policy['gmm'],out/'matrix.dat')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Saved',out,flush=True)
    for r in csv.DictReader((out/'slopes.tsv').open(),delimiter='\t'):
        if r['ndel']=='-1':print(r['subset'],r['n'],r['slope_mrad_per_cm'],r['slope_error'])
    print((out/'x4_y4_ytar.tsv').read_text())
    name_outputs(out,policy)

if __name__=='__main__':main()
