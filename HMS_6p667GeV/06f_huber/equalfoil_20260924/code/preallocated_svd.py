#!/usr/bin/env python3
"""Verify an immutable core build and run explicit SVD sample acceptance (default: no solve)."""
import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path
import numpy as np
import uproot
from core_sample import PROJECT,digest,read_tsv,write_tsv


def membership_manifest(path,ids,settings,allocation_hash):
    optics={r['rungroup']:int(r['optics_id']) for r in settings}
    with path.open('w') as out:
        out.write('# core_preallocated_v1\n# allocation_sha256='+allocation_hash+'\n')
        out.write('rungroup\tentry\toptics_id\tfoil\tndel\tzfoil\tdelta_low\tdelta_high\txscol\tyscol\n')
        for r in sorted(ids,key=lambda r:(r['rungroup'],int(r['entry']))):
            values=[r['rungroup'],r['entry'],optics[r['rungroup']],r['foil'],r['ndel'],r['zfoil'],r['delta_low'],r['delta_high'],r['xscol'],r['yscol']]
            out.write('\t'.join(map(str,values))+'\n')


def verify_build(build):
    manifest=json.loads((build/'build.json').read_text())
    if manifest.get('schema')!='core_preallocated_v1' or manifest['sample']!='fit':
        raise ValueError('Explicit preallocated mode requires a training allocation/build manifest')
    for rel,h in manifest['outputs'].items():
        if digest(build/rel)!=h: raise ValueError(f'Build output changed: {rel}')
    allocated=json.loads((build/'allocation_manifest.json').read_text())
    if digest(build/'allocation_manifest.json')!=manifest['sample_manifest'] or allocated.get('mode')!='event_allocation' or allocated.get('status')!='quotas_ready':
        raise ValueError('Missing successful event-level allocation manifest')
    ids=read_tsv(build/'selected_ids.tsv')
    expected={(r['rungroup'],int(r['entry'])) for r in ids}
    if len(expected)!=len(ids) or len(ids)!=allocated['totals']['fit'] or not ids:
        raise ValueError('Invalid intended membership')
    return manifest,ids


def run(build,out,fit=False,max_events=200000,max_bytes=1073741824):
    manifest,ids=verify_build(build)
    if max_events<=0 or len(ids)>max_events or max_bytes<=0:
        raise ValueError('Resource limit exceeded; never truncate preallocated input')
    if out.exists(): raise FileExistsError(f'Immutable solver output exists: {out}')
    exe=shutil.which(os.environ.get('HCANA','hcana'))
    if not exe: raise ValueError('Load ROOT/HCANA or set HCANA=root')
    out.mkdir(parents=True)
    # Keep failure audits, but never certify a failed handoff.
    table,=build.glob('rungroups_*_inputs.tsv')
    args=[json.dumps('preallocated'),'-1',str(max_events),'-1',json.dumps(str(build/'root')),json.dumps(str(out)),
          json.dumps(str(build/'oldfit.dat')),json.dumps(str(table)),json.dumps(str(build/'optics.dat')),
          json.dumps(str(build/'solver_input.tsv')),'false' if fit else 'true',str(max_bytes)]
    with (out/'solver.log').open('w') as log:
        status=subprocess.run([exe,'-b','-l','-q',str(PROJECT/'fit_opt_matrix_gmm.C')+'('+','.join(args)+')'],cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT)
    usedpath=out/'solver_used_ids.tsv'
    used=read_tsv(usedpath) if usedpath.exists() else []
    expected={(r['rungroup'],int(r['entry'])) for r in ids}
    actual={(r['rungroup'],int(r['entry'])) for r in used}
    rejects=read_tsv(out/'solver_rejected_ids.tsv') if (out/'solver_rejected_ids.tsv').exists() else []
    if status.returncode or rejects or len(actual)!=len(used) or actual!=expected:
        write_tsv(out/'handoff_missing_ids.tsv',[dict(rungroup=r,entry=e) for r,e in sorted(expected-actual)])
        raise RuntimeError(f'Failed preallocated handoff; inspect {out}/solver.log and membership audits')
    (out/'acceptance.json').write_text(json.dumps(dict(status='verified',events=len(actual),fit_run=fit,
        build_sha256=digest(build/'build.json'),solver_sha256=digest(PROJECT/'fit_opt_matrix_gmm.C'),
        admission_sha256=digest(PROJECT/'preallocated_sample.h'),used_ids_sha256=digest(usedpath)),indent=2)+'\n')
    print(f'Exact membership verified: {len(actual)} events. Production fit run: {fit}')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path);p.add_argument('tag');p.add_argument('output',type=Path)
    p.add_argument('--fit',action='store_true',help='Explicitly solve after acceptance; omit for handoff review')
    p.add_argument('--max-events',type=int,default=200000);p.add_argument('--max-design-bytes',type=int,default=1073741824)
    a=p.parse_args()
    run((a.campaign/'06c_core_ntuple'/a.tag/'fit').resolve(),a.output.resolve(),a.fit,a.max_events,a.max_design_bytes)
