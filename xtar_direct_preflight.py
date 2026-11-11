#!/usr/bin/env python3
"""Read-only ifarm evidence collection for the separate xtar predictor.

No fitting, selection changes, transport writes, or replay installation. Works
with the current checkout: no new numerical modules are imported. Prints a
bounded JSON report suitable for pasting back into the chat.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import traceback

DEFAULT_CAMPAIGN = '/u/group/nps/cploen/ML_HMS_OPTICS_JOINT_SURFACE/HMS_LowMom_LowDelta_2018_Accepted'
BASELINE = '/u/home/cploen/nps_replay_coeff2018_zero_20261110/DATFILES/hms_recon_coeff2018_new_zero_offsets.dat'
BASELINE_HASH = 'f715879845f777f2da2aa62b7e66756907929bdf4b7b40809ac55370156f8801'
OLD_HASH = 'b2036a935db3f17189fcf489d1017da030a8bb58de7d5a091384ef3d3b60c92c'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def excerpt(path, words, context=4, limit=200):
    lines = path.read_text().splitlines()
    indices = set()
    for i, line in enumerate(lines):
        if any(word.lower() in line.lower() for word in words):
            indices.update(range(max(0, i-context), min(len(lines), i+context+1)))
    indices = sorted(indices)
    return dict(path=str(path), sha256=sha(path), matched_lines=len(indices),
                truncated=len(indices)>limit,
                lines=[f'{i+1}: {lines[i]}' for i in indices[:limit]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkout', type=Path, default=Path.cwd())
    parser.add_argument('--campaign', type=Path, default=Path(DEFAULT_CAMPAIGN))
    parser.add_argument('--hcana-source', type=Path, default=Path('/u/group/halla/apps/hcana/1.1/src/src'))
    args = parser.parse_args()
    c, checkout, src = args.campaign, args.checkout, args.hcana_source
    report = dict(schema='separate_xtar_preflight_v1', mode='read only; no fit', errors=[])
    def collect(name, operation):
        try: report[name] = operation()
        except Exception as exc:
            report['errors'].append(dict(stage=name, error=str(exc), type=type(exc).__name__))
    def git():
        return {name:subprocess.check_output(['git', '-C', str(checkout), *command], text=True).strip()
                for name,command in [('head',['rev-parse','HEAD']),('status',['status','--short'])]}
    collect('checkout', git)
    report['python'] = sys.version
    report['versions'] = {}
    for module in ('numpy','uproot','scipy','scikit-learn','threadpoolctl'):
        try: report['versions'][module] = importlib.metadata.version(module)
        except importlib.metadata.PackageNotFoundError: report['versions'][module] = 'not installed'
    report['matrices'] = []
    for path, expected in [(c/'config/oldfit.dat',OLD_HASH),
                           (c/'06c_core_ntuple/inner20/fit/oldfit.dat',OLD_HASH),
                           (Path(BASELINE),BASELINE_HASH)]:
        def matrix(path=path,expected=expected):
            actual = sha(path)
            header = []
            for line in path.read_text().splitlines():
                header.append(line)
                if line.strip().startswith('---'): break
            report['matrices'].append(dict(path=str(path), sha256=actual, expected_match=actual==expected,
                                            header=header[:50], header_truncated=len(header)>50))
            if actual != expected: raise ValueError('Matrix does not match recorded preflight/handoff hash')
        collect('matrix_check_'+str(path), matrix)
    for name, words in {
        'THcHallCSpectrometer.cxx': ['CalculateTargetQuantities','hut[','sum[','fThetaOffset','fPhiOffset','fDeltaOffset','fXOffset','fYOffset','fXRot','fYRot'],
        'THcExtTarCor.cxx': ['xtar','x_tg','vertex','pointing','CalculateTargetQuantities','GetTheta','GetPhi'],
        'THcReactionPoint.cxx': ['vertex','GetPosition','GetDirection','pointing','beam->'],
    }.items():
        collect(name, lambda name=name,words=words: excerpt(src/name,words))
    for filename in ('config/comparison.json','config/core_balance.json','05c_core_sample/inner20/config.json'):
        collect(filename, lambda filename=filename: dict(sha256=sha(c/filename),content=json.loads((c/filename).read_text())))
    collect('ntuple_exporter', lambda: excerpt(checkout/'make_fit_ntuple_from_gmm.C',
        ['reactycalc','reactxcalc','raster','SetBranchAddress','xtarT=','preallocated','GMM','gmm'], context=1, limit=100))
    collect('allocation_implementation',lambda:excerpt(checkout/'allocate_inner20.py',
        ['quality','core_pool','core_keep','inputs','manifest'],context=2,limit=90))
    sys.path.insert(0,str(checkout))
    def membership():
        import numpy as np
        from fit_elastic import load_sample
        from elastic_convergence import saved_folds
        from core_sample import read_tsv
        a, hashes = load_sample(c, 'inner20', 'fit', xtar_coordinate='xtarT')
        path = c/'06e_beam_search/beam30_inner20/tsv/folds.tsv'
        rows = read_tsv(path)
        folds = saved_folds(a, rows, len({int(r['fold']) for r in rows}))
        coverage = []
        for rg in np.unique(a['rungroup']):
            for z in np.unique(a['ztarT']):
                mask = (a['rungroup']==rg)&(a['ztarT']==z)
                if mask.any():
                    coverage.append(dict(run=str(rg),foil_cm=float(z),n=int(mask.sum()),
                        core=int(np.sum(mask&(a['quality']==2))),inner_shoulder=int(np.sum(mask&(a['quality']==1))),
                        angle_deg=float(np.rad2deg(a['angle'][mask][0])),
                        delta_min=float(a['delta'][mask].min()),delta_max=float(a['delta'][mask].max())))
        keys = sorted((str(r),int(e)) for r,e in zip(a['rungroup'],a['entry']))
        return dict(n=len(folds), retained_branches=sorted(a), input_files_checked=len(hashes),
            input_hashes=hashes, fold_sha256=sha(path),
            membership_sha256=hashlib.sha256(json.dumps(keys,separators=(',',':')).encode()).hexdigest(),
            fold_counts={str(f):int(np.sum(folds==f)) for f in np.unique(folds)},
            coverage=coverage,
            note='Original-oldfit historical provenance and full HCANA configuration require further evidence; hash equality alone is not provenance.')
    collect('frozen_training_validation',membership)
    print(json.dumps(report,indent=2,allow_nan=False))
    return int(bool(report['errors']))


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    raise SystemExit(main())
