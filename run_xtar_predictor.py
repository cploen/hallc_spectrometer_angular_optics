#!/usr/bin/env python3
"""Separate x-target proof of principle on the frozen 2018 allocation.

All transport coefficients stay fixed. Offline raw-coordinate polynomial
comparison, not historical replay reproduction. No HCANA installation.
"""
import argparse
import itertools
import importlib.metadata
import json
from pathlib import Path
import re
import shutil
import sys
import numpy as np
from threadpoolctl import threadpool_limits
from core_sample import digest, read_tsv, write_tsv
from fit_elastic import load_sample, design, matrix_rows, checked
from elastic_convergence import saved_folds
from xtar_predictor import fit_training_design, export_model, predict
from xtar_extended_heatmaps import metrics, write_report, MODELS, SCALES

REFERENCE_HASH = 'b2036a935db3f17189fcf489d1017da030a8bb58de7d5a091384ef3d3b60c92c'
COMPLETED_HASH = 'f715879845f777f2da2aa62b7e66756907929bdf4b7b40809ac55370156f8801'
BASELINE = Path('/u/home/cploen/nps_replay_coeff2018_zero_20261110/DATFILES/hms_recon_coeff2018_new_zero_offsets.dat')
LABELS = dict(zip(MODELS, ['oldfit reference (2018 header; historical provenance uncertain)',
    'Completed 2018 refit: zero-offset export', 'Completed refit + separate xtar predictor']))
ASSUMPTIONS = '''OFFLINE PROOF OF PRINCIPLE, NOT A NEW OR REPRODUCED HCANA REPLAY.
The same saved raw focal-plane measurements feed each coefficient polynomial.
Position inputs are converted cm -> m. No focus shifts, detector offsets,
angular offsets/slopes or external target offsets are applied in this declared
comparison convention. Existing matrix constants are evaluated exactly as stored.
Reference and completed-refit models use the same saved replay xtar input.
The candidate substitutes a separate measured-input xtar prediction once;
there is no reconstructed-vertex update or historical afterburner iteration.
All matrix coefficients remain fixed, including xtar-dependent coefficients
and the delta column. Selection, reserved IDs and delta binning are unchanged.
Geometric truth uses assigned foil/hole labels only for targets/diagnostics.
Previously inspected reserved events are development evidence, not a fresh test.
Historical original-2018 matrix provenance and historical replay closure remain
uncertified. Raw-polynomial changes can demonstrate sensitivity/improvement
under this convention; eventual NPS adoption requires matched HCANA replays.'''
ANCHOR = dict(expression='-reactycalc - 0.1*(2.37-0.086*d+0.0012*d*d)',
              d='min(abs(angle_deg),50)', unit='cm',
              source='saved H.react.y beam y; current geometric HMS pointing profile',
              limitation='profile is declared for this offline experiment, not certified historical replay pointing')


def measured_problem(a):
    """20 degree-1/2 monomials of measured inputs times cos(angle).

    Variables: xfp_m,xpfp_rad,yfp_m,ypfp_rad,beam_y_m. No angle-only or
    foil/hole/ztar/xtar/transport-output features. The cos factor follows the
    longitudinal projection. Predictor intercept is free and separate.
    """
    degree = np.minimum(np.abs(np.rad2deg(a['angle'])), 50.)
    pointing = .1*(2.37-.086*degree+.0012*degree**2)
    anchor = -np.asarray(a['reactycalc'])-pointing
    names = ['xfp_m','xpfp_rad','yfp_m','ypfp_rad','beam_y_m']
    units = ['m','rad','m','rad','m']
    sources = ['H.dc.x_fp/100','H.dc.xp_fp','H.dc.y_fp/100','H.dc.yp_fp','H.react.y/100']
    u = np.column_stack([a['xfp']/100.,a['xpfp'],a['yfp']/100.,a['ypfp'],a['reactycalc']/100.])
    columns, features = [], []
    for degree in (1,2):
        for idx in itertools.combinations_with_replacement(range(len(names)),degree):
            columns.append(np.prod(u[:,idx],axis=1)*np.cos(a['angle']))
            features.append(dict(name='cos_angle__'+'__'.join(names[i] for i in idx),
                expression='cos(angle_rad)*'+'*'.join(names[i] for i in idx),
                unit='*'.join(units[i] for i in idx),
                source='; '.join(sources[i] for i in idx)))
    x = np.column_stack(columns)
    if not np.isfinite(x).all() or not np.isfinite(anchor).all():
        raise ValueError('Nonfinite measured feature/anchor; no event dropping')
    return x,anchor,features


def polynomial_targets(a, rows, xtar_cm):
    # Explicit raw-coordinate convention. Not a substitute for HCANA's optional
    # focus/offset/slope transformations when they are configured in a replay.
    coordinates = dict(a,xtar=np.asarray(xtar_cm))
    values = design(coordinates,[t for c,t in rows]) @ np.array([c for c,t in rows])
    return np.column_stack((values[:,0],values[:,2],100*values[:,1],xtar_cm))


def cells(a):
    return np.array([f'{z}:{d}' for z,d in zip(a['ztarT'],a['ndel'])])


def support_keys(a):
    return list(zip(map(float,a['ztarT']),map(int,a['ndel']),map(int,a['xscol']),map(int,a['yscol'])))


def run(args):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',args.name):
        raise ValueError('Invalid output name')
    reference = args.campaign/'config/oldfit.dat'
    hashes = {}
    checked(reference,args.reference_sha256,hashes)
    checked(args.baseline,args.baseline_sha256,hashes)
    a, inputs = load_sample(args.campaign,args.tag,'fit',xtar_coordinate='xtarT')
    hashes.update(inputs)
    folds_path = args.campaign/'06e_beam_search/beam30_inner20/tsv/folds.tsv'
    if args.folds_from is not None: folds_path = args.folds_from
    hashes[str(folds_path)] = digest(folds_path)
    table = read_tsv(folds_path)
    folds = saved_folds(a,table,len({int(r['fold']) for r in table}))
    angle_folds=[]
    for f in sorted(set(folds)):
        if set(a['angle'][folds==f]) != set(a['angle']):
            raise ValueError('Each saved training fold must represent every accepted angle')
        for angle in sorted(set(a['angle'])):
            angle_folds.append(dict(fold=int(f),angle_deg=float(np.rad2deg(angle)),
                                    n=int(np.sum((folds==f)&(a['angle']==angle)))))
    x,anchor,features = measured_problem(a)
    print(f'{len(x):,} verified training rows; {len(features)} measured feature columns; beam width 30',flush=True)
    print('Separate cm-valued xtar predictor. Frozen transport. Declared offline raw-polynomial comparison.',flush=True)
    if not args.fit:
        print('Check only: no fit, reserved evaluation or output writes.',flush=True)
        return
    out = args.campaign/'06h_xtar_predictor'/args.name
    out.mkdir(parents=True,exist_ok=False)
    manifest = dict(schema='separate_xtar_experiment_v1',status='running',tag=args.tag,
        training=len(x),inputs=hashes,assumptions=ASSUMPTIONS,model_labels=LABELS,
        basis='20 linear/quadratic measured monomials times cos(angle), separate free intercept',
        versions={m:importlib.metadata.version(m) for m in ('numpy','scipy','scikit-learn','uproot','threadpoolctl')},
        anchor=ANCHOR,holdout_status='not opened',historical_closure='not required; unresolved')
    def save(): (out/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    save()
    try:
        for name in ['run_xtar_predictor.py','xtar_predictor.py','xtar_extended_heatmaps.py',
                     'elastic_net.py','beam_search.py','fit_elastic.py','elastic_refit.py',
                     'elastic_convergence.py','core_sample.py','preallocated_svd.py','mixed_sample.py',
                     'spectrometer_config.py','spectrometer_profiles.def']:
            source=Path(__file__).with_name(name)
            (out/'code').mkdir(exist_ok=True)
            shutil.copyfile(source,out/'code'/name)
        shutil.copyfile(reference,out/'oldfit_reference.dat')
        shutil.copyfile(args.baseline,out/'completed_refit_fixed.dat')
        write_tsv(out/'angle_folds.tsv',angle_folds)
        write_tsv(out/'folds.tsv',[dict(rungroup=r,entry=int(e),fold=int(f)) for r,e,f in zip(a['rungroup'],a['entry'],folds)])
        en_config={}
        policy=args.campaign/'config/elastic.json'
        if policy.exists():
            en_config=json.loads(policy.read_text());hashes[str(policy)]=digest(policy)
        beam_config=dict(threads=args.threads,steps=args.steps)
        with threadpool_limits(limits=args.threads):
            fitted=fit_training_design(x,a['xtarT'],folds,cells(a),en_config=en_config,
                                      beam_config=beam_config,known_offset_cm=anchor)
        validity=dict(training_measured_min=x.min(axis=0).tolist(),training_measured_max=x.max(axis=0).tolist(),
                      training_angle_deg=sorted(set(map(float,np.rad2deg(a['angle'])))),
                      training_delta_percent=[float(a['delta'].min()),float(a['delta'].max())],
                      extrapolation='No automatic clipping; out-of-range inputs must be flagged by a replay adapter')
        model=export_model(out/'xtar_predictor.json',fitted,features=features,provenance=hashes,
                           validity=validity,anchor=ANCHOR)
        prediction=predict(model,x,anchor_cm=anchor)
        native=anchor+np.c_[np.ones(len(x)),x]@fitted['coefficients_cm']
        np.testing.assert_allclose(prediction,native,rtol=1e-12,atol=1e-11)
        write_tsv(out/'enet_cv.tsv',fitted['enet_cv']);write_tsv(out/'enet_summary.tsv',fitted['enet_summary'])
        write_tsv(out/'conditioning.tsv',fitted['fold_fits'])
        (out/'search.json').write_text(json.dumps(fitted['search'],indent=2)+'\n')
        np.savez_compressed(out/'training_development.npz',rungroup=a['rungroup'],entry=a['entry'],
            fold=folds,truth_xtar_cm=a['xtarT'],saved_xtar_cm=a['xtar'],oof_xtar_cm=fitted['oof_prediction_cm'])
        manifest.update(selected_columns=fitted['selected_columns'],final_fit=fitted['final_fit'],
                        enet_choices=fitted['enet_choices'],full_design_svd=fitted['full_design_svd'])
        if args.evaluate:
            h,hold_hashes=load_sample(args.campaign,args.tag,'holdout',xtar_coordinate='xtarT')
            train_ids=set(zip(a['rungroup'],map(int,a['entry'])))
            hold_ids=set(zip(h['rungroup'],map(int,h['entry'])))
            if train_ids&hold_ids: raise ValueError('Training/reserved ID overlap')
            hx,ha,hfeatures=measured_problem(h)
            if hfeatures!=features: raise ValueError('Evaluation feature definitions changed')
            candidate_x=predict(model,hx,anchor_cm=ha)
            matrices=[matrix_rows(out/'oldfit_reference.dat'),matrix_rows(out/'completed_refit_fixed.dat')]
            bundles={m:dict(rungroup=h['rungroup'],entry=h['entry'],values=p) for m,p in zip(MODELS,[
                polynomial_targets(h,matrices[0],h['xtar']),polynomial_targets(h,matrices[1],h['xtar']),
                polynomial_targets(h,matrices[1],candidate_x)])}
            truth=np.column_stack([h['xptarT'],h['yptarT'],h['ytarT'],h['xtarT']])
            training_support=set(support_keys(a))
            supported=np.array([k in training_support for k in support_keys(h)])
            rows,counts=metrics(h,truth,bundles,physical_foils=np.unique(a['ztarT']),supported=supported)
            write_report(out/'reserved_comparison',rows,counts,labels=LABELS,assumptions=ASSUMPTIONS)
            summary=[]
            for population in ('core','supported_noncore'):
                for target in ('xptar','yptar','ytar','xtar'):
                    for model_name in MODELS:
                        group=[r for r in rows if r['population']==population and r['target']==target and r['model']==model_name and r['n']]
                        n=sum(r['n'] for r in group)
                        if not n: continue
                        summary.append(dict(population=population,target=target,model=model_name,n=n,unit=group[0]['unit'],
                            bias=sum(r['n']*r['bias'] for r in group)/n,
                            rms_including_bias=float(np.sqrt(sum(r['n']*r['rms_including_bias']**2 for r in group)/n))))
            lookup={(r['population'],r['target'],r['model']):r for r in summary}
            for row in summary:
                for reference_name,suffix in [(MODELS[0],'oldfit'),(MODELS[1],'refit')]:
                    reference_row=lookup[row['population'],row['target'],reference_name]
                    row['rms_change_vs_'+suffix]=row['rms_including_bias']-reference_row['rms_including_bias']
            if summary: write_tsv(out/'reserved_summary.tsv',summary)
            np.savez_compressed(out/'reserved_predictions.npz',rungroup=h['rungroup'],entry=h['entry'],
                delta=h['delta'],ztarT=h['ztarT'],quality=h['quality'],excluded=h['excluded'],
                truth=truth,**{m:b['values'] for m,b in bundles.items()})
            manifest.update(holdout_status='evaluated once after selection; reused development holdout',
                reserved_counts=counts,reserved_inputs=hold_hashes,
                reserved_feature_outside_training_box=int(np.sum(np.any((hx<x.min(axis=0))|(hx>x.max(axis=0)),axis=1))))
            hashes.update(hold_hashes)
        for p,expected in hashes.items(): checked(Path(p),expected,{})
        manifest.update(status='complete',outputs={str(p.relative_to(out)):digest(p)
                        for p in out.rglob('*') if p.is_file() and p.name!='manifest.json'})
        save()
    except BaseException as exc:
        manifest.update(status='incomplete',error=str(exc));save();raise
    print('Complete:',out,flush=True)
    print('No matrix/predictor installed in replay. No GMM fitting performed.',flush=True)
    if args.evaluate:
        print('RESERVED POPULATION COUNTS:',json.dumps(manifest['reserved_counts']),flush=True)
        summary_path=out/'reserved_summary.tsv'
        if summary_path.exists(): print(summary_path.read_text(),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path);p.add_argument('--tag',default='inner20')
    p.add_argument('--name',default='direct_quadratic_beam30_v1')
    p.add_argument('--baseline',type=Path,default=BASELINE)
    p.add_argument('--reference-sha256',default=REFERENCE_HASH)
    p.add_argument('--baseline-sha256',default=COMPLETED_HASH)
    p.add_argument('--folds-from',type=Path)
    p.add_argument('--threads',type=int,default=4);p.add_argument('--steps',type=int,default=20)
    p.add_argument('--fit',action='store_true');p.add_argument('--evaluate',action='store_true')
    args=p.parse_args()
    if args.evaluate and not args.fit: p.error('--evaluate requires --fit')
    if args.threads<1 or args.steps<1:p.error('threads/steps must be positive')
    sys.dont_write_bytecode=True
    run(args)


if __name__=='__main__':main()
