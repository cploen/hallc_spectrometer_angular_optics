#!/usr/bin/env python3
"""Standard fold reports and matched-event afterburner matrix comparisons."""
import argparse
import json
from pathlib import Path
import re
import shutil
import tempfile

import numpy as np
import uproot
from beam_holdout import prepare as protected_inputs
from compare_matrices import read_matrix
from comparison_plots import heatmaps, stats
from core_sample import PROJECT, digest, read_tsv, write_tsv
from elastic_diagnostics import groups, residuals
from fit_elastic import checked
from sieve_afterburner import reconstruct

TARGETS = ('xptar', 'ytar', 'yptar', 'ztar')
LABELS = ('x′tar (mrad)', 'ytar (cm)', 'y′tar (mrad)', 'ztar (cm)')
NAMES = dict(old='2024 replay matrix', svd='Full-basis SVD', start='EN-selected SVD basis', beam='Final beam matrix')


def references(saved, beam, matrices, hashes):
    checked(saved/'seed.dat', beam['outputs']['seed.dat'], hashes)
    models = dict(old=read_matrix(saved/'seed.dat'), **matrices)
    files = dict(old=saved/'seed.dat', **{m:saved/f'matrices/{m}.dat' for m in matrices})
    offsets = {m:[0., 0., 0.] for m in models}
    notes = {m:'Complete frozen matrix; zero external offsets' for m in models}
    notes['old'] = 'Frozen 2024 replay seed matrix; zero external offsets'
    return models, files, offsets, dict(NAMES), notes


def afterburner_inputs(campaign, tag, a):
    a = dict(a)
    a['ry'], a['rz'] = np.empty(len(a['entry'])), np.empty(len(a['entry']))
    source = campaign/'05c_core_sample'/tag
    table, = source.glob('rungroups_*_inputs.tsv')
    settings = {r['rungroup']:r for r in read_tsv(table)}
    for name in np.unique(a['rungroup']):
        chosen = np.flatnonzero(a['rungroup'] == name)
        path = campaign/'06c_core_ntuple'/tag/'holdout/root'/f'Optics_{settings[name]["optics_id"]}_-1_fit_tree_gmm.root'
        with uproot.open(path) as root:
            raw = root['TFit'].arrays(['entry','reactycalc','ztar'], library='np')
        order = np.argsort(raw['entry'])
        if not np.array_equal(raw['entry'][order], a['entry'][chosen]):
            raise ValueError(f'{name}: fixed-vertex input IDs disagree with protected export')
        a['ry'][chosen], a['rz'][chosen] = raw['reactycalc'][order], raw['ztar'][order]
    if not np.isfinite(a['ry']).all() or not np.isfinite(a['rz']).all():
        raise ValueError('Nonfinite fixed-vertex afterburner inputs')
    return a


def predictions(a, models, offsets):
    predicted, iteration_rows, closure = {}, [], []
    for model, matrix in models.items():
        values = np.empty((len(a['entry']), 3))
        for name in np.unique(a['rungroup']):
            idx = np.flatnonzero(a['rungroup'] == name)
            angle = float(np.rad2deg(a['angle'][idx[0]]))
            xy, iterations, active, target, xtar = reconstruct(
                {k:a[k][idx] for k in ('xfp','xpfp','yfp','ypfp','ry','rz')},
                matrix, angle, offsets[model], return_target=True)
            values[idx] = target
            iteration_rows.append(dict(model=model,rungroup=name,n=len(idx),
                                       iteration_limit=int(active.sum()),max_updates=int(iterations.max())))
            if model == 'old':
                difference = np.sqrt((xy[:,0]-a['xs'][idx])**2 + (xy[:,1]-a['ys'][idx])**2)
                closure.append(dict(rungroup=name,n=len(idx),rms_xy_cm=float(np.sqrt(np.mean(difference**2))),
                                    p99_xy_cm=float(np.percentile(difference,99)),max_xy_cm=float(difference.max())))
        if not np.isfinite(values).all():
            raise ValueError(f'Nonfinite afterburner predictions: {model}')
        predicted[model] = values
    return predicted, iteration_rows, closure


def fold_plots(out, pooled, saved):
    import matplotlib.pyplot as plt
    rows = read_tsv(pooled/'fold_scores.tsv')
    cases = sorted({(float(r['alpha']),float(r['l1'])) for r in rows})
    if not cases:
        raise ValueError('No completed reduced-basis fold results')
    folds = sorted({int(r['fold']) for r in rows})
    fig, axes = plt.subplots(4,3,figsize=(15,13),constrained_layout=True)
    keys = ('rms','p90','terms','condition')
    for j,target in enumerate(TARGETS[:3]):
        for fold in folds:
            for model,style,marker in [('refit_svd','-','o'),('full_svd',':','x')]:
                selected = [next((r for r in rows if r['target']==target and int(r['fold'])==fold
                            and (float(r['alpha']),float(r['l1']))==case and r['model']==model),None) for case in cases]
                for i,key in enumerate(keys):
                    values = [np.nan if r is None else np.sqrt(float(r['mse'])) if key=='rms' else float(r[key]) for r in selected]
                    axes[i,j].plot(range(len(cases)),values,linestyle=style,marker=marker,color=f'C{fold%10}',
                                   label=f'Fold {fold+1}: '+('reduced' if model=='refit_svd' else 'full SVD'))
        for i,key in enumerate(keys):
            unit = 'cm' if target=='ytar' else 'mrad'
            axes[i,j].set(title=target,ylabel={'rms':f'Validation RMS ({unit})','p90':f'P90 |residual| ({unit})',
                              'terms':'Terms including constant','condition':'Scaled-design condition'}[key],
                          xticks=range(len(cases)),xticklabels=[f'α={a:g}\nL1={l:g}' for a,l in cases])
            if key=='condition': axes[i,j].set_yscale('log')
        axes[0,j].legend(fontsize=7)
    fig.suptitle('Saved training folds: reduced-basis options versus full SVD; no fits rerun')
    fig.savefig(out/'plots/reduced_basis_folds.png',dpi=140);plt.close(fig)
    rows = read_tsv(saved/'tsv/fold_scores.tsv')
    fig,axes = plt.subplots(2,4,figsize=(16,7),constrained_layout=True)
    for j,target in enumerate(TARGETS):
        for model in ('svd','start','beam'):
            selected = sorted([r for r in rows if r['target']==target and r['model']==model],key=lambda r:int(r['fold']))
            for i,key in enumerate(('rms','bias')):
                axes[i,j].plot([int(r['fold'])+1 for r in selected],[float(r[key]) for r in selected],'.-',label=NAMES[model])
                axes[i,j].set(title=LABELS[j],xlabel='Validation fold',ylabel=key.upper(),xticks=[f+1 for f in folds])
                axes[i,j].legend(fontsize=7)
    fig.suptitle('Beam-search fold behavior: adaptive development diagnostics')
    fig.savefig(out/'plots/beam_folds.png',dpi=140);plt.close(fig)


def prepare(campaign, tag, source):
    saved, beam, a, matrices, pool, counts, hashes = protected_inputs(campaign,tag,source)
    pooled = campaign/'06d_elastic_net'/beam['source']
    pm = json.loads((pooled/'manifest.json').read_text())
    if (pm.get('schema'),pm.get('status'),pm.get('sample')) != ('elastic_pool_v1','complete',tag):
        raise ValueError('Missing completed pooled study for these beam results')
    hashes[str(pooled/'manifest.json')] = digest(pooled/'manifest.json')
    for file in ('fold_scores.tsv','pooled.tsv','terms.tsv'):
        checked(pooled/file,pm['outputs'][file],hashes)
    for file in ('tsv/fold_scores.tsv','tsv/conditioning.tsv'):
        checked(saved/file,beam['outputs'][file],hashes)
    models,files,offsets,names,notes = references(saved,beam,matrices,hashes)
    a = afterburner_inputs(campaign,tag,a)
    return saved,pooled,beam,a,models,files,offsets,names,notes,pool,counts,hashes


def run(campaign,tag,name,source,check=False):
    target = campaign/'07_diagnostics'/name
    if not check and target.exists(): raise FileExistsError(f'Output exists: {target}')
    saved,pooled,beam,a,models,files,offsets,names,notes,pool,counts,hashes = prepare(campaign,tag,source)
    print('Frozen comparison matrices: '+', '.join(names[m] for m in models),flush=True)
    print('Protected pools: '+str(counts),flush=True)
    if check:
        print('Saved folds, matrix inputs, offsets and protected IDs verified. No residuals evaluated.');return
    target.parent.mkdir(parents=True,exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.standard_diagnostics_',dir=target.parent))
    try:
        for folder in ('plots','tsv','matrices','code'): (stage/folder).mkdir()
        predicted,iterations,closure = predictions(a,models,offsets)
        truth = np.column_stack((a['xptarT'],a['ytarT']/100,a['yptarT']))
        rr = {m:residuals(a,truth,v,np.zeros_like(truth)) for m,v in predicted.items()}
        details=[]
        for group,mask in groups(a,pool):
            for m in models:
                for j,t in enumerate(TARGETS):
                    details.append(dict(group,model=m,target=t,qa_low=int(mask.sum())<10,**stats(rr[m][mask,j])))
        write_tsv(stage/'tsv/residuals.tsv',details)
        write_tsv(stage/'tsv/afterburner_iterations.tsv',iterations)
        write_tsv(stage/'tsv/replay_closure.tsv',closure)
        write_tsv(stage/'tsv/matrices.tsv',[dict(model=m,label=names[m],path=str(files[m]),sha256=digest(files[m]),
                      xptar_offset_mrad=offsets[m][0],ytar_offset_cm=offsets[m][1],yptar_offset_mrad=offsets[m][2],source=notes[m]) for m in models])
        write_tsv(stage/'tsv/ids.tsv',[dict(rungroup=r,entry=int(e),pool=p) for r,e,p in zip(a['rungroup'],a['entry'],pool)])
        pairs=[(m,'old') for m in models if m!='old']+[('beam','svd'),('beam','start')]
        cell_rows=[r for r in details if r['level']=='foil_delta']
        lookup={(r['pool'],r['model'],r['target'],r['zfoil'],r['ndel']):r for r in cell_rows}
        changes=[]
        for row in cell_rows:
            for model,reference in pairs:
                if row['model']!=model:continue
                base=lookup[(row['pool'],reference,row['target'],row['zfoil'],row['ndel'])]
                ratio=row['rms']/base['rms'] if base['rms']>0 else None
                changes.append(dict(pool=row['pool'],model=model,reference=reference,target=row['target'],
                              zfoil=row['zfoil'],ndel=row['ndel'],n=row['n'],rms=row['rms'],reference_rms=base['rms'],
                              rms_ratio=ratio,rms_percent_change=100*(ratio-1) if ratio is not None else ''))
        write_tsv(stage/'tsv/cell_comparisons.tsv',changes)
        for p in ('core','noncore','unsupported'):
            heatmaps(stage,p,details,models=tuple(models),names=names,pairs=pairs,absolute_format='.4g')
            heatmaps(stage,p,details,models=('beam',),names=names,pairs=(),prefix='final_beam_',absolute_format='.4g')
        fold_plots(stage,pooled,saved)
        for file in ('fold_scores.tsv','pooled.tsv','terms.tsv'):shutil.copy2(pooled/file,stage/'tsv'/f'pooled_{file}')
        for file in ('fold_scores.tsv','conditioning.tsv'):shutil.copy2(saved/'tsv'/file,stage/'tsv'/f'beam_{file}')
        for m,path in files.items():shutil.copy2(path,stage/f'matrices/{m}.dat')
        np.savez_compressed(stage/'residuals.npz',entry=a['entry'],rungroup=a['rungroup'],pool=pool,zfoil=a['ztarT'],
                            ndel=a['ndel'],xscol=a['xscol'],yscol=a['yscol'],targets=np.array(TARGETS),**rr)
        for file in ('matrix_diagnostic_suite.py','comparison_plots.py','beam_holdout.py','sieve_afterburner.py',
                     'elastic_diagnostics.py','fit_elastic.py','elastic_net.py','core_sample.py','mixed_sample.py',
                     'preallocated_svd.py','compare_matrices.py','spectrometer_config.py','spectrometer_profiles.def','run_standard_diagnostics.sh'):
            shutil.copy2(PROJECT/file,stage/'code'/file)
        (stage/'README.md').write_text(
            '# Standard matrix diagnostics\n\n'
            'reduced_basis_folds.png: individual saved folds for each elastic-net-selected SVD basis; RMS, P90, term counts and conditioning. '
            'beam_folds.png: full-SVD, seed and beam validation fold RMS and bias. These use saved training-validation results; no fits were rerun.\n\n'
            '*_foil_delta.png: absolute RMS for every matrix on identical reserved events, with shared target color scales. '
            '*_change.png: RMS percentage change relative to the initial matrix and modern references; negative/blue means smaller RMS. '
            'final_beam_*_foil_delta.png: absolute RMS for the final beam matrix alone, not a ratio. '
            'Angles are mrad; ytar and derived ztar are cm. N is identical across matrices within each cell. '
            'Blank = no events; asterisk = N<10. RMS includes bias; spread is separately tabulated.\n\n'
            'tsv/cell_comparisons.tsv records foil/delta RMS ratios and percentages. tsv/residuals.tsv includes bias, spread, RMS, '
            'median and P90 down to setting/hole. Core, supported noncore and unsupported pools are separate; blocked labels have counts only.\n\n'
            'Every matrix uses the same fixed-vertex iterative-xtar afterburner and unchanged event identities. Complete matrix constants '
            'and documented external offsets are retained; no coefficients or offsets are fitted. ztar is derived from the target identity; '
            'this is not a new replay. tsv/replay_closure.tsv checks the initial matrix against saved sieve coordinates. '
            'Afterburner iteration limits are recorded and no events are silently dropped.\n\n'
            'The new matrices exclude these holdout IDs from fitting. Historical matrix training membership is unknown unless separately '
            'documented; historical comparisons describe behavior on the same sample, not an independent test of every old fit.\n')
        for path,sha in hashes.items():checked(Path(path),sha,{})
        manifest=dict(schema='standard_matrix_diagnostics_v1',status='complete',campaign=campaign.name,sample=tag,source=source,
                      inputs=hashes,counts=counts,models=names,offsets=offsets,offset_sources=notes,
                      convention='Common fixed-vertex iterative-xtar afterburner; saved delta; no fitting')
        manifest['outputs']={str(p.relative_to(stage)):digest(p) for p in stage.rglob('*') if p.is_file()}
        (stage/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        stage.rename(target)
    except BaseException:
        shutil.rmtree(stage);raise
    print(f'Standard diagnostics complete: {target}',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path);p.add_argument('tag');p.add_argument('name');p.add_argument('--source',required=True)
    p.add_argument('--check',action='store_true')
    a=p.parse_args()
    if not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',v) for v in (a.tag,a.name,a.source)):p.error('Invalid tag/name/source')
    campaign=(a.campaign if a.campaign.is_absolute() else PROJECT/a.campaign).resolve()
    try:run(campaign,a.tag,a.name,a.source,a.check)
    except (OSError,ValueError,KeyError,RuntimeError) as exc:p.exit(1,f'ERROR: {exc}\n')

if __name__=='__main__':main()
