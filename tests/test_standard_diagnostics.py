"""Standard fold displays and common-afterburner, cell-by-cell comparisons."""
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import uproot
from test_beam_holdout import fixture as protected_fixture, write_root
from core_sample import digest, read_tsv, write_tsv
from matrix_diagnostic_suite import afterburner_inputs, prepare, predictions, references, run


def fixture(directory):
    campaign,saved,beam=protected_fixture(directory)
    holdout=campaign/'06c_core_ntuple/balanced/holdout'
    path=holdout/'root/Optics_999_-1_fit_tree_gmm.root'
    with uproot.open(path) as f:a=f['TFit'].arrays(library='np')
    a['ztar']=a['ztarT'].copy()
    for field in ('xptar','ytar','yptar'):a[field]=np.zeros(len(a['entry']))
    write_root(path,'TFit',a)
    hm=json.loads((holdout/'build.json').read_text())
    hm['outputs']['root/Optics_999_-1_fit_tree_gmm.root']=digest(path)
    (holdout/'build.json').write_text(json.dumps(hm))
    shutil.copy2(campaign/'06c_core_ntuple/balanced/fit/oldfit.dat',saved/'seed.dat')
    pooled=campaign/'06d_elastic_net/pooled_balanced';pooled.mkdir(parents=True)
    scores=[]
    for target in ('xptar','ytar','yptar'):
        for fold in range(3):
            for alpha in (.0001,.00003):
                for model in ('full_svd','refit_svd'):
                    scores.append(dict(target=target,fold=fold,alpha=alpha,l1=.5,model=model,
                                       mse=.04+.01*fold,cell_mse=.05,p90=.3,median=.1,n=100,
                                       terms=6 if model=='full_svd' else 4,condition=10+fold))
    write_tsv(pooled/'fold_scores.tsv',scores)
    write_tsv(pooled/'pooled.tsv',[dict(target='ytar',model='refit_svd',mse=.04)])
    write_tsv(pooled/'terms.tsv',[dict(target='ytar',fold=0,term='10000',selected=True)])
    pm=dict(schema='elastic_pool_v1',status='complete',sample='balanced',
            outputs={p.name:digest(p) for p in pooled.iterdir()})
    (pooled/'manifest.json').write_text(json.dumps(pm))
    fold_scores=[];conditioning=[]
    for target in ('xptar','ytar','yptar','ztar'):
        for model in ('svd','start','beam'):
            for fold in range(3):
                fold_scores.append(dict(target=target,model=model,fold=fold,rms=.2+.01*fold,bias=.005))
                if target!='ztar':conditioning.append(dict(target=target,model=model,fold=fold,rank=3,condition=10))
    write_tsv(saved/'tsv/fold_scores.tsv',fold_scores)
    write_tsv(saved/'tsv/conditioning.tsv',conditioning)
    beam.update(source='pooled_balanced')
    beam['inputs'][str(pooled/'manifest.json')]=digest(pooled/'manifest.json')
    beam['outputs']={str(p.relative_to(saved)):digest(p) for p in saved.rglob('*') if p.is_file() and p.name!='manifest.json'}
    (saved/'manifest.json').write_text(json.dumps(beam))
    return campaign,saved


class StandardTests(unittest.TestCase):
    def test_all_reports_without_refitting(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign,saved=fixture(Path(temporary))
            before={p:digest(p) for p in campaign.rglob('*') if p.is_file()}
            run(campaign,'balanced','standard','beam20',check=True)
            self.assertFalse((campaign/'07_diagnostics').exists())
            with patch('numpy.linalg.lstsq',side_effect=AssertionError('No refitting')):
                run(campaign,'balanced','standard','beam20')
            out=campaign/'07_diagnostics/standard'
            manifest=json.loads((out/'manifest.json').read_text())
            self.assertEqual(set(manifest['models']),{'old','svd','start','beam'})
            for rel,sha in manifest['outputs'].items():self.assertEqual(digest(out/rel),sha)
            self.assertEqual(len(list((out/'plots').glob('*.png'))),11)
            rows=read_tsv(out/'tsv/residuals.tsv')
            self.assertNotIn('blocked',{r['pool'] for r in rows})
            for pool in ('core','noncore','unsupported'):
                self.assertTrue((out/f'plots/final_beam_{pool}_foil_delta.png').exists())
            comparisons=read_tsv(out/'tsv/cell_comparisons.tsv')
            self.assertIn(('beam','old'),{(r['model'],r['reference']) for r in comparisons})
            for row in comparisons:
                self.assertAlmostEqual(float(row['rms_ratio']),float(row['rms'])/float(row['reference_rms']))
                self.assertAlmostEqual(float(row['rms_percent_change']),100*(float(row['rms_ratio'])-1))
            with np.load(out/'residuals.npz') as residuals:
                self.assertGreater(int(residuals['entry'].min()),2**54)
                for row in rows:
                    if row['level']!='foil_delta':continue
                    mask=(residuals['pool']==row['pool'])&(residuals['zfoil']==float(row['zfoil']))&(residuals['ndel']==int(row['ndel']))
                    j=list(residuals['targets']).index(row['target'])
                    self.assertAlmostEqual(float(row['rms']),float(np.sqrt(np.mean(residuals[row['model']][mask,j]**2))))
            self.assertEqual(before,{p:digest(p) for p in before})
            with self.assertRaises(FileExistsError):run(campaign,'balanced','standard','beam20')
            if os.environ.get('STANDARD_PREVIEW'):
                preview=Path(os.environ['STANDARD_PREVIEW']);preview.mkdir(parents=True,exist_ok=True)
                for name in ('reduced_basis_folds.png','core_change.png','final_beam_core_foil_delta.png'):
                    shutil.copy2(out/'plots'/name,preview/name)

    def test_2024_reference_matches_legacy_afterburner(self):
        from legacy_afterburner import export_tree
        with tempfile.TemporaryDirectory() as temporary:
            campaign,saved=fixture(Path(temporary))
            _,_,_,a,models,files,offsets,names,notes,_,_,_=prepare(campaign,'balanced','beam20')
            self.assertEqual(names['old'], '2024 replay matrix')
            predicted,iterations,closure=predictions(a,models,offsets)
            self.assertEqual(len(iterations),4)
            self.assertEqual(len(closure),1)
            raw=campaign/'06c_core_ntuple/balanced/holdout/root/Optics_999_-1_fit_tree_gmm.root'
            for model in ('old','beam'):
                target=Path(temporary)/f'{model}.root'
                export_tree(raw,target,models[model],12.5,[0,0,0])
                with uproot.open(target) as root: expected=root['TFit'].arrays(library='np')
                np.testing.assert_array_equal(expected['entry'],a['entry'])
                for j,field in enumerate(('xptar','ytar','yptar')):
                    np.testing.assert_allclose(predicted[model][:,j],expected[field]/(100 if field=='ytar' else 1),rtol=1e-13,atol=1e-13)
            with (saved/'seed.dat').open('a') as stream:stream.write('! changed\n')
            with self.assertRaisesRegex(ValueError,'Frozen input changed'):prepare(campaign,'balanced','beam20')

if __name__=='__main__':unittest.main()
