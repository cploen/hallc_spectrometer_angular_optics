"""Targeted elastic solver and exact-membership integration checks."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import uproot
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fit_elastic as fe
from elastic_net import DEFAULTS, folds_for, scaling, transform, raw_coefficients, tune, validate
from core_sample import digest, write_tsv
from spectrometer_config import from_campaign


def fixture(root):
    campaign = root/'HMS_fixture'
    source = campaign/'05c_core_sample/balanced'
    build = campaign/'06c_core_ntuple/balanced'
    for p in (source/'metadata', source/'root', source/'tsv', campaign/'config', build/'fit/root', build/'holdout/root'):
        p.mkdir(parents=True, exist_ok=True)
    table = 'rungroups_fixture_inputs.tsv'
    write_tsv(source/table, [dict(rungroup='rg01_test', optics_id=999)])
    (source/'metadata/optics.dat').write_text('999,test,12.5,5,1,3\n-8,-3,0,3,8\n-10,0,10\n')
    (source/'metadata/sieve_mask.json').write_text(json.dumps(dict(blocked=[[2,3]])))
    seed = build/'fit/oldfit.dat'
    seed.write_text('! fixture\n --------------------------------------------\n'
                    ' 0 0 0 .25 10000\n 0 0 0 .12 01000\n 0 0 0 0 00100\n'
                    ' 0 0 0 0 00010\n 0 0 0 0 20000\n .01 .02 .03 .04 00001\n --------------------------------------------\n')
    spec = from_campaign(campaign)
    rng = np.random.default_rng(43)
    a = []; m = []; ids = []
    entry = 2**54
    for zidx, z in enumerate([-8, -3, 0, 3, 8]):
        for nd in range(2):
            for x, y in [(2, 2), (3, 4), (5, 6), (6, 3), (2, 3)]:
                for code in (1, 2):
                    if code == 1 and (x, y) == (2, 3): continue
                    for i in range(12 if code == 1 else 9):
                        entry += 1
                        quality = 2 if code == 1 else i % 3
                        xtar = rng.normal(0, .1)
                        xs, ys = spec.xs(x), spec.ys(y)
                        ang = np.deg2rad(12.5); s, c = np.sin(ang), np.cos(ang)
                        ymis = .1*(.52-.012*12.5+.002*12.5**2)
                        xp = (xs-xtar)/spec.sieve_distance
                        yp = (ys+ymis-z*s)/(spec.sieve_distance-z*c)
                        yt = z*(s-yp*c)-ymis
                        r = dict(entry=entry, foil=zidx, ndel=nd, xscol=x, yscol=y, sample=code,
                                 core_keep=int(quality == 2), xtar=xtar, xptarT=xp, yptarT=yp, ytarT=yt,
                                 xsT=xs, ysT=ys, xs=xs+rng.normal(0,.2), ys=ys+rng.normal(0,.2), ztarT=float(z),
                                 xfp=100*xp*.2+rng.normal(0,.01), xpfp=xp*1.7+rng.normal(0,.0002),
                                 yfp=yt*.8+rng.normal(0,.01), ypfp=yp*2+rng.normal(0,.0001), delta=-5.+nd*10.)
                        a.append(r)
                        m.append(dict(entry=entry, foil=zidx, ndel=nd, xscol=x, yscol=y, sample=code,
                                      core_keep=r['core_keep'], quality=quality, xsieve=r['xs'], ysieve=r['ys'], zfoil=float(z)))
                        if code == 1: ids.append(dict(rungroup='rg01_test', entry=entry))
    mask = {k: np.array([r[k] for r in m]) for k in m[0]}
    with uproot.recreate(source/'root/CoreSample_rg01_test.root') as f: f.mktree('CoreSample', mask)
    write_tsv(source/'tsv/selected_ids.tsv', ids)
    manifest = dict(schema='core_balance_v1', mode='event_allocation', status='quotas_ready', totals=dict(fit=len(ids)))
    manifest['outputs'] = {str(p.relative_to(source)): digest(p) for p in source.rglob('*') if p.is_file()}
    (source/'manifest.json').write_text(json.dumps(manifest))
    for sample, code in [('fit', 1), ('holdout', 2)]:
        d = {k: np.array([r[k] for r in a if r['sample'] == code]) for k in a[0]}
        bp = build/sample
        with uproot.recreate(bp/'root/Optics_999_-1_fit_tree_gmm.root') as f: f.mktree('TFit', d)
        if sample == 'fit':
            (bp/'allocation_manifest.json').write_bytes((source/'manifest.json').read_bytes())
            (bp/'selected_ids.tsv').write_bytes((source/'tsv/selected_ids.tsv').read_bytes())
        bm = dict(schema='core_preallocated_v1', sample=sample, sample_manifest=digest(source/'manifest.json'))
        bm['outputs'] = {str(p.relative_to(bp)): digest(p) for p in bp.rglob('*') if p.is_file()}
        (bp/'build.json').write_text(json.dumps(bm))
    (campaign/'config/elastic.json').write_text(json.dumps(dict(alphas=[.03, .003, .0003], l1=[.2, .8], max_iter=10000)))
    return campaign


class ElasticTests(unittest.TestCase):
    def test_fold_order_and_scaling(self):
        ids = [f'rg:{2**54+i}' for i in range(40)]
        strata = [str(i//5) for i in range(40)]
        f = folds_for(ids, strata, 3, 667)
        order = np.random.default_rng(5).permutation(40)
        g = folds_for([ids[i] for i in order], [strata[i] for i in order], 3, 667)
        np.testing.assert_array_equal(f[order], g)
        x = np.array([[1e-12, 4], [2e-12, 4], [4e-12, 4.]])
        y = np.arange(9).reshape(3,3).astype(float)
        state = scaling(x,y)
        z,t = transform(x,y,state)
        beta = np.ones((2,3))
        beta[1] = 0
        np.testing.assert_allclose(np.c_[np.ones(3),x] @ raw_coefficients(beta,state), z @ beta * state[3]+state[2])
        self.assertTrue(state[4][1]); self.assertFalse(state[4][0])
        with self.assertRaises(ValueError): validate(dict(DEFAULTS, alphas=[0]))
        with self.assertRaises(ValueError): validate(dict(DEFAULTS, folds=1))

    def test_pipeline_and_protected_separation(self):
        with tempfile.TemporaryDirectory() as tmp:
            campaign = fixture(Path(tmp))
            train, hashes = fe.load_sample(campaign,'balanced','fit')
            holdout, _ = fe.load_sample(campaign,'balanced','holdout')
            self.assertEqual(len(train['entry']),480)
            self.assertEqual(len(holdout['entry']),450)
            self.assertGreater(int(train['entry'].min()), 2**54)
            out = fe.fit(campaign,'balanced','first')
            # Fitting never consumes protected export values, even if unavailable.
            hp = campaign/'06c_core_ntuple/balanced/holdout/root/Optics_999_-1_fit_tree_gmm.root'
            content = hp.read_bytes(); hp.unlink()
            with patch('elastic_diagnostics.training_plots'):
                repeat = fe.fit(campaign,'balanced','repeat')
            hp.write_bytes(content)
            with np.load(out/'model.npz') as x, np.load(repeat/'model.npz') as y:
                np.testing.assert_array_equal(x['enet'], y['enet'])
            with self.assertRaises(FileExistsError): fe.fit(campaign,'balanced','first')
            # Render the actual portfolio once on the synthetic campaign.
            fe.evaluate(campaign,'first')
            self.assertTrue((out/'evaluation/plots/slice_0_holes.png').is_file())
            self.assertTrue((out/'evaluation/plots/slice_1_noncore_clouds.png').is_file())
            with self.assertRaises(FileExistsError): fe.evaluate(campaign,'first')
            report = json.loads((out/'evaluation/manifest.json').read_text())
            self.assertEqual(sum(report['counts'].values()),450)
            self.assertEqual(report['counts']['blocked'],90)
            # Native coefficient export leaves delta/fixed terms untouched.
            old = dict((t,c) for c,t in fe.matrix_rows(out/'seed.dat'))
            for c,t in fe.matrix_rows(out/'matrices/enet.dat'):
                self.assertEqual(c[3], old.get(t,np.zeros(4))[3])
                if t[4]: np.testing.assert_array_equal(c,old[t])
            hp.write_bytes(content+b'changed')
            with self.assertRaisesRegex(ValueError,'Frozen input changed'):
                fe.load_sample(campaign,'balanced','holdout')

    def test_nonconvergence_not_selected(self):
        rng=np.random.default_rng(123)
        x=rng.normal(size=(120,8)); x[:,1]=x[:,0]+rng.normal(0,.001,120)
        y=rng.normal(size=(120,3))+x[:,:3]
        cfg=dict(DEFAULTS,alphas=[1e-12],l1=[.5],max_iter=1,tol=1e-14)
        with self.assertRaisesRegex(ValueError,'No converged candidate'):
            tune(x,y,np.arange(120)%3,np.array(['a']*120),cfg)


if __name__=='__main__': unittest.main()
