"""Direct-X equivalence, deterministic add/remove search, and exact campaign handoff."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from beam_search import Folds,DEFAULTS,search
from elastic_net import scaling,transform,macro_mse,SCALES
from core_sample import read_tsv
from test_elastic import fixture
import fit_elastic
import fit_beam
from elastic_convergence import run as pooled_run


class BeamTests(unittest.TestCase):
    def test_qr_matches_direct_x_and_pooled_cell_score(self):
        rng = np.random.default_rng(93)
        x = rng.normal(size=(240,7));x[:,1] = x[:,0]+1e-6*x[:,1]
        y = x[:,:3] @ rng.normal(size=(3,3)) + rng.normal(size=(240,3))*.2
        folds = np.arange(240)%3
        cells = np.array(['a']*130+['b']*70+['sparse']*40)
        prepared = Folds(x,y,folds,cells,1e-12,2)
        for mask in ((),(0,2),(0,1,2,3),tuple(range(7))):
            predicted,stats = prepared.predict(0,mask)
            direct = np.empty(len(x))
            for f in range(3):
                tr,va = folds!=f,folds==f
                state = scaling(x[tr],y[tr]);z,t = transform(x[tr],y[tr],state)
                zv,_ = transform(x[va],y[va],state)
                b,_,rank,s = np.linalg.lstsq(z[:,mask],t[:,0],rcond=1e-12)
                direct[va] = zv[:,mask] @ b*state[3][0]+state[2][0]
                self.assertEqual(stats[f]['rank'],rank)
                if mask: self.assertAlmostEqual(stats[f]['condition']/(s[0]/s[-1]),1.,places=7)
            np.testing.assert_allclose(predicted,direct,rtol=2e-7,atol=2e-7)
            score = prepared.evaluate(0,mask)['score']
            self.assertAlmostEqual(score/float(macro_mse((direct-y[:,0])*SCALES[0],cells)),1.,places=7)

    def test_add_remove_and_thread_determinism(self):
        def evaluate(mask):
            # Replacing term 0 with term 1 requires both addition and removal.
            score = 1.+(0 in mask)+2*(1 not in mask)+sum(i>=2 for i in mask)
            return dict(mask=mask,score=score,terms=len(mask)+1,condition=1.,rank_min=len(mask),valid=True)
        answers = []
        for threads in (1,8):
            log = []
            cfg = dict(DEFAULTS,beam=2,threads=threads,steps=4,patience=3,slack=0)
            result = search(evaluate,[(0,)],4,cfg,log.append)
            self.assertEqual(result['chosen']['mask'],(1,))
            self.assertEqual(len(log),len({r['mask'] for r in log}))
            answers.append((result['chosen'],[(r['mask'],r['score']) for r in log]))
        self.assertEqual(answers[0],answers[1])

    def test_seed_and_output_round_trip_without_holdout(self):
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            p = fixture(Path(tmp))
            with patch('elastic_diagnostics.training_plots'): fit_elastic.fit(p,'balanced','enet')
            (p/'config/elastic_refit.json').write_text(json.dumps(dict(alphas=[.003],l1=[.5],budgets=[10000])))
            pooled_run(p,'balanced','pooled',refit=True,pooled=True)
            (p/'06c_core_ntuple/balanced/holdout/build.json').unlink()
            selectors = {t:dict(alpha=.003,l1=.5) for t in fit_beam.TARGETS}
            cfg = dict(DEFAULTS,beam=2,steps=2,threads=2)
            fit_beam.run(p,'balanced','beam','pooled',cfg,selectors,check=True)
            self.assertFalse((p/'06e_beam_search/beam').exists())
            fit_beam.run(p,'balanced','beam','pooled',cfg,selectors)
            out = p/'06e_beam_search/beam'
            m = json.loads((out/'manifest.json').read_text())
            self.assertEqual(m['status'],'complete');self.assertEqual(m['training'],480)
            self.assertTrue((out/'plots/search.png').exists())
            self.assertTrue((out/'plots/foil_delta.png').exists())
            self.assertEqual(len(list((out/'seeds').glob('fold*.dat'))),3)
            self.assertEqual(set(np.load(out/'residuals.npz')['entry']),set(fit_elastic.load_sample(p,'balanced','fit')[0]['entry']))
            rows = fit_elastic.matrix_rows(out/'seed.dat')
            for model in ('beam','start','svd'):
                got = dict((t,c) for c,t in fit_elastic.matrix_rows(out/f'matrices/{model}.dat'))
                for c,t in rows:
                    self.assertEqual(c[3],got[t][3])
                    if t[4]: np.testing.assert_array_equal(c,got[t])
            for r in read_tsv(out/'tsv/terms.tsv'):
                if r['term']=='00000':self.assertEqual(r['beam'],'True')
            # Frozen input products are checked before a new search.
            (p/'06d_elastic_net/pooled/terms.tsv').write_text('changed')
            with self.assertRaisesRegex(ValueError,'Frozen input changed'):
                fit_beam.run(p,'balanced','second','pooled',cfg,selectors,check=True)

    def test_smallest_within_slack_and_rank_rejection(self):
        def evaluate(mask):
            scores = {():.9,(0,):1.004,(1,):2.,(0,1):1.}
            return dict(mask=mask,score=scores[mask],terms=len(mask)+1,condition=1.,
                        rank_min=len(mask),valid=bool(mask))
        cfg = dict(DEFAULTS,threads=1,steps=1,slack=.005)
        result = search(evaluate,[(0,),(0,1)],2,cfg,lambda r:None)
        self.assertEqual(result['best']['mask'],(0,1))
        self.assertEqual(result['chosen']['mask'],(0,))


if __name__=='__main__':unittest.main()
