"""Frozen ROOT round trip and HMS TFit export fixture; no production data fit."""
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import uproot
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from test_core_sample import fixture
import core_sample as cs
from reallocate_core import run,load_inputs
from build_core_fit import build
from preallocated_svd import verify_build

class FrozenTests(unittest.TestCase):
    def test_frozen_roundtrip_and_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'HMS_fixture';fixture(p)
            with patch('core_sample.plot_slice'):
                cs.run(p,'frozen')
            parent=p/'05c_core_sample/frozen'
            (parent/'metadata').mkdir()
            shutil.copy2(cs.PROJECT/'DATfiles/list_of_optics_run.dat',parent/'metadata/optics.dat')
            m=json.loads((parent/'manifest.json').read_text())
            m['outputs']['metadata/optics.dat']=cs.digest(parent/'metadata/optics.dat')
            (parent/'manifest.json').write_text(json.dumps(m))
            (p/'config/sieve_mask.json').write_text(json.dumps(dict(blocked=[[2,3],[5,5]])))
            (p/'config/oldfit.dat').write_text('! fixture\n 0 0 0 0 10000\n ----\n')
            with patch('core_sample.select_slice',side_effect=AssertionError('Density rerun forbidden')):
                out=run(p,'frozen','balanced')
            self.assertTrue((out/'plots/foil_0_delta2_training_cloud.png').exists())
            self.assertGreater(json.loads((out/'diagnostic_provenance.json').read_text())['cloud_scale'][1],1)
            with patch('balance_diagnostics.diagnostics'):
                out2=run(p,'balanced','repeat')
                preview=run(p,'frozen','preview',True)
                with self.assertRaises(FileExistsError):run(p,'frozen','balanced')
            a=json.loads((out/'allocation.json').read_text());b=json.loads((preview/'allocation.json').read_text())
            self.assertEqual([r['quota'] for r in a['leaves']],[r['quota'] for r in b['leaves']])
            self.assertEqual((out/'tsv/selected_ids.tsv').read_bytes(),(out2/'tsv/selected_ids.tsv').read_bytes())
            name='CoreSample_rg01_test.root'
            with uproot.open(parent/'root'/name) as f:original=f['CoreSample'].arrays(library='np')
            with uproot.open(out/'root'/name) as f:balanced=f['CoreSample'].arrays(library='np')
            for k in original:
                if k!='sample':np.testing.assert_array_equal(original[k],balanced[k])
            np.testing.assert_array_equal(original['sample']==2,balanced['sample']==2)
            for path in (parent/'models').glob('*'):
                self.assertEqual(cs.digest(path),cs.digest(out/'models'/path.name))
            if shutil.which('root'):
                with patch.dict(os.environ,HCANA='root'):
                    build(p,'balanced','fit')
                manifest,ids=verify_build(p/'06c_core_ntuple/balanced/fit')
                self.assertEqual(len(ids),a['training_total'])
            # Zero-capacity report is retained, but exporter cannot accept it.
            (p/'config/core_balance.json').write_text(json.dumps(dict(balance=dict(foil='equal',hole=.0001,delta=1.5))))
            with patch('balance_diagnostics.diagnostics'):
                zero=run(p,'frozen','zero')
            self.assertEqual(json.loads((zero/'allocation.json').read_text())['training_total'],0)
            with self.assertRaises(ValueError):build(p,'zero','fit')
            # Metadata boundary disagreement is rejected, not rebinned.
            metadata=parent/'metadata/optics.dat'
            lines=metadata.read_text().splitlines()
            i=next(i for i,line in enumerate(lines) if line.startswith('666701,'))
            edges=[float(v) for v in lines[i+2].split(',') if v.strip()]
            edges[2]=-4.
            lines[i+2]=','.join(map(str,edges));text='\n'.join(lines)+'\n'
            metadata.write_text(text)
            m['outputs']['metadata/optics.dat']=cs.digest(metadata)
            (parent/'manifest.json').write_text(json.dumps(m))
            with self.assertRaises(ValueError):load_inputs(p,'frozen')

if __name__=='__main__':unittest.main()
