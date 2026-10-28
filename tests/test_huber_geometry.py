"""Exercise the frozen-event Huber adapter independently of TFit export."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import uproot

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from core_sample import digest, write_tsv
from fit_huber import load_frozen
from spectrometer_config import HMS_XPTAR_DEFINITION


class HuberGeometryTests(unittest.TestCase):
    def test_frozen_vertex_geometry_and_reference_definition(self):
        with tempfile.TemporaryDirectory() as folder:
            campaign=Path(folder)/'HMS_fixture'
            source=campaign/'05c_core_sample/min10'
            (source/'root').mkdir(parents=True)
            table=source/'inputs.tsv';table.write_text('fixture')
            metadata=source/'optics.dat';metadata.write_text('fixture')
            maskpath=source/'mask.json';maskpath.write_text('{}')
            (source/'manifest.json').write_text('{}')
            a={k:np.zeros(3) for k in ('xfp','xpfp','yfp','ypfp','xtar','delta','xbpm_tar')}
            a.update(entry=np.arange(3,dtype=np.int64),run=np.ones(3,dtype=np.int32),
                     sample=np.ones(3,dtype=np.int32),quality=np.full(3,2,dtype=np.int32),
                     foil=np.arange(3,dtype=np.int32),ndel=np.zeros(3,dtype=np.int32),
                     xscol=np.array([0,4,8],dtype=np.int32),yscol=np.full(3,4,dtype=np.int32),
                     zfoil=np.array([-8.,0.,8.]),reacty=np.array([-.1,0.,.08]),
                     delta_low=np.full(3,-10.),delta_high=np.full(3,10.))
            tree=source/'root/CoreSample_rg.root'
            with uproot.recreate(tree) as f:f.mktree('CoreSample',a)
            inp=dict(source=source,table=table,metadata=metadata,maskpath=maskpath,
                     settings=[dict(rungroup='rg',optics_id='1')],blocked=set(),
                     meta=dict(rg=dict(foils=[-8.,0.,8.],edges=[-10.,10.],angle_deg=12.49)),
                     manifest=dict(outputs={str(tree.relative_to(source)):digest(tree)}))
            beam=campaign/'06e_beam_search/beam10';(beam/'tsv').mkdir(parents=True)
            for rel in ('seed.dat','matrices/svd.dat','matrices/start.dat','matrices/beam.dat'):
                path=beam/rel;path.parent.mkdir(exist_ok=True);path.write_text('fixture')
            ids=[dict(rungroup='rg',entry=i) for i in range(3)]
            write_tsv(beam/'tsv/folds.tsv',ids)
            selected=campaign/'05c_core_sample/equal15/tsv/selected_ids.tsv'
            selected.parent.mkdir(parents=True);write_tsv(selected,ids)
            manifest=dict(xptar_definition=HMS_XPTAR_DEFINITION,
                          outputs={str(p.relative_to(beam)):digest(p) for p in beam.rglob('*') if p.is_file()})
            (beam/'manifest.json').write_text(json.dumps(manifest))
            with patch('fit_huber.load_inputs',return_value=inp):
                result,_,_=load_frozen(campaign)
                xmis=.1*(2.37-.086*12.49+.0012*12.49**2)
                vertex=-a['reacty']-xmis
                intercept=vertex-result['xptarT']*a['zfoil']*np.cos(np.deg2rad(12.49))
                np.testing.assert_allclose(intercept+168*result['xptarT'],result['xsT'],atol=1e-13)
                del manifest['xptar_definition']
                (beam/'manifest.json').write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError,'old HMS xptarT definition'):
                    load_frozen(campaign)


if __name__=='__main__':
    unittest.main()
