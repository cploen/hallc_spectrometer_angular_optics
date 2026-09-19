"""ROOT integration check: known slopes, core/foil isolation and arm selection."""
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

PROJECT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('root'), 'ROOT is required')
class GeometryStudyTest(unittest.TestCase):
    def test_reference_rays(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            fixture = base/'fixture.C'
            fixture.write_text(r'''
#include <TFile.h>
#include <TTree.h>
#include <cmath>
void fixture(const char* path, bool shms) {
 TFile file(path,"RECREATE");TTree tree("CoreSample","");
 double entry=0,run=1,core_keep=1,zfoil=0,xscol=shms?5:4,yscol=xscol,ndel=0;
 double delta=0,delta_low=-10,delta_high=10,xfp=0,xpfp=0,yfp=0,ypfp=0;
 double reactx=0,reacty=0,reactz=0,xptar=0,sumnpe=5,etracknorm=1,core_score=0;
#define B(v) tree.Branch(#v,&v);
 B(entry) B(run) B(core_keep) B(zfoil) B(xscol) B(yscol) B(ndel)
 B(delta) B(delta_low) B(delta_high) B(xfp) B(xpfp) B(yfp) B(ypfp)
 B(reactx) B(reacty) B(reactz) B(xptar) B(sumnpe) B(etracknorm) B(core_score)
 double angle=12.49*acos(-1.)/180., L=shms?253.:168.;
 double xmis=shms?-.126:.1*(2.37-.086*12.49+.0012*12.49*12.49);
 for(double foil:{-8.,0.,8.}) for(int i=-100;i<=100;++i) {
   zfoil=reactz=foil;reacty=.0005*i;core_score=1.-std::abs(i)/200.;
   xfp=sin(i*.13);xpfp=cos(i*.17);yfp=sin(i*.19);ypfp=cos(i*.23);
   xptar=(reacty+xmis)/(L-zfoil*cos(angle))+.0002*cos(i*.37);
   core_keep=1;tree.Fill();++entry;
   core_keep=0;xptar=99;tree.Fill();++entry;
 }
 tree.Write();
}
''')
            metadata=base/'metadata.dat'
            metadata.write_text('1,fixture,12.49,3,1,2\n-8,0,8\n-10,10\n')
            for arm in ('HMS','SHMS'):
                campaign=base/f'{arm}_fixture';(campaign/'config').mkdir(parents=True)
                core=campaign/'core.root'
                command=f'{fixture}({json.dumps(str(core))},{str(arm=="SHMS").lower()})'
                subprocess.run(['root','-l','-b','-q',command],check=True,capture_output=True,text=True)
                (campaign/'config/rungroups_fixture_inputs.tsv').write_text(
                    'rungroup\toptics_id\tangle_deg\trootfile\nrg01\t1\t12.49\tmissing.root\n')
                config=dict(core_file='core.root',output='results',holes=[[4,4] if arm=='HMS' else [5,5]],min_events=30)
                (campaign/'config/geometry_study.json').write_text(json.dumps(config))
                for foil in (0,8):
                    run=subprocess.run([sys.executable,str(PROJECT/'geometry_study.py'),str(campaign),'--foil',str(foil),'--central-dense-half'],
                        env={**os.environ,'OPTICS_METADATA':str(metadata)},capture_output=True,text=True)
                    self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                    with (campaign/f'results/foil_{foil}cm/central_dense_half/rg01/slopes.tsv').open() as stream:
                        rows=list(csv.DictReader(stream,delimiter='\t'))
                    rows=[r for r in rows if r['ndel']=='-1']
                    self.assertEqual(len(rows),2)
                    for row in rows:
                        self.assertEqual(int(row['n']),201 if row['subset']=='core' else 101)
                        expected=0 if arm=='SHMS' else 1000/(168-foil*math.cos(math.radians(12.49)))
                        self.assertAlmostEqual(float(row['slope_mrad_per_cm']),expected,places=6)
                        self.assertAlmostEqual(float(row['predicted_slope']),expected,places=6)


if __name__=='__main__':
    unittest.main()
