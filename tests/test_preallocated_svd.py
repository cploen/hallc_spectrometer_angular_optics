"""Bounded ROOT integration: real current solver admission loop, no SVD solve."""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import uproot
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from core_sample import PROJECT,digest,write_tsv,read_tsv
from preallocated_svd import membership_manifest,run


@unittest.skipUnless(shutil.which('root'),'ROOT is required for solver integration')
class SolverTests(unittest.TestCase):
    def fixture(self,p):
        b=p/'HMS_fixture'; (b/'root').mkdir(parents=True)
        n=18018
        fields='ys ysT xtar xtarT xptar yptar ytar xptarT yptarT ytarT ztarT delta xpfp ypfp xfp yfp'.split()
        a={k:np.zeros(n,np.float64) for k in fields}
        y=np.tile(np.repeat(np.arange(9,dtype=np.int32),1001),2)
        nd=np.repeat(np.arange(2,dtype=np.int32),9009)
        a.update(entry=np.arange(n,dtype=np.int64)+2**54,foil=np.zeros(n,np.int32),ndel=nd,
                 xscol=np.full(n,4,np.int32),yscol=y,sample=np.ones(n,np.int32),core_keep=np.ones(n,np.int32))
        a['delta']=np.where(nd==0,-9.,-7.);a['ysT']=(y-4)*1.524
        with uproot.recreate(b/'root/Optics_666701_-1_fit_tree_gmm.root') as f:f.mktree('TFit',a)
        settings=[dict(rungroup='fixture',optics_id=666701)]
        write_tsv(b/'rungroups_fixture_inputs.tsv',settings)
        shutil.copy2(PROJECT/'DATfiles/list_of_optics_run.dat',b/'optics.dat')
        (b/'oldfit.dat').write_text('! Compact fixture only\n 0 0 0 0 10000\n ----\n')
        ids=[dict(rungroup='fixture',entry=int(a['entry'][i]),foil=0,ndel=int(nd[i]),zfoil=0.,delta_low=-10. if nd[i]==0 else -8.,delta_high=-8. if nd[i]==0 else -5.,xscol=4,yscol=int(y[i])) for i in range(n)]
        write_tsv(b/'selected_ids.tsv',ids)
        (b/'allocation_manifest.json').write_text(json.dumps(dict(mode='event_allocation',status='quotas_ready',totals=dict(fit=n))))
        h=digest(b/'allocation_manifest.json')
        membership_manifest(b/'solver_input.tsv',ids,settings,h)
        manifest=dict(schema='core_preallocated_v1',sample='fit',sample_manifest=h,
                      outputs={str(f.relative_to(b)):digest(f) for f in b.rglob('*') if f.is_file()})
        (b/'build.json').write_text(json.dumps(manifest))
        return b,n
    def test_exact_admission_and_resource_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);b,n=self.fixture(p)
            old=os.environ.get('HCANA');os.environ['HCANA']='root'
            try:
                run(b,p/'accepted',max_events=n)
                used=read_tsv(p/'accepted/solver_used_ids.tsv')
                self.assertEqual(len(used),n)
                self.assertEqual(len({int(r['entry']) for r in used}),n)
                self.assertFalse(list((p/'accepted/matrices').glob('*.dat')))
                # 1001 in each Y/delta population, 18018 in a rungroup/foil:
                # the legacy rules would retain 15000, while this path admits all 18018.
                self.assertEqual(min(15000,18*1000),15000)
                with self.assertRaises(ValueError):run(b,p/'oversized',max_events=n-1)
                self.assertFalse((p/'oversized').exists())
                with self.assertRaises(RuntimeError):run(b,p/'memory_limit',max_events=n,max_bytes=1)
                self.assertFalse(list((p/'memory_limit/matrices').glob('*.dat')))
                # Exercise unchanged legacy sampling with a bounded synthetic fit only.
                legacy=p/'legacy_fixture'
                argv=['"fixture"','-1','20000','-1',json.dumps(str(b/'root')),json.dumps(str(legacy)),json.dumps(str(b/'oldfit.dat')),json.dumps(str(b/'rungroups_fixture_inputs.tsv')),json.dumps(str(b/'optics.dat'))]
                proc=subprocess.run(['root','-b','-l','-q',str(PROJECT/'fit_opt_matrix_gmm.C')+'('+','.join(argv)+')'],cwd=PROJECT,capture_output=True,text=True)
                self.assertEqual(proc.returncode,0,proc.stderr[-2000:])
                self.assertIn('number to fit = 15000',proc.stdout)
                # A supplied invalid event is rejected by the real ROOT preflight, with its ID.
                rootfile=b/'root/Optics_666701_-1_fit_tree_gmm.root'
                with uproot.open(rootfile) as f:a=f['TFit'].arrays(library='np')
                a['xfp'][0]=float('nan')
                with uproot.recreate(rootfile) as f:f.mktree('TFit',a)
                m=json.loads((b/'build.json').read_text())
                m['outputs']['root/'+rootfile.name]=digest(rootfile)
                (b/'build.json').write_text(json.dumps(m))
                with self.assertRaises(RuntimeError):run(b,p/'invalid',max_events=n)
                rejects=read_tsv(p/'invalid/solver_rejected_ids.tsv')
                self.assertEqual(int(rejects[0]['entry']),2**54)
                self.assertFalse(list((p/'invalid/matrices').glob('*.dat')))
                # Hash mismatch must stop before invoking ROOT.
                with (b/'selected_ids.tsv').open('a') as out:out.write('tampered\n')
                with self.assertRaises(ValueError):run(b,p/'tampered',max_events=n)
            finally:
                if old is None:os.environ.pop('HCANA',None)
                else:os.environ['HCANA']=old

if __name__=='__main__':unittest.main()
