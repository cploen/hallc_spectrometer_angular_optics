"""Targeted annotation checks; these do not claim native Geant4 rendering."""
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

HERE=Path(__file__).resolve().parents[1]


class VisualizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler=shutil.which(os.environ.get('CXX','c++'))
        if not compiler: raise unittest.SkipTest('C++ compiler required for annotation checks')
        cls.tmp=tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.exe=Path(cls.tmp.name)/'annotation_probe'
        result=subprocess.run([compiler,'-std=c++17','-I',str(HERE/'src'),
            str(HERE/'tests/annotation_probe.cc'),'-o',str(cls.exe)],text=True,capture_output=True)
        if result.returncode: raise RuntimeError(result.stderr)
        cls.layers=cls.probe('--layers').splitlines()

    @classmethod
    def probe(cls,*args,input=''):
        return subprocess.run([str(cls.exe),*args],input=input,text=True,capture_output=True,check=True).stdout

    def test_existing_primitives_partition_without_duplication(self):
        names=[line.split()[1] for line in (HERE/'generated/scene.tsv').read_text().splitlines()]
        rows=[line.split('\t') for line in self.probe(input='\n'.join(names)).splitlines()]
        self.assertEqual([r[0] for r in rows],names)
        self.assertEqual(Counter(r[1] for r in rows),{
            'HMS_base':98,'HMS_lab_frame':3,'HMS_transport_frame':3,'HMS_sieve_frame':3,'HMS_rays':7})

    def test_saved_aliases_do_not_merge_constructed_quantities(self):
        rows=[line.split('\t') for line in self.probe('--aliases').splitlines()]
        aliases={r[1]:r[2] for r in rows}
        self.assertEqual(aliases,{
            'xtar':'H.gtr.x','ytar':'H.gtr.y','xptar':'H.gtr.th','yptar':'H.gtr.ph',
            'xsieve':'H.extcor.xsieve','ysieve':'H.extcor.ysieve',
            'reactx':'H.react.x','reacty':'H.react.y','reactz = ztar':'H.react.z',
            'xbpm_tar':'H.rb.raster.fr_xbpm_tar','ybpm_tar':'H.rb.raster.fr_ybpm_tar'})

    def test_presets_switch_layers_without_camera_changes(self):
        expected={
            'geometry':{'HMS_labels'}, 'frames':{'HMS_labels','HMS_lab_frame','HMS_transport_frame','HMS_sieve_frame'},
            'lab_frame':{'HMS_labels','HMS_lab_frame'}, 'hms_frame':{'HMS_labels','HMS_transport_frame'},
            'sieve_frame':{'HMS_labels','HMS_sieve_frame'},'hcana':{'HMS_code_names'},
            'branches':{'HMS_branch_names'},'aliases':{'HMS_code_names','HMS_branch_names'},
            'rays':{'HMS_labels','HMS_rays'},'everything':set(self.layers)-{'HMS_base'}}
        def apply(name,state):
            for raw in (HERE/'macros'/name).read_text().splitlines():
                parts=raw.split()
                if not parts or parts[0].startswith('#'): continue
                command=parts[0]
                if command=='/control/execute':
                    self.assertEqual(parts[1],'reset_layers.mac');apply(parts[1],state)
                elif command=='/vis/scene/activateModel':
                    matches=[key for key in self.layers if parts[1] in key]
                    self.assertEqual(len(matches),1,raw)
                    state[matches[0]]=parts[2]=='true'
                else:
                    self.assertIn(command,{'/vis/viewer/set/autoRefresh','/vis/viewer/refresh'})
        for name,on in expected.items():
            state={name:True for name in self.layers}
            apply('view_'+name+'.mac',state)
            expected_visible=on if name in {'hcana','branches','aliases'} else on|{'HMS_base'}
            self.assertEqual({key for key,value in state.items() if value},expected_visible)
        cmake=(HERE/'CMakeLists.txt').read_text()
        for name in expected: self.assertIn('view_'+name+'.mac',cmake)

    def test_geometry_sources_match_frozen_validation_provenance(self):
        # This visualization task must not alter coordinate/source equations.
        report=json.loads((HERE/'generated/validation.json').read_text())
        repo=HERE.parents[1]
        for relative,expected in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256((repo/relative).read_bytes()).hexdigest(),expected,relative)
        self.assertEqual(hashlib.sha256((HERE/'config/geometry.json').read_bytes()).hexdigest(),report['config_sha256'])


if __name__=='__main__': unittest.main()
