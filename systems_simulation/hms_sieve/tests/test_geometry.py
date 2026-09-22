import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from geometry import *
from comparisons import fit_target_ray, hcana_geometry_step, extcor_xtar_update
from build_geometry import compare, scene


class GeometryTests(unittest.TestCase):
    def setUp(self):
        self.cfg=load_config(); self.frame=configured_frame(self.cfg)

    def test_rotation_is_right_handed_and_inverse_closes(self):
        for angle in [-40.,-12.490,0.,12.490]:
            for phi in [-1.,0.,1.]:
                for mode in ['transport_components','podd_lab_components']:
                    f=transport_frame(angle,phi,(.15,.07,0.),mode)
                    self.assertAlmostEqual(dot(cross(f.axes_lab[0],f.axes_lab[1]),f.axes_lab[2]),1.,14)
                    for p in [(0,0,0),(.2,-.3,8.),(2,3,168)]:
                        for a,b in zip(p,f.from_lab(f.to_lab(p))): self.assertAlmostEqual(a,b,12)

    def test_lab_signs_from_source_rotation(self):
        f=transport_frame(-12.490,0,(0,0,0),'transport_components')
        self.assertEqual(f.from_lab((0,1,0)),(-1.,0.,0.))
        self.assertGreater(f.from_lab((0,0,1))[1],0)
        self.assertLess(f.axes_lab[2][0],0)  # HMS central ray points toward lab -X
        self.assertNotEqual(f.from_lab((.2,0,0))[2],0)  # horizontal beam affects z_v

    def test_reference_rays_close_all_requested_cases(self):
        rows=compare(self.cfg)
        refs=[r for r in rows if r['construction']=='A_endpoint']
        self.assertEqual(len(refs),84)
        for r in refs:
            for k in ['dx_cm','dy_cm']: self.assertLess(abs(r[k]),1e-12)
        self.assertTrue(any(r['lab_X_cm'] for r in refs))
        self.assertTrue(any(r['lab_Y_cm'] for r in refs))
        self.assertEqual({r['foil_lab_z_cm'] for r in refs},{-8.,0.,8.})

    def test_planes_are_not_interchangeable(self):
        f=self.frame
        # Lab Z=0 contains a horizontally displaced point, generally not on z_T=0.
        self.assertGreater(abs(f.from_lab((.2,0,0))[2]),.01)
        # Transverse pointing changes the coordinate origin, not this plane normal.
        self.assertAlmostEqual(dot(f.origin_lab_cm,f.axes_lab[2]),0.,14)

    def test_unmodified_fit_mismatch_is_preserved(self):
        mx,my,_=hms_default_mispointing(12.490)
        for z in [-8,0,8]:
            for bx,by in [(0,0),(.2,.2)]:
                r=fit_target_ray(12.490,z,2.54,1.524,by,-bx,168.)
                self.assertAlmostEqual(r.at_z(168.)[0]-2.54,-by-mx,13)
                self.assertAlmostEqual(r.at_z(168.)[1]-1.524,bx*r.q*math.sin(math.radians(12.490)),13)

    def test_actual_cpp_fit_targets_match_comparison(self):
        compiler=shutil.which(os.environ.get('CXX', 'c++'))
        if not compiler: self.skipTest('C++ compiler unavailable for source parity check')
        repo=HERE.parents[1]
        if not (repo/'spectrometer_config.h').exists(): self.skipTest('Parent optics source unavailable')
        cases=[(12.490,z,2.54,1.524, .7, by,-bx)
               for z in [-8,0,8] for bx,by in [(0,0),(.2,.2),(-.2,-.2)]]
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'fit_probe'
            built=subprocess.run([compiler,'-std=c++17','-I',str(repo),str(HERE/'tests/fit_probe.cc'),'-o',str(exe)],text=True,capture_output=True)
            self.assertEqual(built.returncode, 0, built.stderr)
            result=subprocess.run([str(exe)],input='\n'.join(' '.join(map(str,c)) for c in cases),text=True,capture_output=True,check=True)
        for c,line in zip(cases,result.stdout.splitlines()):
            a,z,x,y,rx,ry,bpm=c;r=fit_target_ray(a,z,x,y,ry,bpm,168.)
            for got,want in zip(map(float,line.split()),[r.xt_cm,r.yt_cm,r.p,r.q]): self.assertAlmostEqual(got,want,14)
        self.assertEqual(len(result.stdout.splitlines()),len(cases))

    def test_hcana_interface_old_new_slopes_and_units(self):
        v=(.2,.2,8);a=endpoint_ray(self.frame.from_lab(v),hole(self.cfg,6,6))
        theta=math.radians(self.cfg['scenario']['central_angle_magnitude_deg'])
        mx=configured_pointing(self.cfg)[0]
        for pold in [a.p,a.p+.001]:
            c=hcana_geometry_step(v,mx,theta,pold,a.p,a.yt_cm/100,a.q)
            self.assertAlmostEqual(c.at_z(168)[0]-a.at_z(168)[0],
                8*math.cos(theta)*(a.p-pold)-.2*a.p*math.sin(theta),13)
            self.assertAlmostEqual(c.at_z(168)[1],a.at_z(168)[1],13)
        self.assertAlmostEqual(extcor_xtar_update(.2,8,.15,.001,0),-.358)

    def test_one_offset_per_setting_does_not_hide_displacement_terms(self):
        for name in self.cfg['settings']:
            cfg=load_config(setting=name);rows=compare(cfg)
            for label in ['B_fit_target','C_HCANA_coordinate_step']:
                selected=[r for r in rows if r['construction']==label]
                self.assertEqual(len({(r['setting_offset_x_cm'],r['setting_offset_y_cm']) for r in selected}),1)
                for r in selected:
                    L=cfg['sieve']['projection_distance_cm']
                    self.assertAlmostEqual(r['dx_cm']-r['dxt_cm'],L*r['dp_slope'],12)
                    self.assertAlmostEqual(r['dy_cm']-r['dyt_cm'],L*r['dq_slope'],12)
                    if label=='B_fit_target':
                        self.assertAlmostEqual(r['dx_after_offset_cm'],-r['lab_Y_cm'],12)
            self.assertLess(max(abs(r[k]) for r in rows if r['construction']=='A_endpoint' for k in ['dx_cm','dy_cm']),1e-12)

    def test_pointing_can_be_adjusted_independently_per_setting(self):
        cfg=load_config();cfg['scenario']['mispointing_model']='explicit_translation'
        cfg['scenario']['pointing_translation_cm']=[.3,-.1,0]
        f=configured_frame(cfg)
        self.assertEqual(configured_pointing(cfg),(.3,-.1,0))
        self.assertAlmostEqual(f.from_lab((0,0,0))[0],-.3,14)
        self.assertNotEqual(configured_pointing(load_config(setting='HMS_15p195_code_defaults')),(.3,-.1,0))
        cfg['scenario']['pointing_translation_cm']=None
        with self.assertRaises(ValueError): configured_frame(cfg)

    def test_scene_points_use_validated_frame(self):
        data=scene(self.cfg)
        centers={m['name']:m for m in data['markers']}
        expected=self.frame.to_lab(hole(self.cfg,4,4))
        self.assertEqual(tuple(centers['nominal_hole_4_4']['position_lab_cm']),expected)
        self.assertIsNone(data['physical_sieve_solid'])
        self.assertEqual(len([m for m in data['markers'] if m['group']=='holes']),81)

    def test_refuse_unsupported_physical_geometry(self):
        cfg=json.loads(json.dumps(self.cfg));cfg['sieve']['material_solid_enabled']=True
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'bad.json';p.write_text(json.dumps(cfg))
            with self.assertRaises(ValueError): load_config(p)
        with self.assertRaises(ValueError): endpoint_ray((0,0,1),(1,0,1))

if __name__=='__main__': unittest.main()
