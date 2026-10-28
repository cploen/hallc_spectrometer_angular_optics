"""Targeted campaign configuration and global-label controls (no ROOT)."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from joint_surface.campaign import apply_map, choose_global_map, label_group, load_config
from setup_joint_surface_campaign import create_joint_campaign


class CampaignTests(unittest.TestCase):
    def test_configuration_changes_cannot_publish_stale_predictions(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign = Path(temporary) / 'HMS_joint'
            base = campaign / '00_joint_surface_work/test'
            base.mkdir(parents=True)
            (base / 'geometry_frozen.json').write_text(json.dumps({'sha256': {}}))
            (base / 'fit_started.json').write_text(json.dumps({'config': {'coverage_min_events': 10}}))
            with self.assertRaisesRegex(ValueError, 'configuration changed'):
                label_group(campaign, {'rungroup': 'test'}, {'foils': [0.]},
                            {'coverage_min_events': 5}, None)

    def test_global_map_is_whole_family(self):
        latent = np.repeat(np.arange(5), 20)
        reconstructed = (6 - latent - 4) * 1.524 + .03
        mapping = choose_global_map(latent, reconstructed, 9, 1.524, lambda j: (j - 4) * 1.524)
        self.assertEqual((mapping['orientation'], mapping['origin']), (-1, 6))
        np.testing.assert_array_equal(apply_map(latent, mapping), 6 - latent)
        self.assertTrue(mapping['reconstructed_sieve_used'])
        manual = choose_global_map(latent, np.full(len(latent), np.nan), 9, 1.524,
                                   lambda j: (j - 4) * 1.524, {'orientation': 1, 'origin': 2})
        np.testing.assert_array_equal(apply_map(np.array([-1, 0, 4]), manual), [-1, 2, 6])
        self.assertFalse(manual['reconstructed_sieve_used'])
        with self.assertRaises(ValueError):
            choose_global_map(latent, reconstructed, 9, 1., float, {'orientation': 1, 'origin': 6})

    def test_fresh_setup_and_plan_need_no_replay_mount(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for arm in ('HMS', 'SHMS'):
                source, destination = root / (arm + '_stable'), root / (arm + '_joint')
                (source / 'config').mkdir(parents=True)
                table = source / 'config/rungroups_trial_inputs.tsv'
                table.write_text('rungroup\toptics_id\tangle\tfoils\tnruns\truns\trootfile\n'
                                 'test\t99002\t12.5\t0\t1\t123\t/not/mounted/replay.root\n')
                (source / 'config/oldfit.dat').write_text('seed\n')
                (source / 'config/pre_gmm_veto.tsv').write_text('historical reference\n')
                inherited = {
                    'core_sample.json': {'min_events': 10, 'seed': 12345, 'holdout': .25},
                    'core_balance.json': {'balance': {'foil': 'equal', 'hole': 1.5, 'delta': 1.5},
                                          'qa_min': 10},
                    'sieve_mask.json': {'blocked': [[2, 3], [5, 5]]},
                }
                for name, value in inherited.items():
                    (source / 'config' / name).write_text(json.dumps(value) + '\n')
                (source / 'config/optics.dat').write_text('99002,test,12.5,1,1,3,0\n0\n-10,0,22\n')
                before = {path.relative_to(source): path.read_bytes()
                          for path in source.rglob('*') if path.is_file()}
                marker = create_joint_campaign(source, destination)
                spec, config = load_config(destination)
                self.assertEqual(spec.name, arm)
                self.assertEqual(config['global_maps'], {})
                for name in (*inherited, 'optics.dat'):
                    expected = before[Path('config') / name]
                    self.assertEqual((destination / 'config' / name).read_bytes(), expected)
                    self.assertIn(name, marker['copied_configs'])
                    self.assertEqual(marker['source_sha256'][name], hashlib.sha256(expected).hexdigest())
                self.assertFalse((destination / 'config/geometric_selection.json').exists())
                self.assertFalse((destination / 'config/pre_gmm_veto.tsv').exists())
                self.assertTrue((destination / 'config/previous_selection_reference/pre_gmm_veto.tsv').exists())
                result = subprocess.run([sys.executable, str(REPO / 'run_joint_surface_campaign.py'),
                                         str(destination), '--plan'],
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('[candidates]', result.stdout)
                self.assertIn(str(destination / 'config/optics.dat'), result.stdout)
                self.assertNotIn('[gmm', result.stdout)
                self.assertFalse((destination / '00_joint_surface_work').exists())
                metadata = root / (arm + '_metadata.dat')
                metadata.write_text('99002,test,12.5,1,1,3,0\n0\n-10,0,23\n')
                rejected = subprocess.run([sys.executable, str(REPO / 'run_joint_surface_campaign.py'),
                                           str(destination), '--plan', '--metadata', str(metadata)],
                                          text=True, capture_output=True)
                self.assertNotEqual(rejected.returncode, 0)
                self.assertIn('snapshot', rejected.stderr.lower())
                self.assertEqual({path.relative_to(source): path.read_bytes()
                                  for path in source.rglob('*') if path.is_file()}, before)
                (destination / 'config/pre_gmm_veto.tsv').write_text('# Legacy cut names are incompatible\nold_cut\n')
                with self.assertRaisesRegex(ValueError, 'legacy TCutG vetoes'):
                    load_config(destination)
                with self.assertRaises(FileExistsError):
                    create_joint_campaign(source, destination)

    def test_setup_validates_explicit_metadata_before_creating_destination(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, destination = root / 'HMS_stable', root / 'HMS_joint'
            (source / 'config').mkdir(parents=True)
            (source / 'config/rungroups_trial_inputs.tsv').write_text(
                'rungroup\toptics_id\tangle\tfoils\tnruns\truns\trootfile\n'
                'test\t99002\t12.5\t0\t1\t123\t/not/mounted/replay.root\n')
            (source / 'config/oldfit.dat').write_text('seed\n')
            (source / 'config/optics.dat').write_text('99001,other,12.5,1,1,3,0\n0\n-10,0,22\n')
            with self.assertRaises(ValueError):
                create_joint_campaign(source, destination)
            self.assertFalse(destination.exists())
            metadata = root / 'explicit.dat'
            metadata.write_text('99002,test,12.5,1,1,3,0\n0\n-10,0,22\n')
            create_joint_campaign(source, destination, metadata=metadata)
            self.assertEqual((destination / 'config/optics.dat').read_bytes(), metadata.read_bytes())


if __name__ == '__main__':
    unittest.main()
