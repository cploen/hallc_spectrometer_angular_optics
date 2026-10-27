"""Targeted campaign configuration and global-label controls (no ROOT)."""
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
                before = table.read_bytes()
                create_joint_campaign(source, destination)
                spec, config = load_config(destination)
                self.assertEqual(spec.name, arm)
                self.assertEqual(config['global_maps'], {})
                self.assertFalse((destination / 'config/geometric_selection.json').exists())
                self.assertFalse((destination / 'config/pre_gmm_veto.tsv').exists())
                self.assertTrue((destination / 'config/previous_selection_reference/pre_gmm_veto.tsv').exists())
                metadata = root / (arm + '_metadata.dat')
                metadata.write_text('99002,test,12.5,1,1,3,0\n0\n-10,0,22\n')
                result = subprocess.run([sys.executable, str(REPO / 'run_joint_surface_campaign.py'),
                                         str(destination), '--plan', '--metadata', str(metadata)],
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('[candidates]', result.stdout)
                self.assertNotIn('[gmm', result.stdout)
                self.assertFalse((destination / '00_joint_surface_work').exists())
                self.assertEqual(table.read_bytes(), before)
                (destination / 'config/pre_gmm_veto.tsv').write_text('# Legacy cut names are incompatible\nold_cut\n')
                with self.assertRaisesRegex(ValueError, 'legacy TCutG vetoes'):
                    load_config(destination)
                with self.assertRaises(FileExistsError):
                    create_joint_campaign(source, destination)


if __name__ == '__main__':
    unittest.main()
