import sys
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from case import load_case  # noqa: E402


class ProvisionalCaseTests(unittest.TestCase):
    def test_overlay_preserves_reference_and_applies_per_angle_code_default(self):
        case, frame, pointing = load_case()
        self.assertEqual(case["coordinate_reference"]["config"], "../../hms_sieve/config/geometry.json")
        self.assertAlmostEqual(case["setting"]["nominal_central_angle_deg"], 12.5)
        self.assertAlmostEqual(pointing[0], 0.14825)
        self.assertAlmostEqual(pointing[1], 0.06825)
        self.assertAlmostEqual(frame.from_lab((0, 0, 0))[0], -pointing[0], places=12)

    def test_provisional_carbon_slab_preserves_reported_areal_density(self):
        case, _, _ = load_case()
        target = case["target"]
        self.assertEqual(target["name"], "Carbon 0.5%")
        self.assertAlmostEqual(target["thickness_in"], 0.032)
        self.assertAlmostEqual(target["effective_density_g_cm3"] * target["thickness_cm"],
                               target["areal_density_g_cm2"], places=12)
        self.assertIn("PROVISIONAL GUESS", target["thickness_status"])

    def test_beam_and_sieve_inputs_keep_unresolved_values_explicit(self):
        case, _, _ = load_case()
        self.assertEqual(case["beam"]["target_centroid_lab_xy_mm"], [0.22, -0.31])
        self.assertEqual(case["beam"]["intrinsic_sigma_lab_xy_mm"], [None, None])
        self.assertEqual(case["sieve"]["hole_pattern"]["orientation_to_hms_axes"],
                         "TBD; retain drawing front-view coordinates until the installed/replay view is verified")
        self.assertIsNone(case["sieve"]["longitudinal_placement_cm"])


if __name__ == "__main__":
    unittest.main()
