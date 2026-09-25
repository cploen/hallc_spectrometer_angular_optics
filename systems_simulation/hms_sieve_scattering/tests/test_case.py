import sys
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from case import load_case, write_geant4_config  # noqa: E402


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
        self.assertEqual(case["sieve"]["longitudinal_placement_cm"], 168.0)

    def test_runtime_inputs_use_provisional_beam_sieve_and_q1_dimensions(self):
        case, frame, _ = load_case()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "geant4.tsv"
            write_geant4_config(path)
            values = dict(line.split() for line in path.read_text().splitlines())
        self.assertEqual(values["events"], "10000")
        self.assertEqual(values["source_z_cm"], "-5")
        self.assertEqual(values["beam_energy_GeV"], "10.6")
        expected = frame.to_lab((0.0, 0.0, 168.0))
        for axis, value in zip(("x", "y", "z"), expected):
            self.assertAlmostEqual(float(values[f"sieve_center_{axis}_lab_cm"]), value, places=9)
        self.assertEqual(values["q1_aperture_radius_cm"], "20.05")
        self.assertEqual(case["momentum_selection"]["hms_central_momentum_GeV_c"], 6.667)
        self.assertEqual(values["blocked_2_3"], "1")
        self.assertEqual(values["blocked_5_5"], "1")

    def test_undeflected_lab_beam_misses_the_hms_sieve(self):
        case, frame, _ = load_case()
        x = case["beam"]["target_centroid_lab_xy_mm"][0] / 10.0
        y = case["beam"]["target_centroid_lab_xy_mm"][1] / 10.0
        normal = frame.axes_lab[2]
        plane_origin = frame.to_lab((0.0, 0.0, case["sieve"]["longitudinal_placement_cm"]))
        z = sum(normal[k] * plane_origin[k] for k in range(3))
        z = (z - normal[0] * x - normal[1] * y) / normal[2]
        local = frame.from_lab((x, y, z))
        self.assertAlmostEqual(local[2], 168.0, places=10)
        self.assertGreater(abs(local[1]), case["sieve"]["plate_envelope_vertical_cm"] / 2)


if __name__ == "__main__":
    unittest.main()
