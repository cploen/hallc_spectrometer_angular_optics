"""Provisional material-simulation overlay on the tested HMS coordinate model."""
import hashlib
import json
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
COORD_DIR = REPO / "systems_simulation" / "hms_sieve"
BASE_CONFIG = COORD_DIR / "config" / "geometry.json"
sys.path.insert(0, str(COORD_DIR))
from geometry import hms_default_mispointing, transport_frame  # noqa: E402


def load_case(path=HERE / "config" / "case_12p5_provisional.json"):
    case_path = Path(path).resolve()
    case = json.loads(case_path.read_text())
    reference_path = (case_path.parent / case["coordinate_reference"]["config"]).resolve()
    if reference_path != BASE_CONFIG.resolve():
        raise ValueError("The overlay must use the established HMS coordinate configuration")
    expected = case["coordinate_reference"]["sha256"]
    actual = hashlib.sha256(reference_path.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError("The tested coordinate reference changed; review and revalidate before using this overlay")
    if case["schema"] != "hms_sieve_scattering_case_v1":
        raise ValueError("Unsupported scattering case schema")
    target = case["target"]
    if abs(target["areal_density_g_cm2"] / target["thickness_cm"] -
           target["effective_density_g_cm3"]) > 1e-12:
        raise ValueError("Target effective density no longer preserves the declared areal density")
    angle = case["setting"]["nominal_central_angle_deg"]
    offset = hms_default_mispointing(angle)
    frame = transport_frame(-abs(angle), case["setting"]["out_of_plane_deg"],
                            offset, case["setting"]["pointing_mode"])
    return case, frame, offset


def write_geant4_config(output, case_path=HERE / "config" / "case_12p5_provisional.json"):
    case, frame, _ = load_case(case_path)
    beam, target, sieve = case["beam"], case["target"], case["sieve"]
    q1, selection = case["q1_entrance"], case["momentum_selection"]
    center = frame.to_lab((0.0, 0.0, sieve["longitudinal_placement_cm"]))
    score = frame.to_lab((0.0, 0.0, q1["scoring_plane_hms_z_cm"]))
    values = {
        "events": beam["event_count"],
        "beam_x_cm": beam["target_centroid_lab_xy_mm"][0] / 10.0,
        "beam_y_cm": beam["target_centroid_lab_xy_mm"][1] / 10.0,
        "raster_x_half_cm": beam["raster_half_width_xy_mm"][0] / 10.0,
        "raster_y_half_cm": beam["raster_half_width_xy_mm"][1] / 10.0,
        "source_z_cm": beam["source_plane_lab_z_cm"],
        "beam_energy_GeV": beam["incident_energy_GeV"],
        "target_thickness_cm": target["thickness_cm"],
        "target_density_g_cm3": target["effective_density_g_cm3"],
        "target_half_x_cm": target["simulation_patch_half_width_cm"][0],
        "target_half_y_cm": target["simulation_patch_half_width_cm"][1],
        "sieve_center_x_lab_cm": center[0],
        "sieve_center_y_lab_cm": center[1],
        "sieve_center_z_lab_cm": center[2],
        "sieve_half_x_cm": sieve["plate_envelope_horizontal_cm"] / 2.0,
        "sieve_half_y_cm": sieve["plate_envelope_vertical_cm"] / 2.0,
        "sieve_half_z_cm": sieve["plate_thickness_cm"] / 2.0,
        "sieve_density_g_cm3": sieve["material_density_g_cm3"],
        "pitch_x_cm": sieve["hole_pattern"]["horizontal_pitch_cm"],
        "pitch_y_cm": sieve["hole_pattern"]["vertical_pitch_cm"],
        "hole_diameter_cm": sieve["hole_pattern"]["large_hole_diameter_cm"],
        "center_hole_diameter_cm": sieve["hole_pattern"]["reduced_center_hole_diameter_cm"],
        "grid_top_offset_cm": sieve["hole_pattern"]["drawing_grid_center_offset_magnitude_cm"],
        "q1_score_x_lab_cm": score[0],
        "q1_score_y_lab_cm": score[1],
        "q1_score_z_lab_cm": score[2],
        "q1_aperture_radius_cm": q1["vacuum_vessel_inner_radius_cm"],
        "hms_p0_GeV_c": selection["hms_central_momentum_GeV_c"],
        "hms_delta_half_width": selection["delta_half_width"],
    }
    for axis_name, axis in zip(("hms_x", "hms_y", "hms_z"), frame.axes_lab):
        for component, value in zip(("x", "y", "z"), axis):
            values[f"{axis_name}_{component}"] = value
    for i, j in sieve["blocked_positions_grid_indices"]:
        values[f"blocked_{i}_{j}"] = 1
    Path(output).write_text("".join(f"{key} {value:.12g}\n" for key, value in values.items()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Geant4 inputs from the frozen geometry plus provisional case")
    parser.add_argument("--write-g4-config", required=True, help="output flat key/value configuration file")
    args = parser.parse_args()
    write_geant4_config(args.write_g4_config)
