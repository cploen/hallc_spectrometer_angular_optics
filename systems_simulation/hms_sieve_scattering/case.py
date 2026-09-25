"""Provisional material-simulation overlay on the tested HMS coordinate model."""
import hashlib
import json
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
