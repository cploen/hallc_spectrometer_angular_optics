"""Independent geometric rays; all vectors are ordered (x,y,z), lengths cm.

The JSON defines a declared benchmark, not a silently completed survey.
Fit and saved-replay constructions live in comparisons.py, not in this model.
"""
from dataclasses import dataclass
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def add(a, b): return tuple(x+y for x, y in zip(a, b))
def sub(a, b): return tuple(x-y for x, y in zip(a, b))
def scale(a, k): return tuple(k*x for x in a)
def dot(a, b): return sum(x*y for x, y in zip(a, b))
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a): return math.sqrt(dot(a, a))


def load_config(path=HERE/'config/geometry.json', setting=None):
    cfg = json.loads(Path(path).read_text())
    chosen=setting or cfg['default_setting']
    if chosen not in cfg['settings']:
        raise ValueError('Unknown setting: '+chosen)
    cfg['scenario']=cfg['settings'][chosen]
    cfg['selected_setting']=chosen
    if cfg['schema'] != 'hms_coordinate_audit_v1' or cfg['length_unit'] != 'cm' or cfg['angle_unit'] != 'deg':
        raise ValueError('Require documented cm/degree coordinate-audit schema')
    if cfg['scenario']['pointing_mode'] not in ('transport_components', 'podd_lab_components'):
        raise ValueError('Pointing frame must be explicit')
    if cfg['scenario']['mispointing_model'] not in ('HMS_code_defaults','explicit_translation'):
        raise ValueError('Unrecognized pointing model; do not silently supply defaults')
    if cfg['scenario']['nominal_target_lab_cm'] != [0., 0., 0.]:
        raise ValueError('This schema defines the nominal target as lab origin; do not shift it silently')
    if cfg['validation']['fit_benchmark_lab_x_equals_minus_xbpm'] is not True:
        raise ValueError('Synthetic B comparison requires the explicit X=-xbpm hypothesis')
    if cfg['sieve']['material_solid_enabled']:
        raise ValueError('Physical sieve solid is gated on unresolved geometry')
    if cfg['sieve']['projection_distance_cm'] <= max(cfg['validation']['foil_test_z_cm']):
        raise ValueError('Sieve must be downstream of test vertices')
    configured_pointing(cfg)  # Validate explicit per-setting translations now.
    return cfg


def hms_default_mispointing(angle_magnitude_deg):
    """Source DEFAULT translations, not measured hardware offsets."""
    ax = min(abs(angle_magnitude_deg), 50.)
    ay = min(abs(angle_magnitude_deg), 40.)
    return (.1*(2.37-.086*ax+.0012*ax*ax),
            .1*(.52-.012*ay+.002*ay*ay), 0.)


@dataclass(frozen=True)
class Frame:
    name: str
    origin_lab_cm: tuple
    axes_lab: tuple  # columns of the active local-to-lab rotation

    def vector_to_lab(self, v):
        return tuple(sum(self.axes_lab[j][i]*v[j] for j in range(3)) for i in range(3))

    def to_lab(self, p):
        return add(self.origin_lab_cm, self.vector_to_lab(p))

    def from_lab(self, p):
        d = sub(p, self.origin_lab_cm)
        return tuple(dot(a, d) for a in self.axes_lab)


def transport_frame(alpha_deg, phi_deg, pointing_cm, pointing_mode):
    """Podd SetCentralAngles rotation, bend_down=false, |phi|<90 degrees.

    alpha: signed geographic angle toward +lab X (HMS negative).
    Pointing interpretation is explicitly selected; no frame is inferred.
    """
    if not all(math.isfinite(x) for x in (alpha_deg, phi_deg, *pointing_cm)) or abs(phi_deg) >= 90:
        raise ValueError('Invalid frame parameters')
    a, f = math.radians(alpha_deg), math.radians(phi_deg)
    s, c, sf, cf = math.sin(a), math.cos(a), math.sin(f), math.cos(f)
    axes = ((s*sf, -cf, c*sf), (c, 0., -s), (s*cf, sf, c*cf))
    zero = Frame('HMS_TRANSPORT', (0., 0., 0.), axes)
    if pointing_mode == 'transport_components':
        origin = zero.vector_to_lab(pointing_cm)
    elif pointing_mode == 'podd_lab_components':
        origin = tuple(pointing_cm)
    else:
        raise ValueError('Unknown pointing interpretation')
    return Frame('HMS_TRANSPORT', origin, axes)


def configured_frame(cfg):
    s = cfg['scenario']
    angle = s['central_angle_magnitude_deg']
    return transport_frame(-abs(angle), s['out_of_plane_deg'],
                           configured_pointing(cfg), s['pointing_mode'])


def configured_pointing(cfg):
    s=cfg['scenario']
    if s['mispointing_model']=='HMS_code_defaults':
        return hms_default_mispointing(s['central_angle_magnitude_deg'])
    m=s['pointing_translation_cm']
    if m is None or len(m)!=3 or not all(math.isfinite(v) for v in m):
        raise ValueError('Explicit per-setting translation requires three finite cm values')
    if m[2]!=0:
        raise ValueError('Longitudinal pointing is not supported by the inspected HCANA setup')
    return tuple(m)


def hole(cfg, i, j):
    s = cfg['sieve']
    if not 0 <= i < s['nx'] or not 0 <= j < s['ny']:
        raise ValueError('Hole index outside nominal grid')
    return ((i-(s['nx']-1)/2)*s['pitch_x_cm'],
            (j-(s['ny']-1)/2)*s['pitch_y_cm'], s['projection_distance_cm'])


@dataclass(frozen=True)
class Ray:
    xt_cm: float
    yt_cm: float
    p: float
    q: float

    def at_z(self, z_cm):
        return (self.xt_cm+z_cm*self.p, self.yt_cm+z_cm*self.q, z_cm)


def endpoint_ray(vertex, endpoint):
    """Direct Euclidean construction; does not import any optics equations."""
    if not all(math.isfinite(x) for x in (*vertex, *endpoint)):
        raise ValueError('Nonfinite endpoint')
    direction = sub(endpoint, vertex)
    if abs(direction[2]) < 1e-12:
        raise ValueError('Line is parallel to target plane')
    p, q = direction[0]/direction[2], direction[1]/direction[2]
    return Ray(vertex[0]-vertex[2]*p, vertex[1]-vertex[2]*q, p, q)


def line_plane_intersection(point, direction, plane_origin, normal):
    """Independent vector-plane form used to check the slope/intercept form."""
    denom = dot(direction, normal)
    if abs(denom) < 1e-12:
        raise ValueError('Line parallel to plane')
    return add(point, scale(direction, dot(sub(plane_origin, point), normal)/denom))
