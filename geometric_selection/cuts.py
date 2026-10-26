"""Convert numerical partitions to editable, raw-coordinate TCutG polygons.

X curves are sampled ordered strips. Y uses contours of each Gaussian's
log-density advantage over every competitor, clipped to observed support.
Disconnected regions remain separate components; holes use a doubled bridge
(the even/odd polygon rule used by ROOT). Export is checked on actual events.
"""
import numpy as np
from matplotlib.path import Path
from .x_geometry import global_centers, global_coordinates
from .y_geometry import _rotate, _logp


def to_raw(model, qp):
    theta = np.deg2rad(model['theta_deg'])
    c, s = np.cos(theta), np.sin(theta)
    return np.asarray(qp) @ np.array([[c, s], [-s, c]]) * model['scale'] + model['mean']


def build_x_polygons(geometry, occupancy, raw, samples=513, boundaries_u=None):
    indices = np.asarray(occupancy['accepted_indices'], int)
    if not len(indices) or np.any(np.diff(indices) <= 0):
        raise ValueError('X occupancy must contain increasing accepted indices')
    count = len(geometry['possible_indices'])
    inner = (indices[1:] + indices[:-1])/2 if boundaries_u is None else np.asarray(boundaries_u, float)
    bounds = np.r_[-.5, inner, count-.5]
    if len(bounds) != len(indices)+1 or not np.isfinite(bounds).all() or np.any(np.diff(bounds) <= 0):
        raise ValueError('X boundaries_u must be increasing internal dividers, one fewer than bands')
    p = global_coordinates(geometry, raw)['p']
    padding = max(1e-8, float(np.ptp(p))*1e-6)
    pp = np.linspace(float(p.min())-padding, float(p.max())+padding, samples)
    curves = global_centers(geometry, pp, indices=bounds)
    if np.any(np.diff(curves, axis=0) <= 0):
        raise ValueError('Exported X boundaries cross within the target support')
    polygons = []
    for k, index in enumerate(indices):
        qp = np.r_[np.c_[curves[k], pp], np.c_[curves[k+1, ::-1], pp[::-1]]]
        points = to_raw(geometry, np.r_[qp, qp[:1]])
        polygons.append(dict(index=int(index), points=points.tolist()))
    return polygons


def _bridge_rings(path):
    # Path's rendering simplifier uses coordinate units; it must never alter cuts.
    path.should_simplify = False
    rings = path.to_polygons(closed_only=True)
    if not rings:
        return None
    # Matplotlib returns the exterior followed by its oppositely oriented holes.
    joined = list(rings[0])
    for hole in rings[1:]:
        joined.extend(hole)
        joined.append(rings[0][0])
    return np.asarray(joined)


def build_y_polygons(model, raw, grid=401):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    qp = _rotate(raw, np.asarray(model['mean']), np.asarray(model['scale']), model['theta_deg'])
    lo, hi = qp.min(0), qp.max(0)
    pad = np.maximum((hi-lo)*.01, 1e-6)
    q, p = [np.linspace(lo[k]-pad[k], hi[k]+pad[k], grid) for k in (0, 1)]
    Q, P = np.meshgrid(q, p)
    logs = _logp(np.c_[Q.ravel(), P.ravel()], np.asarray(model['means']),
                 np.asarray(model['covariances']), np.asarray(model['weights']))
    polygons = []
    fig, ax = plt.subplots()
    try:
        for index in range(model['band_count']):
            margin = logs[:, index] - np.max(np.delete(logs, index, axis=1), axis=1)
            if margin.max() <= 0:
                continue
            contour = ax.contourf(q, p, margin.reshape(Q.shape), levels=[0., float(margin.max())+1.])
            # Matplotlib 3.x/3.10 expose paths through different interfaces.
            paths = contour.collections[0].get_paths() if hasattr(contour, 'collections') else contour.get_paths()
            for path in paths:
                path.should_simplify = False
                # A newer Matplotlib may put several exteriors into one Path.
                # Split at MOVETO only when that subpath has the exterior sign.
                rings = path.to_polygons(closed_only=True)
                current = []
                for ring in rings:
                    area = np.sum(ring[:-1, 0]*ring[1:, 1]-ring[1:, 0]*ring[:-1, 1])
                    if current and area > 0:
                        points = _bridge_rings(Path.make_compound_path(*[Path(r) for r in current]))
                        polygons.append(dict(index=index, points=to_raw(model, points).tolist()))
                        current = []
                    current.append(ring)
                if current:
                    points = _bridge_rings(Path.make_compound_path(*[Path(r) for r in current]))
                    polygons.append(dict(index=index, points=to_raw(model, points).tolist()))
            ax.clear()
    finally:
        plt.close(fig)
    if not polygons:
        raise ValueError('No Y decision regions could be exported')
    return polygons


def polygon_predictions(polygons, raw):
    """Evaluate exported polygons; same-index pieces are ORed for model QA."""
    raw = np.asarray(raw)
    membership = {}
    for polygon in polygons:
        index = int(polygon['index'])
        inside = Path(np.asarray(polygon['points'])).contains_points(raw)
        membership[index] = membership.get(index, np.zeros(len(raw), bool)) | inside
    labels = np.full(len(raw), -1, int)
    matches = np.zeros(len(raw), int)
    for index, inside in membership.items():
        labels[inside] = index
        matches += inside
    return labels, matches
