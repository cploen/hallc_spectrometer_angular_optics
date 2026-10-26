"""One standardized rotation and an ordered polynomial X band family.

The input columns are ``[xpfp, xfp]``. Integer latent coordinates identify
possible sieve bands; they do not assert that a band is populated.
"""
import numpy as np
from scipy.optimize import least_squares


def _slot_count(max_slots):
    try:
        count = int(max_slots)
    except (TypeError, ValueError, OverflowError):
        raise ValueError('max_slots must be an integer >= 3') from None
    if isinstance(max_slots, (bool, np.bool_)) or count != max_slots or count < 3:
        raise ValueError('max_slots must be an integer >= 3')
    return count


def _raw_array(raw, minimum=0):
    raw = np.asarray(raw, dtype=float)
    if raw.ndim != 2 or raw.shape[1] != 2 or len(raw) < minimum:
        raise ValueError('Expected a finite N by 2 sample with >= %d events' % minimum)
    if not np.isfinite(raw).all():
        raise ValueError('Input coordinates must all be finite')
    return raw


def _design(j, p, j_degree, p_degree, p_scale, max_slots=9):
    half_span = (_slot_count(max_slots) - 1.) / 2.
    if not np.isfinite(p_scale) or p_scale <= 0:
        raise ValueError('The longitudinal coordinate scale must be positive and finite')
    jj = (np.asarray(j) - half_span) / half_span
    pp = np.asarray(p) / p_scale
    return np.column_stack([jj**r * pp**s
                            for r in range(j_degree + 1)
                            for s in range(p_degree + 1)])


def _fit(j, p, q, w, j_degree, p_degree, p_scale, max_slots=9):
    design = _design(j, p, j_degree, p_degree, p_scale, max_slots)
    initial = np.linalg.lstsq(design * w[:, None], q * w, rcond=None)[0]
    residual = design @ initial - q
    robust_scale = max(.002, float(1.4826 * np.median(np.abs(residual - np.median(residual)))))
    fitted = least_squares(lambda x: (design @ x - q) * w, initial,
                           loss='soft_l1', f_scale=robust_scale, max_nfev=250)
    if not np.isfinite(fitted.x).all():
        raise ValueError('Nonfinite polynomial family coefficients')
    return fitted.x


def _count(fan):
    return _slot_count(fan.get('max_slots', len(fan.get('possible_indices', range(9)))))


def global_centers(fan, p, indices=None):
    """Return q(j,p), shape (number of indices, number of p coordinates)."""
    p = np.asarray(p, float).reshape(-1)
    count = _count(fan)
    j = np.asarray(np.arange(count) if indices is None else indices, float).reshape(-1)
    jj, pp = np.broadcast_arrays(j[:, None], p[None, :])
    result = _design(jj.ravel(), pp.ravel(), fan['j_degree'], fan['p_degree'],
                     fan['p_scale'], count) @ np.asarray(fan['coefficients'])
    return result.reshape(jj.shape)


def _evaluate_pairwise(fan, j, p):
    return _design(j, p, fan['j_degree'], fan['p_degree'], fan['p_scale'],
                   _count(fan)) @ np.asarray(fan['coefficients'])


def global_coordinates(fan, raw):
    """Rotate and invert q(j,p) to u, flagging points outside ordered support.

    The inverse domain includes half a band beyond each outer center. The
    ``finite_geometry`` flag is separate from the occupancy support flag.
    """
    raw = _raw_array(raw)
    mean = np.asarray(fan['mean'], float)
    scale = np.asarray(fan['scale'], float)
    if (mean.shape != (2,) or scale.shape != (2,) or not np.isfinite(mean).all()
            or not np.isfinite(scale).all() or np.any(scale <= 0)
            or not np.isfinite(fan['theta_deg'])):
        raise ValueError('Invalid fitted X standardization or rotation')
    z = (raw - mean) / scale
    theta = np.deg2rad(fan['theta_deg'])
    q = z @ np.array([np.cos(theta), np.sin(theta)])
    p = z @ np.array([-np.sin(theta), np.cos(theta)])
    lo = np.full(len(p), -.5)
    hi = np.full(len(p), _count(fan) - .5)
    step = .001
    dl = (_evaluate_pairwise(fan, lo + step, p) - _evaluate_pairwise(fan, lo - step, p)) / (2 * step)
    dh = (_evaluate_pairwise(fan, hi + step, p) - _evaluate_pairwise(fan, hi - step, p)) / (2 * step)
    valid = ((dl > 0) & (dh > 0) & (_evaluate_pairwise(fan, lo, p) <= q)
             & (_evaluate_pairwise(fan, hi, p) >= q))
    for _ in range(40):
        middle = (lo + hi) / 2
        right = _evaluate_pairwise(fan, middle, p) < q
        lo = np.where(right, middle, lo)
        hi = np.where(right, hi, middle)
    return dict(q=q, p=p, u=np.where(valid, (lo + hi) / 2, np.nan), finite_geometry=valid)
