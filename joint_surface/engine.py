"""Feature-only periodic discovery and regularized ordered joint surfaces.

Measured inputs are ``raw[:,0] = fp slope`` (radians), ``raw[:,1] = fp
position`` (cm), and event-level ``delta`` (percent). No reconstructed sieve
coordinates, physical row/column labels, or reference selections are accepted.
The learned integer coordinate is an ordered *latent* band index. Assigning its
physical sieve name is a separate, later operation.

This is a portable version of the successful 2026-10-27 feature-only pilot.
Its priors are smooth low-order surfaces, approximately periodic ordered bands,
one fixed standardized rotation, and separate observed-population checks.
"""
from copy import deepcopy
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.optimize import minimize
from scipy.special import ndtr, logsumexp

FAMILIES = ['bilinear', 'curved', 'static', 'shift', 'spacing', 'tilt', 'divergence']
DEFAULT_CONFIG = dict(
    theta_deg=135., phase_bins=[70, 70, 24], fft_bins=100, fft_padding=400,
    phase_starts=10, phase_iterations=100, phase_extra_penalty=.002,
    phase_curve_penalty=.00002, min_q_frequency=1.5, max_q_frequency=12.,
    minimum_validation_phase_power=.03, maximum_latent_span=24, max_slots=11,
    iterations=4, max_fit_events=18000, irls=6, extra_penalty=.002,
    high_order_penalty=.00002, seed_radius=.24, refine_radius=.30,
    coverage_min_events=10, coverage_bin_width=1., coverage_padding=.04,
    maximum_bridge=2., coverage_quantiles=[.005, .995],
    occupancy_cap=2000, occupancy_min_events=100, mixture_iterations=120,
    occupancy_likelihood_gain=10., occupancy_min_effective_events=8.,
    occupancy_min_coherent_intervals=2, occupancy_seed=32191,
    model_choice_fraction=.2, model_choice_tolerance=.01,
    density_q_bounds=[-8., 8.], families=FAMILIES, seed=86237)


def _configuration(config):
    cfg = deepcopy(DEFAULT_CONFIG)
    if config:
        unknown = set(config) - set(cfg)
        if unknown:
            raise ValueError('Unknown joint-surface configuration keys: ' + ', '.join(sorted(unknown)))
        cfg.update(deepcopy(config))
    if not 0 < cfg['model_choice_fraction'] < .5:
        raise ValueError('model_choice_fraction must be between 0 and .5')
    if cfg['coverage_bin_width'] <= 0 or cfg['maximum_bridge'] < 0:
        raise ValueError('Coverage bin width must be positive; maximum bridge nonnegative')
    if cfg['coverage_min_events'] < 1:
        raise ValueError('coverage_min_events must be positive')
    if cfg['max_slots'] < 2 or cfg['max_slots'] != int(cfg['max_slots']):
        raise ValueError('max_slots is a maximum hardware count, integer >=2')
    if not cfg['families'] or any(x not in FAMILIES for x in cfg['families']):
        raise ValueError('families must be a nonempty subset of FAMILIES')
    cfg['families'] = [f for f in FAMILIES if f in cfg['families']]
    return cfg


def _arrays(raw, delta, minimum=0):
    raw, delta = np.asarray(raw, float), np.asarray(delta, float)
    if raw.ndim != 2 or raw.shape[1] != 2 or delta.shape != (len(raw),) or len(raw) < minimum:
        raise ValueError('Expected raw N by 2 [slope,position] and delta length N, N >= %d' % minimum)
    return raw, delta


def _edges(edges):
    e = np.asarray(edges, float)
    if e.ndim != 1 or len(e) < 2 or not np.isfinite(e).all() or np.any(np.diff(e) <= 0):
        raise ValueError('delta_edges must be finite strictly increasing limits')
    return e


def make_frame(raw, theta_deg=135.):
    """Learn one standardization/rotation on geometry-training events only."""
    raw = np.asarray(raw, float)
    if raw.ndim != 2 or raw.shape[1] != 2 or len(raw) < 2 or not np.isfinite(raw).all():
        raise ValueError('Frame needs finite N by 2 coordinates, N >=2')
    scale = raw.std(0, ddof=1)
    if np.any(scale <= 0):
        raise ValueError('Both focal-plane input coordinates need nonzero variation')
    return dict(mean=raw.mean(0).tolist(), scale=scale.tolist(), theta_deg=float(theta_deg))


def rotate(model, raw):
    raw = np.asarray(raw, float)
    mean, scale = np.asarray(model['mean'], float), np.asarray(model['scale'], float)
    if mean.shape != (2,) or scale.shape != (2,) or np.any(scale <= 0) or not np.isfinite(np.r_[mean, scale, model['theta_deg']]).all():
        raise ValueError('Invalid fixed standardization/rotation')
    z = (raw - mean) / scale
    theta = np.deg2rad(model['theta_deg'])
    return z @ np.array([np.cos(theta), np.sin(theta)]), z @ np.array([-np.sin(theta), np.cos(theta)])


def _phase_design(q, r, d, extended):
    cols = [q, r, r*r, q*r, q*r*r, q*q]
    if extended:
        cols += [d, q*d, r*d]
    return np.column_stack(cols)


def _phase_values(seed, q, r, d):
    return (_phase_design((q-seed['qc'])/seed['qs'], (r-seed['rc'])/seed['rs'],
                          (d-seed['dcenter'])/seed['dscale'], seed['extended'])
            @ np.asarray(seed['coefficients']) - seed['phase'] - seed['offset'])


def _phase_power(seed, q, r, d):
    return float(abs(np.mean(np.exp(2j*np.pi*_phase_values(seed, q, r, d))))**2)


def _discover(raw, delta, frame, edges, cfg, extended):
    q, r = rotate(frame, raw)
    qlo, qhi = np.quantile(q, [.002, .998]); rlo, rhi = np.quantile(r, [.002, .998])
    qc, qs, rc, rs = (qlo+qhi)/2, (qhi-qlo)/2, (rlo+rhi)/2, (rhi-rlo)/2
    if min(qs, rs) <= 0:
        raise ValueError('Insufficient coordinate variation for periodic discovery')
    dc, ds = (edges[0]+edges[-1])/2, (edges[-1]-edges[0])/2
    Q, R, D = (q-qc)/qs, (r-rc)/rs, (delta-dc)/ds
    h, _, _ = np.histogram2d(Q, R, bins=cfg['fft_bins'], range=[[-1, 1], [-1, 1]])
    h = gaussian_filter(h, .35)
    fft = abs(np.fft.fftshift(np.fft.fft2(h, s=[cfg['fft_padding']]*2))) / max(h.sum(), 1)
    freq = np.fft.fftshift(np.fft.fftfreq(cfg['fft_padding'], d=2/cfg['fft_bins']))
    good = ((freq[:, None] >= cfg['min_q_frequency']) & (freq[:, None] <= cfg['max_q_frequency'])
            & (abs(freq[None, :]) <= 8))
    peaks = np.argwhere(good & (fft == maximum_filter(fft, size=5)))
    peaks = sorted(peaks, key=lambda z: fft[tuple(z)], reverse=True)[:cfg['phase_starts']]
    hist, bins = np.histogramdd(np.c_[Q, R, D], bins=cfg['phase_bins'], range=[[-1, 1]]*3)
    grids = np.meshgrid(*[(x[1:]+x[:-1])/2 for x in bins], indexing='ij')
    take = hist.ravel() > 0
    w = hist.ravel()[take]; w /= w.sum()
    X = _phase_design(*[g.ravel()[take] for g in grids], extended)
    full = _phase_design(Q, R, D, extended)
    penalty = np.full(X.shape[1], cfg['phase_curve_penalty']); penalty[:2] = 0
    if extended:
        penalty[6:] = cfg['phase_extra_penalty']

    def objective(b):
        E = np.exp(2j*np.pi*(X @ b)); C = np.sum(w*E)
        grad = 2*np.real(np.conj(C)*(2j*np.pi*np.sum((w*E)[:, None]*X, axis=0)))
        return -abs(C)**2 + np.sum(penalty*b*b), -grad + 2*penalty*b

    solutions = []
    for ai, bi in peaks:
        initial = np.zeros(X.shape[1]); initial[:2] = [freq[ai], freq[bi]]
        bounds = [(cfg['min_q_frequency'], cfg['max_q_frequency']), (-10, 10), (-5, 5), (-5, 5), (-3, 3), (-3, 3)]
        if extended:
            bounds += [(-3, 3)]*3
        opt = minimize(objective, initial, jac=True, method='L-BFGS-B', bounds=bounds,
                       options=dict(maxiter=cfg['phase_iterations'], ftol=1e-10))
        b = opt.x; v = full @ b
        deriv = b[0] + b[3]*R + b[4]*R*R + 2*b[5]*Q
        if extended:
            deriv += b[7]*D
        if np.quantile(deriv, .001) <= .35:
            continue
        phase = float(np.angle(np.mean(np.exp(2j*np.pi*v)))/(2*np.pi)); v -= phase
        near = abs(v-np.rint(v)) < cfg['seed_radius']
        if near.sum() < 100:
            continue
        low = int(np.floor(np.quantile(v[near], .002)+.5)); high = int(np.floor(np.quantile(v[near], .998)+.5))
        count = high-low+1
        if not 2 <= count <= min(cfg['maximum_latent_span'], cfg['max_slots']):
            continue
        solutions.append(dict(coefficients=b.tolist(), phase=phase, offset=low, count=count,
            power=float(abs(np.mean(np.exp(2j*np.pi*v)))**2), objective=float(objective(b)[0]),
            optimizer_success=bool(opt.success), min_derivative_quantile=float(np.quantile(deriv, .001))))
    if not solutions:
        raise ValueError('No monotone periodic proposal within the maximum hardware count; inspect rotation, density and frequency limits')
    solutions.sort(key=lambda x: x['objective'])
    seed = dict(extended=bool(extended), qc=float(qc), qs=float(qs), rc=float(rc), rs=float(rs),
                dcenter=float(dc), dscale=float(ds), **solutions[0], alternatives=solutions[1:])
    u = _phase_values(seed, q, r, delta); j = np.rint(u).astype(int)
    keep = (j >= 0) & (j < seed['count']) & (abs(u-j) < cfg['seed_radius'])
    return seed, j, keep


def design(j, r, delta, rscale, family, count, dcenter=0., dscale=10.):
    """Design for q=F(j,r,delta); low-order nested models, at most13 terms."""
    half = (count-1)/2
    t, R, D = np.broadcast_arrays((np.asarray(j)-half)/half, np.asarray(r)/rscale,
                                  (np.asarray(delta)-dcenter)/dscale)
    cols = [np.ones(t.shape), R, t, t*R]
    order = FAMILIES.index(family)
    if order >= 1:
        cols += [R*R, t*R*R]
    if order >= 2:
        cols += [t*t, t*t*R, t*t*R*R]
    cols += [D, t*D, R*D, t*R*D][:max(0, order-2)]
    return np.column_stack([x.ravel() for x in cols])


def _design(model, j, r, d):
    return design(j, r, d, model['rscale'], model['family'], model['count'], model['dcenter'], model['dscale'])


def surface_values(model, j, r, delta):
    """Evaluate fitted q at broadcastable latent j, rotated r, measured delta."""
    j, r, delta = np.broadcast_arrays(j, r, delta)
    return (_design(model, j, r, delta) @ np.asarray(model['coefficients'])).reshape(j.shape)


def surface_derivative(model, j, r, delta):
    return (surface_values(model, np.asarray(j)+.001, r, delta)
            - surface_values(model, np.asarray(j)-.001, r, delta))/.002


def _order_rows(model, cells, re, de, coefficients=None):
    rr, dd, jj = [], [], []
    for j in [-.5, model['count']-.5]:
        for edge in [0, 1]:
            d = de[cells[:, 1]+edge]
            for ridx in [0, 1]:
                rr.extend(re[cells[:, 0]+ridx]); dd.extend(d); jj.extend(np.full(len(d), j))
            if coefficients is not None:
                def row(rv):
                    return (_design(model, j+.001, rv, d)-_design(model, j-.001, rv, d))/.002
                zero = row(np.zeros(len(d))) @ coefficients
                minus = row(np.full(len(d), -model['rscale'])) @ coefficients
                plus = row(np.full(len(d), model['rscale'])) @ coefficients
                a, b = (plus+minus)/2-zero, (plus-minus)/2
                vertex = -b/(2*np.where(abs(a)>1e-12, a, 1))*model['rscale']
                inside = (a > 1e-12) & (vertex > re[cells[:, 0]]) & (vertex < re[cells[:, 0]+1])
                rr.extend(vertex[inside]); dd.extend(d[inside]); jj.extend(np.full(inside.sum(), j))
    rr, dd, jj = np.asarray(rr), np.asarray(dd), np.asarray(jj)
    return (_design(model, jj+.001, rr, dd)-_design(model, jj-.001, rr, dd))/.002


def _coverage_records(r, delta, j, count, edges, cfg):
    # Store every occupied bin, including low-count bins. Changing the minimum
    # later only changes this independent envelope gate, never any fit.
    width = cfg['coverage_bin_width']
    bins = np.r_[np.arange(edges[0], edges[-1], width), edges[-1]]
    result = {}
    for k in range(count):
        rows = []
        for low, high in zip(bins[:-1], bins[1:]):
            take = (j == k) & (delta >= low) & (delta < high)
            if take.any():
                a, b = np.quantile(r[take], cfg['coverage_quantiles'])
                rows.append([float((low+high)/2), float(a), float(b), int(take.sum()), float(low), float(high)])
        result[str(k)] = rows
    return result


def _domain(model, j, r, delta, minimum=None):
    cfg = model['configuration']; minimum = cfg['coverage_min_events'] if minimum is None else int(minimum)
    if minimum < 1:
        raise ValueError('coverage_min_events must be >=1')
    j, r, delta = np.broadcast_arrays(j, r, delta)
    ok = np.zeros(j.shape, bool); bridge = ok.copy()
    for k in range(model['count']):
        rows = np.asarray([x for x in model['coverage'].get(str(k), []) if x[3] >= minimum], float)
        if len(rows) < 2:
            continue
        x, lo, hi = rows[:, :3].T; ix = np.clip(np.searchsorted(x, delta), 1, len(x)-1)
        interp = (delta >= rows[ix-1, 5]) & (delta < rows[ix, 4])
        use = ((j == k) & (delta >= rows[0, 4]) & (delta < rows[-1, 5])
               & ((x[ix]-x[ix-1] <= cfg['maximum_bridge']) | ~interp))
        use &= ((r >= np.interp(delta, x, lo)-cfg['coverage_padding'])
                & (r <= np.interp(delta, x, hi)+cfg['coverage_padding']))
        ok |= use; bridge |= use & interp
    return ok, bridge


def _coordinates(model, raw, delta):
    q, r = rotate(model, raw)
    lo, hi = np.full(len(q), -.5), np.full(len(q), model['count']-.5)
    finite = ((surface_derivative(model, lo, r, delta) > 0)
              & (surface_derivative(model, hi, r, delta) > 0)
              & (surface_values(model, lo, r, delta) <= q)
              & (surface_values(model, hi, r, delta) >= q))
    for _ in range(27):
        mid = (lo+hi)/2; right = surface_values(model, mid, r, delta) < q
        lo, hi = np.where(right, mid, lo), np.where(right, hi, mid)
    return dict(q=q, r=r, u=np.where(finite, (lo+hi)/2, np.nan), finite=finite)


def _fit_surface(raw, delta, seed_j, seed_keep, frame, family, count, edges, cfg):
    q, r = rotate(frame, raw); rs = float(np.quantile(abs(r), .99))
    if rs <= 0:
        raise ValueError('Rotated r scale is zero')
    j, keep = seed_j.copy(), seed_keep.copy()
    h, re, de = np.histogram2d(r, delta, bins=[30, 40], range=[[r.min(), r.max()], [edges[0], edges[-1]]])
    cells = np.argwhere(h >= 3)
    if not len(cells):
        raise ValueError('Insufficient events to constrain ordered joint domain')
    rng = np.random.RandomState(cfg['seed']); history = []
    model = dict(schema_version=1, format='feature_only_joint_surface_v1', family=family, count=int(count),
                 **frame, rscale=rs, dcenter=float((edges[0]+edges[-1])/2), dscale=float((edges[-1]-edges[0])/2),
                 delta_edges=edges.tolist(), configuration=deepcopy(cfg))
    for step in range(cfg['iterations']):
        ix = np.flatnonzero(keep)
        if len(ix) < 100:
            raise ValueError('Fewer than100 geometry-inferred near-center events')
        if len(ix) > cfg['max_fit_events']:
            ix = rng.choice(ix, cfg['max_fit_events'], replace=False)
        X, y = _design(model, j[ix], r[ix], delta[ix]), q[ix]
        bin_id = np.clip(((delta[ix]-edges[0])/(edges[-1]-edges[0])*10).astype(int), 0, 9)
        cell = bin_id*count+j[ix]; freq = np.bincount(cell, minlength=10*count)
        w = 1/np.maximum(freq[cell], 20); w /= w.mean()
        penalty = np.full(X.shape[1], cfg['high_order_penalty']); penalty[:4] = 1e-9
        penalty[9:] = cfg['extra_penalty']
        c = np.linalg.solve((X.T*w) @ X/len(X)+np.diag(penalty), X.T @ (w*y)/len(X))
        A = _order_rows(model, cells, re, de); margin = max(.002, .10*float(np.median(abs(A @ c))))
        for _ in range(cfg['irls']):
            e = y-X @ c; sig = max(.005, 1.4826*np.median(abs(e-np.median(e))))
            ww = w/np.sqrt(1+(e/(2*sig))**2)
            H = (X.T*ww) @ X/len(X)+np.diag(penalty); b = X.T @ (ww*y)/len(X)
            c = np.linalg.solve(H, b)
            for _ in range(12):
                A = _order_rows(model, cells, re, de, c)
                if np.min(A @ c) >= margin-1e-9:
                    break
                opt = minimize(lambda z: (.5*z @ H @ z-b @ z, H @ z-b), c, jac=True, method='SLSQP',
                    constraints=[dict(type='ineq', fun=lambda z: A @ z-margin, jac=lambda z: A)],
                    options=dict(maxiter=120, ftol=1e-11))
                if not opt.success:
                    raise ValueError('Order fit failed: '+opt.message)
                c = opt.x
            else:
                raise ValueError('Ordering checks did not converge')
        model['coefficients'] = c.tolist()
        co = _coordinates(model, raw, delta); u = co['u']
        j = np.rint(np.nan_to_num(u, nan=-99)).astype(int)
        keep = np.isfinite(u) & (j >= 0) & (j < count) & (abs(u-j) < cfg['refine_radius'])
        normalized = X/np.sqrt(np.mean(X*X, axis=0))
        history.append(dict(iteration=step, fit_events=len(ix), near_center_events=int(keep.sum()),
                            condition=float(np.linalg.cond(normalized)), rms_q=float(np.sqrt(np.mean((X @ c-y)**2))),
                            minimum_order=float(np.min(A @ c))))
    model['coverage'] = _coverage_records(r[keep], delta[keep], j[keep], count, edges, cfg)
    model['diagnostics'] = dict(history=history, seed_events=int(seed_keep.sum()), minimum_order=float(np.min(A @ c)),
        reconstruction_inputs_used=False, reference_inputs_used=False, backup_used=False)
    return model


def _mixture(v, centers, bounds, width, cfg, active=None):
    centers = np.asarray(centers, float)
    active = np.arange(len(centers)) if active is None else np.asarray(active, int)
    weights = np.full(len(active)+1, 1/(len(active)+1)); sigma = width; last = -np.inf
    for _ in range(cfg['mixture_iterations']):
        normal = np.exp(-.5*((v[:, None]-centers[active])/sigma)**2)/(sigma*np.sqrt(2*np.pi))
        density = np.column_stack([normal, np.full(len(v), 1/(bounds[1]-bounds[0]))])
        joint = density*weights; total = np.maximum(joint.sum(1), 1e-300); resp = joint/total[:, None]
        weights = resp.mean(0); signal = resp[:, :-1]
        variance = np.sum(signal*(v[:, None]-centers[active])**2)/max(signal.sum(), 1e-12)
        sigma = float(np.clip(np.sqrt(variance), .6*width, 1.6*width))
        ll = float(np.log(total).sum())
        if abs(ll-last) < 1e-6*len(v):
            break
        last = ll
    return dict(active=active.tolist(), weights=weights.tolist(), sigma=sigma, bounds=list(bounds), centers=centers.tolist())


def _evaluate_mixture(mix, u):
    centers = np.asarray(mix['centers'])[np.asarray(mix['active'], int)]; sigma = mix['sigma']
    normal = np.exp(-.5*((u[:, None]-centers)/sigma)**2)/(sigma*np.sqrt(2*np.pi))
    joint = np.column_stack([normal, np.full(len(u), 1/np.diff(mix['bounds'])[0])])*mix['weights']
    total = np.maximum(joint.sum(1), 1e-300)
    return np.log(total), joint/total[:, None]


def _occupancy(model, raw, delta, cfg):
    co = _coordinates(model, raw, delta); edges = np.asarray(model['delta_edges'])
    nd = np.searchsorted(edges, delta, side='right')-1
    rng = np.random.RandomState(cfg['seed']); result = {}
    centers = np.arange(model['count'])
    for sl in range(len(edges)-1):
        ix = np.flatnonzero((nd == sl) & co['finite'])
        if len(ix) > cfg['occupancy_cap']:
            ix = rng.choice(ix, cfg['occupancy_cap'], replace=False)
        if len(ix) < cfg['occupancy_min_events']:
            result[str(sl)] = dict(accepted_indices=[], mixture_model=None, events=len(ix), candidates=[],
                                   error='Fewer than occupancy_min_events finite events')
            continue
        u, r = co['u'][ix], co['r'][ix]
        width = float(np.clip(1.4826*np.median(abs(u-np.rint(u))), .035, .3))
        bounds = (min(float(u.min()), float(centers.min()-4*width)), max(float(u.max()), float(centers.max()+4*width)))
        order = np.random.RandomState(cfg['occupancy_seed']).permutation(len(u)); tr, va = order[::2], order[1::2]
        full = _mixture(u[tr], centers, bounds, width, cfg)
        full_ll, _ = _evaluate_mixture(full, u[va]); _, resp = _evaluate_mixture(full, u)
        accepted, candidates = [], []
        for k in centers:
            absent = _mixture(u[tr], centers, bounds, width, cfg, active=np.delete(centers, k))
            absent_ll, _ = _evaluate_mixture(absent, u[va])
            gain = float(2*(full_ll-absent_ll).sum()-np.log(len(tr))); effective = float(resp[:, k].sum())
            near = (resp[:, k] > .5) & (abs(u-k) < 2*full['sigma']); coherent = 0
            if near.sum() >= 6:
                re = np.quantile(r[near], np.linspace(0, 1, 4))
                for a, b in zip(re[:-1], re[1:]):
                    coherent += int(np.sum(near & (r >= a) & (r <= b)) >= 2 and b > a)
            detected = (gain > cfg['occupancy_likelihood_gain'] and effective >= cfg['occupancy_min_effective_events']
                        and coherent >= cfg['occupancy_min_coherent_intervals'])
            candidates.append(dict(index=int(k), likelihood_gain=gain, effective_events=effective,
                                   coherent_intervals=int(coherent), detected=bool(detected)))
            if detected:
                accepted.append(int(k))
        mix = _mixture(u, centers, bounds, width, cfg, active=accepted) if accepted else None
        result[str(sl)] = dict(accepted_indices=accepted, mixture_model=mix, fitted_width=float(full['sigma']),
                              events=len(ix), candidates=candidates)
    return result


def _logdensity(model, raw, delta):
    co = _coordinates(model, raw, delta); q, r = co['q'], co['r']
    low, high = model['configuration']['density_q_bounds']; nd = np.searchsorted(model['delta_edges'], delta, side='right')-1
    out = np.full(len(q), -np.log(high-low))
    for sl in range(len(model['delta_edges'])-1):
        ix = np.flatnonzero(nd == sl); mix = model['occupancy'][str(sl)]['mixture_model']
        if not len(ix) or mix is None:
            continue
        slots = np.asarray(mix['centers'])[np.asarray(mix['active'], int)]
        j, rr, dd = np.broadcast_arrays(slots[None, :], r[ix, None], delta[ix, None])
        mu = surface_values(model, j, rr, dd); factor = surface_derivative(model, j, rr, dd)
        valid = np.isfinite(mu) & np.isfinite(factor) & (factor > 0)
        ordered = (surface_derivative(model, -.5, r[ix], delta[ix]) > 0) & (surface_derivative(model, model['count']-.5, r[ix], delta[ix]) > 0)
        valid &= ordered[:, None]; sigma = float(mix['sigma'])*factor
        valid &= np.isfinite(sigma) & (sigma > 1e-8); safe = np.where(valid, sigma, 1.)
        norm = ndtr((high-mu)/safe)-ndtr((low-mu)/safe); valid &= norm > 1e-14
        if len(slots) > 1:
            crossing = np.diff(mu, axis=1) <= 0; bad = np.zeros(mu.shape, bool)
            bad[:, :-1] |= crossing; bad[:, 1:] |= crossing; valid &= ~bad
        weights = np.asarray(mix['weights'], float); weights /= weights.sum()
        background = weights[-1]+np.sum((~valid)*weights[:-1], axis=1)
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            ll = (-.5*((q[ix, None]-mu)/safe)**2-np.log(safe*np.sqrt(2*np.pi))
                  -np.log(np.where(valid, norm, 1.))+np.log(weights[:-1]))
            ll[~valid] = -np.inf
            out[ix] = logsumexp(np.c_[ll, np.log(background/(high-low))], axis=1)
    out[(q < low) | (q > high) | ~np.isfinite(q)] = np.log(np.finfo(float).tiny)
    scores = [-float(out[nd == sl].mean()) for sl in range(len(model['delta_edges'])-1) if np.any(nd == sl)]
    return float(np.mean(scores)), out


def fit_axis(raw, delta, delta_edges, config=None, frame=None):
    """Fit one axis; return a serializable ``(model, metadata)`` pair.

    A deterministic random subset (20% by default) is used only to choose the
    phase initializer and simplest surface within0.01nat of the best balanced
    validation negative log likelihood. Those events do not fit coefficients,
    the standardization, occupancy, or coverage. A caller can reserve a further
    independent test sample or fit an independent audit model. ``frame`` permits
    the audit to use the primary fit's *fixed* coordinate system.

    The maximum hardware count rejects implausible counts; it never supplies a
    target band count or creates detections. Delta remains continuous in every
    phase/surface evaluation. Reporting intervals are used for occupancy and
    balanced scoring only. The returned index has no physical sieve name.
    """
    raw, delta = _arrays(raw, delta, minimum=500); edges = _edges(delta_edges); cfg = _configuration(config)
    if not np.isfinite(np.c_[raw, delta]).all():
        raise ValueError('Fitting inputs must all be finite')
    if np.any((delta < edges[0]) | (delta >= edges[-1])):
        raise ValueError('Fitting events must lie within half-open delta_edges range')
    order = np.random.RandomState(cfg['seed']).permutation(len(raw))
    nv = max(100, int(len(raw)*cfg['model_choice_fraction'])); va, tr = order[:nv], order[nv:]
    transform = make_frame(raw[tr], cfg['theta_deg']) if frame is None else {k: deepcopy(frame[k]) for k in ['mean', 'scale', 'theta_deg']}
    rotate(transform, raw[:1])
    initializers, candidates, seeds = {}, {}, {}
    for extended in [False, True]:
        name = 'global_delta' if extended else 'global_static'
        try:
            seed, j, keep = _discover(raw[tr], delta[tr], transform, edges, cfg, extended)
            m = _fit_surface(raw[tr], delta[tr], j, keep, transform, 'divergence', seed['count'], edges, cfg)
            m['occupancy'] = _occupancy(m, raw[tr], delta[tr], cfg)
            nll, _ = _logdensity(m, raw[va], delta[va]); power = _phase_power(seed, *rotate(transform, raw[va]), delta[va])
            initializers[name] = dict(count=seed['count'], training_power=seed['power'], validation_power=power, validation_nll=nll)
            candidates[name] = m; seeds[name] = (seed, j, keep)
        except (ValueError, np.linalg.LinAlgError) as exc:
            initializers[name] = dict(error=str(exc))
    eligible = [k for k in candidates if initializers[k]['validation_power'] >= cfg['minimum_validation_phase_power']]
    if not eligible:
        raise ValueError('No reproducible periodic initializer on independent model-choice events: '+str(initializers))
    best = min(initializers[k]['validation_nll'] for k in eligible)
    init = next(k for k in ['global_static', 'global_delta'] if k in eligible and initializers[k]['validation_nll'] <= best+cfg['model_choice_tolerance'])
    seed, j, keep = seeds[init]; models, scores, errors = {}, {}, {}
    for family in cfg['families']:
        try:
            if family == 'divergence':
                m = candidates[init]
            else:
                m = _fit_surface(raw[tr], delta[tr], j, keep, transform, family, seed['count'], edges, cfg)
                m['occupancy'] = _occupancy(m, raw[tr], delta[tr], cfg)
            scores[family] = _logdensity(m, raw[va], delta[va])[0]; models[family] = m
        except (ValueError, np.linalg.LinAlgError) as exc:
            errors[family] = str(exc)
    if not models:
        raise ValueError('No ordered family fit succeeded: '+str(errors))
    best = min(scores.values()); choice = next(f for f in cfg['families'] if f in scores and scores[f] <= best+cfg['model_choice_tolerance'])
    model = models[choice]; model['phase_seed'] = seed
    metadata = dict(training_events=int(len(tr)), model_choice_events=int(len(va)), input_events=int(len(raw)),
        training_indices=tr.tolist(), model_choice_indices=va.tolist(), initializers=initializers, initializer=init,
        selected_family=choice, validation_nll=scores, family_errors=errors, latent_count=model['count'],
        physical_names_assigned=False, reconstructed_sieve_inputs_used=False, reference_labels_used=False,
        central_cross_backup_used=False, frame_source='geometry_training_events' if frame is None else 'caller_supplied_fixed_frame')
    model['selection'] = {k: v for k, v in metadata.items() if k not in ['training_indices', 'model_choice_indices']}
    return model, metadata


def predict_axis(model, raw, delta, coverage_min_events=None):
    """Predict latent memberships with separate inverse, occupancy, coverage gates.

    ``labels`` is -1 for rejected events. ``proposed_label`` retains the nearest
    occupied ridge where one exists. Nearest ridges are measured in q at the
    event's own (r,delta), as in the pilot. Posterior is diagnostic only: this
    routine does not introduce an unrequested Gaussian/core cleanup cut.
    ``coverage_min_events`` changes only the independent coverage count rule.
    """
    raw, delta = _arrays(raw, delta); n = len(delta)
    good = np.isfinite(np.c_[raw, delta]).all(1); safe_raw = np.where(np.isfinite(raw), raw, 0.)
    safe_delta = np.where(np.isfinite(delta), delta, 0.)
    co = _coordinates(model, safe_raw, safe_delta); co['finite'] &= good
    for key in ['q', 'r', 'u']:
        co[key][~good] = np.nan
    edges = np.asarray(model['delta_edges']); nd = np.searchsorted(edges, safe_delta, side='right')-1
    in_delta = good & (safe_delta >= edges[0]) & (safe_delta < edges[-1]); labels = np.full(n, -1, int)
    occupied = np.zeros(n, bool); posterior = np.zeros(n)
    for sl in range(len(edges)-1):
        ix = np.flatnonzero((nd == sl) & co['finite'] & in_delta)
        occ = model['occupancy'][str(sl)]; slots = np.asarray(occ['accepted_indices'], int)
        if not len(ix) or not len(slots):
            continue
        j, rr, dd = np.broadcast_arrays(slots[None, :], co['r'][ix, None], safe_delta[ix, None])
        mu = surface_values(model, j, rr, dd)
        labels[ix] = slots[np.argmin(abs(co['q'][ix, None]-mu), axis=1)]; occupied[ix] = True
        if occ['mixture_model'] is not None:
            _, resp = _evaluate_mixture(occ['mixture_model'], co['u'][ix]); posterior[ix] = resp[:, :-1].max(1)
    coverage, interpolation = _domain(model, labels, co['r'], safe_delta, coverage_min_events)
    accepted = co['finite'] & in_delta & occupied & coverage
    reason = np.full(n, 'accepted', dtype='<U16')
    reason[~coverage] = 'coverage'; reason[~occupied] = 'occupancy'
    reason[~co['finite']] = 'inverse_range'; reason[~in_delta] = 'delta_range'; reason[~good] = 'nonfinite'
    return dict(**co, labels=np.where(accepted, labels, -1), proposed_label=labels,
                accepted=accepted, reason=reason, coverage=coverage, occupancy=occupied,
                posterior=posterior, interpolation=interpolation & accepted, delta_bin=nd)
