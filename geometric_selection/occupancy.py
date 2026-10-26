"""Conditional occupancy evidence for a fitted ordered X geometry.

A possible slot is accepted only with held-out density evidence and longitudinal
support. Mixture foreground flags remain diagnostics; prediction labels are not
silently removed using those flags.
"""
import numpy as np
from .x_geometry import global_coordinates, _raw_array

CONFIG = dict(mixture_iterations=120, heldout_delta_threshold=10.,
              min_effective_events=8., coherence_intervals=3,
              min_coherent_intervals=2, seed=32191)


def mixture(v, centers, bounds, initial_width, active=None):
    centers = np.asarray(centers)
    if active is None:
        active = np.arange(len(centers))
    active = np.asarray(active, dtype=int)
    selected = centers[active]
    weights = np.full(len(active) + 1, 1. / (len(active) + 1))
    sigma = initial_width
    background = 1. / (bounds[1] - bounds[0])
    last = -np.inf
    for iteration in range(CONFIG['mixture_iterations']):
        normal = np.exp(-.5 * ((v[:, None] - selected) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
        density = np.column_stack([normal, np.full(len(v), background)])
        joint = density * weights
        total = np.maximum(joint.sum(1), 1e-300)
        responsibility = joint / total[:, None]
        weights = responsibility.mean(0)
        signal = responsibility[:, :-1]
        variance = np.sum(signal * (v[:, None] - selected) ** 2) / max(signal.sum(), 1e-12)
        sigma = float(np.clip(np.sqrt(variance), .6 * initial_width, 1.6 * initial_width))
        likelihood = float(np.log(total).sum())
        if abs(likelihood - last) < 1e-6 * len(v):
            break
        last = likelihood
    return dict(active=active.tolist(), weights=weights.tolist(), sigma=sigma,
                bounds=list(bounds), centers=centers.tolist())


def evaluate(model, v):
    active = np.asarray(model['active'])
    centers = np.asarray(model['centers'])[active]
    sigma = model['sigma']
    normal = np.exp(-.5 * ((v[:, None] - centers) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    background = np.full((len(v), 1), 1. / np.diff(model['bounds'])[0])
    joint = np.column_stack([normal, background]) * model['weights']
    total = np.maximum(joint.sum(1), 1e-300)
    return np.log(total), joint / total[:, None]


def detect(prior, raw, seed=32191, coordinates=None):
    if coordinates is None:
        raise ValueError('Occupancy detection requires explicit fitted coordinates')
    coord = coordinates
    v, p = coord['v'], coord['p']
    valid = np.isfinite(v) & (coord['expansion'] > 0)
    vv, pp = v[valid], p[valid]
    centers = np.asarray(prior['centers'])
    width = prior['initial_width']
    if len(vv) < 100:
        raise ValueError('Insufficient events inside monotone global geometry (need >=100)')
    if len(centers) < 1 or not np.isfinite(centers).all() or np.any(np.diff(centers) <= 0):
        raise ValueError('Occupancy centers must be finite and ordered')
    if not np.isfinite(width) or width <= 0:
        raise ValueError('The initial occupancy width must be positive and finite')
    lo = min(float(vv.min()), float(centers.min() - 4 * width))
    hi = max(float(vv.max()), float(centers.max() + 4 * width))
    order = np.random.RandomState(seed).permutation(len(vv))
    train, test = order[::2], order[1::2]
    full = mixture(vv[train], centers, (lo, hi), width)
    full_ll, _ = evaluate(full, vv[test])
    _, resp_all = evaluate(full, vv)
    candidates, accepted = [], []
    for k, (j, center) in enumerate(zip(prior['indices'], centers)):
        absent = mixture(vv[train], centers, (lo, hi), width,
                         active=np.delete(np.arange(len(centers)), k))
        absent_ll, _ = evaluate(absent, vv[test])
        delta = float(2 * (full_ll - absent_ll).sum() - np.log(len(train)))
        effective = float(resp_all[:, k].sum())
        near = (resp_all[:, k] > .5) & (abs(vv - center) < 2 * full['sigma'])
        coherent = 0
        if near.sum() >= 6:
            edges = np.quantile(pp[near], np.linspace(0, 1, 4))
            for a, b in zip(edges[:-1], edges[1:]):
                coherent += int(np.sum(near & (pp >= a) & (pp <= b)) >= 2 and b > a)
        accept = (delta > CONFIG['heldout_delta_threshold'] and
                  effective >= CONFIG['min_effective_events'] and
                  coherent >= CONFIG['min_coherent_intervals'])
        candidates.append(dict(index=int(j), center=float(center),
            heldout_delta=delta, effective_events=effective,
            coherent_intervals=int(coherent), detected=bool(accept)))
        if accept:
            accepted.append(k)
    partition = np.full(len(v), -1, dtype=int)
    support = np.zeros(len(v), dtype=bool)
    signal_model = None
    if accepted:
        accepted = np.asarray(accepted)
        signal_model = mixture(vv, centers, (lo, hi), width, active=accepted)
        _, resp = evaluate(signal_model, vv)
        partition[valid] = np.abs(vv[:, None] - centers[accepted]).argmin(1)
        support[valid] = resp[:, :-1].sum(1) > .5
    return dict(prior=prior, candidates=candidates,
        accepted_indices=[int(prior['indices'][k]) for k in accepted],
        partition=partition, supported=support, v=v, p=p,
        q=coord['q'], finite_geometry=valid,
        mixture_model=signal_model,
        fitted_width=float(full['sigma']),
        note='Occupancy split is held out conditional on supplied fitted geometry.')


def detect_x(geometry, target, seed=32191):
    """Evaluate occupancy using target events, independent of the geometry window.

    The held-out split applies to the occupancy test conditional on the supplied
    geometry. It is not an independent-event validation of the geometry fit.
    """
    target = _raw_array(target)
    coordinates = global_coordinates(geometry, target)
    finite = coordinates['finite_geometry']
    u = coordinates['u']
    if finite.sum() < 100:
        raise ValueError('Insufficient events inside monotone global geometry (need >=100)')
    residual = abs(u[finite] - np.rint(u[finite]))
    width = float(np.clip(1.4826 * np.median(residual), .035, .3))
    indices = list(geometry['possible_indices'])
    prior = dict(fan=None, indices=indices, centers=indices,
                 initial_width=width, geometry='single_global_curved_family')
    converted = dict(q=coordinates['q'], p=coordinates['p'], v=u,
                     expansion=np.where(finite, 1., -1.))
    result = detect(prior, target, seed=seed, coordinates=converted)
    result['global_geometry'] = geometry
    return result


def predict_x(geometry, occupancy, target):
    """Assign latent slot indices; -1 denotes points outside fitted geometry.

    These indices still require physical sieve labeling by the next campaign
    stage. Foreground diagnostics do not remove an otherwise assigned event.
    """
    target = _raw_array(target)
    coord = global_coordinates(geometry, target)
    finite = coord['finite_geometry']
    labels = np.full(len(target), -1, dtype=int)
    foreground = np.zeros(len(target), dtype=bool)
    accepted = np.asarray(occupancy['accepted_indices'], dtype=int)
    if len(accepted):
        u = coord['u'][finite]
        labels[finite] = accepted[np.abs(u[:, None] - accepted).argmin(1)]
        _, posterior = evaluate(occupancy['mixture_model'], u)
        foreground[finite] = posterior[:, :-1].sum(1) > .5
    return dict(labels=labels, foreground=foreground, **coord)
