#!/usr/bin/env python3
"""Fixed-cap quota arithmetic shared by preview and frozen event reallocation."""
import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction

LEAF = ('rungroup', 'zfoil', 'delta_low', 'delta_high', 'xscol', 'yscol')


def validate(cfg):
    b = cfg['balance']
    if set(b) != {'foil', 'delta', 'hole'} or b['foil'] != 'equal':
        raise ValueError('balance requires foil=equal and fixed delta/hole factors')
    for k in ('delta', 'hole'):
        if isinstance(b[k], bool) or not isinstance(b[k], (int, float)) or not math.isfinite(b[k]) or b[k] <= 0:
            raise ValueError(f'balance.{k} must be positive and finite')
    for k in ('fit_max', 'qa_min'):
        if type(cfg[k]) is not int or cfg[k] <= 0:
            raise ValueError(f'{k} must be a positive integer')
    if type(cfg['seed']) is not int:
        raise ValueError('seed must be an integer')


def p25(counts):
    a = sorted(int(n) for n in counts if n > 0)
    if not a:
        return Fraction(0)
    i, rem = divmod(len(a) - 1, 4)
    return Fraction(a[i]) + Fraction(rem, 4) * (a[min(i + 1, len(a)-1)] - a[i])


def tie(seed, key):
    canonical = json.dumps(key, separators=(',', ':'), ensure_ascii=True)
    return hashlib.sha256(f'{seed}:quota:{canonical}'.encode()).digest(), canonical


def fair(capacities, budget, keys, seed):
    """Integer water filling in O(n log(max capacity)); keyed residual ties."""
    caps = [int(v) for v in capacities]
    if any(v < 0 for v in caps) or budget < 0 or budget > sum(caps):
        raise ValueError('Infeasible equal-share budget')
    if len(keys) != len(caps) or len({json.dumps(k) for k in keys}) != len(keys):
        raise ValueError('Quota keys must be distinct')
    lo, hi = 0, max(caps, default=0)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if sum(min(c, mid) for c in caps) <= budget:
            lo = mid
        else:
            hi = mid - 1
    out = [min(c, lo) for c in caps]
    active = sorted((i for i, c in enumerate(caps) if c > lo), key=lambda i: tie(seed, keys[i]))
    for i in active[:budget - sum(out)]:
        out[i] += 1
    assert sum(out) == budget and all(0 <= q <= c for q, c in zip(out, caps))
    return out


def allocate_counts(rows, cfg, required_foils, intervals):
    validate(cfg)
    foils = sorted(required_foils)
    if not foils or len(set(foils)) != len(foils):
        raise ValueError('Required physical foils must be explicit and distinct')
    intervals = sorted(tuple(d) for d in intervals)
    if not intervals or len(set(intervals)) != len(intervals) or any(not (math.isfinite(a) and math.isfinite(b) and a < b) for a,b in intervals):
        raise ValueError('Invalid delta intervals')
    if any(a[1] > b[0] for a,b in zip(intervals, intervals[1:])):
        raise ValueError('Overlapping delta intervals')
    leaves = sorted((dict(r) for r in rows), key=lambda r: tuple(r[k] for k in LEAF))
    seen, slices = set(), defaultdict(list)
    for r in leaves:
        key = tuple(r[k] for k in LEAF)
        if key in seen:
            raise ValueError(f'Duplicate allocation leaf: {key}')
        seen.add(key)
        if r['zfoil'] not in foils or (r['delta_low'], r['delta_high']) not in intervals:
            raise ValueError(f'Incompatible foil/delta: {key}')
        if type(r['available']) is not int or r['available'] < 0:
            raise ValueError('Availability must be a nonnegative integer')
        slices[key[:4]].append(r)
    for cell in slices.values():
        ref = p25(r['available'] for r in cell)
        cap = math.floor(Fraction(str(cfg['balance']['hole'])) * ref)
        for r in cell:
            r.update(q_hole=float(ref), hole_cap=cap, capacity=min(r['available'], cap), hole_bound=r['available']>cap, population='blocked' if r.get('blocked') else ('open_populated' if r['available']>0 else 'observed_no_eligible_core'))
    deltas, foil_rows = [], []
    for f in foils:
        cells = [[r for r in leaves if r['zfoil'] == f and (r['delta_low'], r['delta_high']) == d] for d in intervals]
        avail = [sum(r['available'] for r in cell) for cell in cells]
        ref = p25(avail)
        cap = math.floor(Fraction(str(cfg['balance']['delta'])) * ref)
        for (low, high), cell, n in zip(intervals, cells, avail):
            holes = sum(r['capacity'] for r in cell)
            deltas.append(dict(zfoil=f, delta_low=low, delta_high=high, width=high-low,
                               available=n, q_delta=float(ref), delta_cap=cap, hole_capacity=holes,
                               capacity=min(holes, cap), hole_loss=n-holes, delta_loss=max(0, holes-cap), delta_bound=holes>cap))
        foil_rows.append(dict(zfoil=f, available=sum(avail), capacity=sum(r['capacity'] for r in deltas if r['zfoil']==f)))
    budget = min(min(r['capacity'] for r in foil_rows), cfg['fit_max']//len(foils))
    for fr in foil_rows:
        ds = [r for r in deltas if r['zfoil'] == fr['zfoil']]
        quotas = fair([r['capacity'] for r in ds], budget,
                      [(fr['zfoil'], r['delta_low'], r['delta_high']) for r in ds], cfg['seed'])
        for d, quota in zip(ds, quotas):
            d.update(quota=quota, surplus=d['available']-quota, budget_loss=d['capacity']-quota)
            cell = [r for r in leaves if r['zfoil']==d['zfoil'] and (r['delta_low'],r['delta_high'])==(d['delta_low'],d['delta_high'])]
            qs = fair([r['capacity'] for r in cell], quota, [tuple(r[k] for k in LEAF) for r in cell], cfg['seed'])
            for r, q in zip(cell, qs):
                r.update(quota=q, surplus=r['available']-q, qa_low=r['available']>0 and q<cfg['qa_min'],
                         hole_loss=r['available']-r['capacity'], parent_loss=r['capacity']-q)
        fr.update(quota=budget, surplus=fr['available']-budget,
                  hole_loss=sum(d['hole_loss'] for d in ds), delta_loss=sum(d['delta_loss'] for d in ds),
                  budget_loss=fr['capacity']-budget, limiting=fr['capacity']==min(r['capacity'] for r in foil_rows))
    # Add disjoint development categories plus protected counts at every parent level.
    for parents in (deltas, foil_rows):
        for parent in parents:
            children=[r for r in leaves if r['zfoil']==parent['zfoil'] and ('delta_low' not in parent or (r['delta_low'],r['delta_high'])==(parent['delta_low'],parent['delta_high']))]
            for k in ('original_core','excluded_core','total','holdout','holdout_core','unsupported','noncore_development','geometry_excluded_development'):
                parent[k]=sum(r.get(k,0) for r in children)
    zero = [r['zfoil'] for r in foil_rows if r['capacity']==0]
    return dict(leaves=leaves, deltas=deltas, foils=foil_rows, foil_budget=budget,
                training_total=len(foils)*budget, zero_capacity_foils=zero,
                status='zero_budget' if budget==0 else 'quotas_ready',
                budget_constraint='zero_feasible_foil' if zero else ('campaign_ceiling' if cfg['fit_max']//len(foils)<min(r['capacity'] for r in foil_rows) else 'weakest_capped_foil'),
                inactive_legacy_controls=['fit_fraction', 'fit_cap', 'flat_allocation'])
