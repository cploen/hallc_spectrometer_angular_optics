"""Admission checks for frozen, unit-weight core/inner-shoulder allocations."""
import numpy as np


def mixed_fraction(manifest):
    schema = manifest.get('schema')
    if schema not in ('inner10_allocation_v1', 'inner20_allocation_v1'):
        if schema != 'core_balance_v1':
            raise ValueError('Unsupported training allocation schema')
        return None
    fraction = .1 if schema == 'inner10_allocation_v1' else .2
    policy = manifest.get('policy', {})
    required = dict(core_fraction=1-fraction, inner_shoulder_fraction=fraction,
                    event_weight=1, inner_quality=1, q_low_inclusive=.4, q_high_exclusive=.5)
    if any(policy.get(k) != v for k, v in required.items()):
        raise ValueError('Mixed allocation policy disagrees with schema')
    return fraction


def validate_training(a, selected, mixed=False):
    quality, core = a['quality'][selected], a['core_keep'][selected]
    valid = (quality == 2) & (core == 1)
    if mixed:
        score = a['core_score'][selected]
        valid |= ((quality == 1) & (core == 0) & np.isfinite(score)
                  & (score >= .4) & (score < .5))
    excluded = a['balance_excluded'][selected] if 'balance_excluded' in a else a['excluded'][selected]
    if np.any(~valid | (excluded != 0)):
        raise ValueError('Training contains an ineligible or blocked event')


def validate_totals(foil, quality, manifest):
    fraction = mixed_fraction(manifest)
    if fraction is None:
        return
    n = len(quality)
    totals = manifest['totals']
    if (n != totals['fit'] or np.sum(quality == 2) != totals['core']
            or np.sum(quality == 1) != totals['inner_shoulder']):
        raise ValueError('Mixed training totals disagree with allocation')
    for z in np.unique(foil):
        q = quality[foil == z]
        if len(q) != manifest['foil_budget'] or int(np.sum(q == 1)) * 10 != len(q) * round(fraction*10):
            raise ValueError('Mixed allocation violates the per-foil population budget')
