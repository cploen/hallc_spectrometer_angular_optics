# HMS 6.667 — equal15 foil, delta and hole balancing

Conversation recap, 15 September 2026. Exact specification:
[CORE_BALANCE_DESIGN.md](../../docs/CORE_BALANCE_DESIGN.md).
Implementation checks: [CORE_BALANCE_VERIFICATION.md](../../docs/CORE_BALANCE_VERIFICATION.md).
Campaign result: [equal15/HMS_6p667_CORE_BALANCE.md](equal15/HMS_6p667_CORE_BALANCE.md).

## Final rule and why

The user chose a simple, tunable scheme over automatic adaptation. Use accepted
development cores from frozen `min10`, combining its old fit and surplus pools.
Preserve protected membership and quality, and exclude blocked labels from
eligibility. Balance the five **physical positions** −8, −3, 0, +3, +8 cm across
all seven rungroups, not local foil index 0/1 or each run separately.

1. Within each rungroup/physical foil/delta, cap each hole at
   `floor(1.5 × P25(original positive eligible hole counts))`.
2. Within each physical foil, cap each delta at
   `floor(1.5 × P25(original positive eligible delta totals))`.
   Its capacity cannot exceed the sum of its already capped holes.
3. Sum delta capacities for each foil. Give every foil the same budget:
   `min(weakest foil capacity, floor(200000 / 5))`.
4. Fill delta budgets equally up to capacities; redistribute unused shares.
   Fill the distinct setting/hole leaves in each delta the same way. Resolve
   integer ties and select IDs deterministically; unselected eligible events
   become surplus.

P25 is a linearly interpolated percentile, **not** the minimum or the mean of
the lowest quarter. Compute both reference percentiles once, before thinning.
Strict equality is calculated after feasible capacities; there are not three
successive random cuts. The former 80% training fraction and 400-event cap are
inactive. There is **no ten-event floor, automatic cap relaxation, square-root
weighting, duplication, or borrowing from protected/noncore events**.
`qa_min=10` reports sparse allocations only. Zero feasible campaign budget stops
export instead of silently dropping a required foil.

## Campaign outcome

| Quantity | Result |
|---|---:|
| Eligible development cores | 310,513 |
| Training | 68,535 = 5 × 13,707 |
| Surplus cores | 241,978 |
| Protected holdout, unchanged | 101,258 |
| Limiting foil | +3 cm after hole caps |
| Additional loss from delta caps | 0 at these settings |

The 200,000 cap did not bind. Of the available events, hole ceilings removed
131,987 from feasible capacity; foil equality left another 109,991 feasible
events unused. These are allocation losses to training, retained as surplus.
The completed campaign report has `actual_selected=quota` for every foil;
the earlier verification note describes the preview/fixture stage.

Equal foil totals do not guarantee equal counts at a particular delta across
foils, equal setting contributions, or equal design-matrix leverage. Delta 0
training totals are 2,741 / 1,495 / 2,741 / 1,430 / 2,742 in foil order.
Slices have widths 2, 3, 5, 5, 5 percentage points; this policy balances counts,
not counts per unit delta.

## Evidence, reproduction and inputs

- `equal15/pages/delta0.md` through `delta4.md` compare all five foil locations;
  numbered log-circle maps, event clouds and allocation tables retain the
  setting-level views. These clouds contain available/training/surplus **cores**.
- The historical audit in
  `06b_svd_fit/diagnostics/balance/HMS_6p667_SVD_SIEVE_BALANCE.md` instead describes
  GMM events and actual historical SVD admission. That solver selected in input
  order, not a random ROOT draw, and capped rungroup/foil totals rather than
  enforcing physical-foil equality. These two balance schemes are not identical.
- Code: `core_balance.py`, `reallocate_core.py`, `balance_diagnostics.py`.
  Policy: `config/core_balance.json`. Commands are in the campaign analysis index.
- Tests: `test_core_balance.py`, `test_frozen_reallocation.py`, and
  `test_preallocated_svd.py` under `tests/`; the verification note records the
  completed fixture checks and their limitations.
- GitHub contains `equal15` review products, but **not its event-level manifest
  or ROOT masks** at this audit. Keep the complete original tag on ifarm.
  `min10` now has seven tracked masks; do not rerun density merely to obtain
  inputs already frozen there. Recreating a tag is not a substitute for the
  original manifest when downstream runs require its exact checksum.
