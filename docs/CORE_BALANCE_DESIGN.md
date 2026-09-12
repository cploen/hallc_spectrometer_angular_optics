# Strict foil balancing with fixed delta and hole caps

Implementation handoff — revised September 11, 2026. Implementation was subsequently authorized in the handoff conversation. This
replaces the earlier adaptive-cap proposal in full. Do not implement automatic
cap relaxation, a retention target, square-root weights, or a ten-event floor.
Do not run a new production SVD fit or replace a matrix as part of this handoff.

## Agreed objective

Use a clearly explainable selection policy that can be studied and optimized
later. Assign exactly equal training-event totals to the physical target foils,
while limiting dense delta slices and sieve-hole populations. Train on accepted
core events only. Preserve the protected holdout and retain unused eligible
core events as surplus.

For HMS 6.667, balance the five physical positions **-8, -3, 0, +3, +8 cm** across
the campaign. Equal event counts give each physical position equal total count
contribution to an ordinary unweighted least-squares objective. They do not make
its design-matrix rows equally influential: focal-plane coordinates and polynomial
terms still affect leverage and conditioning.

This follows the earlier foil-balancing motivation, but does not exactly recreate
its effective weighting. The historical solver capped events per *rungroup/foil*,
not per physical foil aggregated across the campaign. Repeatedly measuring a
physical foil could therefore give it more total events. Preserve the historical
fit and include its actual totals as a diagnostic comparison.

## Three rules, without automatic exceptions

1. **Hole cap:** within each rungroup / physical foil / delta slice, limit every
   open hole to 1.5 times the 25th percentile of available core-event counts.
2. **Delta cap:** within each physical foil, limit every delta slice to 1.5 times
   the 25th percentile of its available delta-slice core-event counts, combining
   that foil's rungroups.
3. **Strict foil equality:** after accounting for those limits, give every physical
   foil the same feasible number of training events. The weakest feasible foil
   determines that number, subject to the campaign maximum.

The two multipliers remain configurable, but are fixed at 1.5 for the first
review. They are hard ceilings for that run. Never increase them automatically
or pull holdout/noncore events into training to reach a target.

P25 means the linearly interpolated 25th percentile of positive counts. It is
not the average of the lowest quarter, nor the minimum count. Use only eligible
open populations in the reference calculation. A nearly empty population does
not directly dictate the cap for every hole.

## Input definitions and identities

- Start from one frozen core-sample tag, initially HMS 6.667 `min10`. Reuse its
  density regions and quality assignments. Do not rerun density estimation to
  compare allocation policies.
- Eligible events have accepted core quality (`quality == 2`), are outside
  `sample == 2` protected holdout, and are not geometrically excluded. Combine
  the old fit and surplus pools before assigning new quotas.
- Preserve event identity `(rungroup, entry)`. Current candidate run fields use
  the optics ID and merged replay entry; do not invent original-run IDs.
- Map local foil indices to physical ztar using saved metadata. Never treat local
  foil index 0 as a unique physical position across rungroups.
- Identify delta by numerical interval boundaries, not just `ndel`. Initially
  require aligned boundaries across the pool; reject incompatible definitions
  instead of silently merging or rebinning frozen selections.
- Keep each `(rungroup, physical foil, delta interval, xscol, yscol)` as a distinct
  allocation leaf. Same-numbered holes at different settings are not identical
  measurements, even when aggregate foil diagnostics combine their counts.
- Pool only explicitly compatible settings within one spectrometer/campaign.
  Snapshot the grouping. Default to the existing rungroup identities; do not
  infer equivalence from similar names or nearly equal central angles.
- Use the campaign geometry mask. HMS blocked label pairs `[2,3]` and `[5,5]`
  are never eligible. Retain excluded events and reasons in the audit; leave
  protected holdout membership and frozen quality classifications unchanged.
- Distinguish open populated labels, labels with no accepted core events,
  confirmed blocked positions, and absent labels. An absent label does not
  establish whether a physical hole should have been illuminated. The current
  mask is not a full illumination inventory.

## Exact budget calculation

Calculate quotas from counts before choosing any event IDs. There must not be
three successive random cuts that independently discard events.

Let `N[r,f,d,h]` be eligible core events at a leaf, and let `H` and `D` be the
hole and delta multipliers, initially 1.5. Use integer quotas throughout.

### A. Compute references once from original eligible counts

For each `(r,f,d)`, compute:

`Q_hole[r,f,d] = P25(positive N[r,f,d,h])`.

For each physical foil and delta, compute the original available total:

`A[f,d] = sum over r,h of N[r,f,d,h]`.

Then for each physical foil:

`Q_delta[f] = P25(positive A[f,d])`.

Do not recompute the delta reference after hole trimming, or recompute either
reference after a foil budget has been assigned. This avoids a moving reference
that shrinks repeatedly as earlier cuts are applied.

### B. Compute feasible capacities from holes upward

Hole capacity:

`C[r,f,d,h] = min(N[r,f,d,h], floor(H * Q_hole[r,f,d]))`.

Delta capacity:

`C_delta[f,d] = min(sum over r,h of C[r,f,d,h], floor(D * Q_delta[f]))`.

Foil capacity:

`C_foil[f] = sum over d of C_delta[f,d]`.

With `F` required physical foils and campaign ceiling `fit_max`, the equal foil
budget is:

`B_foil = min(min over f of C_foil[f], floor(fit_max / F))`.

The exact campaign training size is `F * B_foil`. Leaving fewer than `F` slots
unused below the campaign ceiling is preferable to breaking foil equality.

**The weakest foil here is the weakest after the delta and hole caps.** Its
feasible count can be smaller than its original available count. Show both in
reports so the source of any statistical loss is explicit.

If a required foil has no eligible events, stop and report it. Do not silently
remove it from the equality requirement. Empty delta slices and holes have zero
capacity and do not enter P25 calculations. Explicitly report zero/very small
sample outcomes; do not bypass a cap to make them look acceptable.

### C. Fill equally within the capacities

For each foil, distribute its `B_foil` budget across delta slices by equal-share
filling: give each slice the same allowance until a slice reaches its capacity,
then redistribute the unused allowance equally among the remaining slices.

Within each `(f,d)`, distribute that slice's assigned budget in the same way
across its distinct `(r,h)` leaves, respecting each hole capacity. This is hole
balancing across the contributing settings, not an additional independent setting
quota. Report contributions by rungroup so a setting cannot vanish unnoticed.

A simple example of equal-share filling: capacities `[100, 300, 600]` and budget
700 produce `[100, 300, 300]`. The smallest population is used in full; it does
not force the other two to stop at 100. Budget 450 gives `[100, 175, 175]`.

Use the existing `fair_counts` concept where appropriate, but resolve integer
remainders through stable seed/key tie-breaking rather than input table order.
Use grouped arithmetic rather than a slow one-event-at-a-time loop. Ensure exact
parent/child sums and `0 <= quota <= capacity <= available` everywhere.

There is no special ten-event reserve and no square-root weighting. Positive
leaves receiving zero or very few events must be explicitly highlighted in QA;
never misreport them as absent or rejected by the core selector. If coverage is
unacceptable, change the named policy parameters in a later reviewed trial,
not through a hidden fallback.

### D. Select events once

Within each leaf, choose the assigned quota with the existing reproducible
`ranks(seed, rungroup, entry, "fit")` ordering. Select among all eligible core
events without favoring higher density scores or earlier replay entries.

Save selected event IDs. Eligible events not selected become surplus. Preserve
noncore/unsupported data and the protected holdout. No duplication, oversampling,
or additional event weights are introduced.

## Short configuration and compatibility

Suggested new mode, in JSON rather than long shell arguments:

```json
{
  "balance": {
    "foil": "equal",
    "delta": 1.5,
    "hole": 1.5
  }
}
```

Retain existing `seed`, `holdout`, and `fit_max` (initially 200000). This mode
replaces `fit_fraction=0.8`, `fit_cap=400`, and the old flat allocation. Explicitly
record that those legacy controls are inactive. If `balance` is absent, preserve
legacy behavior. Validate positive finite factors and a positive integer budget.

Use existing short campaign/tag command conventions. Reallocation of a frozen
tag writes a new tag and snapshots the parent manifest, configuration, geometry
mask, metadata, hashes, and seed. Never overwrite `min10` or historical outputs.

The prototype adapter is currently HMS-limited. Keep spectrometer dispatch
intact, but do not claim end-to-end SHMS export/solver support without checking
that adapter and geometry separately.

## Delta interpretation and cross-foil comparisons

The initial policy balances *counts per configured delta slice*, not events per
unit delta. HMS slice widths are 2, 3, 5, 5, and 5 percentage points. Display the
interval bounds and widths so this choice is explicit. Do not add width weighting
to the first implementation; it can be studied later.

Strict equal foil totals plus within-foil delta caps **do not guarantee identical
counts for a particular delta slice across all foils**. A foil with a sparse
slice supplies fewer events there and can use more in its other slices. This
is deliberate: the first scheme imposes exact foil totals and bounded, partially
equalized delta/hole populations, without requiring the complete foil-by-delta
table to match cell for cell.

The requested cross-foil diagnostics must make this visible. If slice 0 remains
unacceptably imbalanced across physical foils, report the actual counts and
limiting capacities. Do not quietly add a fourth rule enforcing foil-by-delta
equality. That would require a separate choice about further statistics loss.

## Required diagnostics: all slices, all physical foils

Produce the complete campaign diagnostics, not just rg01 or selected examples.
For HMS 6.667 this means **five delta slices at all five physical foil locations**,
plus the underlying rungroup/local-foil breakdown wherever data exist.

### 1. Numerical budget tables

Save original eligibility, P25 references, fixed caps, effective capacities,
assigned quotas, actual selected counts, surplus, and exclusions at every level.
Provide:

- one foil table showing equal final totals;
- a 5-foil by 5-delta matrix of selected counts, with matching available and
  selected/available fraction tables;
- totals per rungroup/physical-foil/delta;
- per-hole tables, including positive available populations receiving zero or
  very few training events;
- the historical SVD's actual foil/delta totals alongside the new selection,
  clearly marked as different input populations (historical GMM versus new cores).

### 2. Log-colored circles with numbers

For every physical foil/delta and every contributing rungroup/foil/delta, create
nominal sieve-plane maps with counts printed inside circles:

**Available development core events | selected training events | surplus core events.**

Use one shared log color scale within each comparison set. Mark blocked positions
with a cross, zero counts with hollow zeros, and absent labels with gray dots.
Show totals, positive-population coverage, selected fraction, and nonzero count
range. Do not label event counts simply as “cores.” Add a companion map of the
per-hole selected fraction on a shared 0–1 linear scale when useful for exposing
strong thinning.

Aggregate physical-foil maps sum same-label counts across settings for the
requested overview, and state that aggregation. Preserve separate setting maps;
aggregate counts alone do not establish identical focal-plane coverage.

### 3. Event clouds

For every physical foil/delta and every contributing rungroup/foil/delta, create
matching reconstructed `(ysieve, xsieve)` cloud panels:

**Available development core events | selected training events | surplus core events.**

Use fixed axes, binning, and shared log counts per bin, as in the accepted rg01
cloud figures. Include event totals and counts outside the displayed frame.
No red density boundaries. Label combined-setting clouds clearly: differences
between settings may broaden the aggregate appearance, so include the individual
setting views before interpreting that broadening as poor core selection.

These new core-allocation clouds are distinct from the already completed
historical sequence **original joint labels → GMM → historical SVD**. Preserve
that sequence and link it as context; do not relabel its GMM events as core events.
Protected holdout is not part of these allocation plots.

### 4. Browse by delta slice across the target

Create a descriptively named index such as `HMS_6p667_CORE_BALANCE.md` and five
cross-foil comparison pages, one per delta slice. Within each page, order the
physical foils **-8, -3, 0, +3, +8 cm**. Include:

- a compact count table for that delta at the five foils;
- a contact sheet of the five selected-training circle maps, on the same scale;
- a contact sheet of the five selected-training event clouds, on the same scale;
- links to each foil's full available/training/surplus comparison and its settings.

The user must be able to inspect delta slice 0 across all five foil positions
without opening unrelated rungroup pages. Missing populations must appear as
explicit empty panels/entries, not disappear from the page. Use large readable
numbers and legends near the relevant panels; avoid tiny five-row figure text.
Export full-resolution individual PNGs as well as the overview sheets.

## Ensure the solver uses the supplied sample

Add an explicit preallocated-input mode to the current SVD workflow. Require
an allocation/build manifest, then admit every valid supplied training event
exactly once. Bypass the old Y-column, per-rungroup/foil, and global sequential
truncation rules in this mode. Check resource limits before fitting and fail
rather than silently truncating an oversized supplied sample.

Do not introduce extra foil weights or another balancing step. Retain normal
validity checks, report any rejected IDs, and treat a mismatch between intended
and admitted membership as a failed handoff. Write the actual solver-used event
IDs and per-leaf counts. Preserve legacy GMM mode unchanged.

No new production fit or matrix replacement is requested here. Verify this path
with a compact integration fixture; return the diagnostics for review first.

## Implementation sequence and reproducibility

1. Implement a counts-only preview from frozen `allocation.tsv`/quality summaries
   and the mask. It must use exactly the same quota logic as event-level allocation.
2. Implement frozen-tag reallocation without re-estimating density or holdout.
3. Extend TFit export and explicit solver sample acceptance as described above.
4. Produce all required numerical tables, maps, clouds, and delta-organized pages.
5. Document one short reproducible campaign/tag invocation in `CORE_SAMPLE.md`.

A counts preview cannot produce event clouds or establish event-level membership.
The local repository currently has full original candidate files only for rg01;
full campaign core masks/replay inputs may need to be used on ifarm. The worker
must check availability early. If data are missing, report exact required paths;
do not substitute historical GMM trees or simulations for missing core samples.

Snapshot all policy inputs and chosen IDs, and hash the originating sample,
geometry, metadata, and code. Record which limits bind and why statistics were
lost. Output tags must be immutable; a changed factor produces a new tag.

## Targeted acceptance checks

- Exact equality of final totals across all required physical foils.
- Exact budget conservation; every quota respects fixed delta/hole capacities.
- P25 references come from original positive eligible counts and exclude blocked
  positions and holdout. No hidden legacy 80%/400-event limits remain active.
- Correct physical-foil mapping and aligned delta boundaries; distinct settings
  and event identities preserved.
- No duplicate, blocked, unsupported, noncore, or protected-holdout events in
  training. Surplus and all event categories reconcile to the frozen input.
- Identical quotas and IDs after reordering input rows/files; stable remainder
  handling with fixed seed. Test empty/tiny populations and tied counts.
- A bounded fixture with dense and sparse foils/deltas/holes demonstrates exact
  foil totals and fixed caps without automatic relaxation or weighted events.
- Event-level quotas agree with the counts preview. A small preallocated SVD
  fixture admits exactly the supplied IDs even when the legacy caps would reject
  some of them; legacy mode remains unchanged.
- The diagnostic index includes all five physical foils for each delta slice,
  labels empty entries, and exposes both aggregate and setting-level views.

Return the implementation, a brief verification report, the full diagnostic
index, and a concise description of the main remaining imbalances. Do not claim
that count equality proves equal leverage, good conditioning, or a stable optics
matrix. Those are assessed in the subsequent fit study after this selection
policy has been reviewed.

## Accepted clarifications

- `qa_min: 10` flags populated leaves with fewer than 10 training events, including zero. Reporting only; never alters quotas.
- Any required foil with zero feasible capacity yields a zero campaign budget. Report its limiting constraints and stop before export or solver use, without removing the foil or relaxing caps.
