# Strict core balance: implementation verification

Implemented in `core_balance.py`, `reallocate_core.py`, `balance_diagnostics.py`,
`build_core_fit.py`, `preallocated_svd.py`, and the explicit preallocated path in
`fit_opt_matrix_gmm.C` / `preallocated_sample.h`. The command and configuration
reference is [CORE_SAMPLE.md](../CORE_SAMPLE.md#strict-frozen-tag-balance-september-2026).

[Complete campaign counts-preview index](../HMS_6p667GeV/05d_core_balance/equal15_review/HMS_6p667_CORE_BALANCE.md)
— all 25 physical-foil/delta and 60 setting/local-foil/delta views, five cross-foil
pages, numerical tables, log circle maps and fraction maps. Event clouds are
explicitly unavailable in this preview; none are inferred from historical GMM.

## Campaign counts result

| Physical foil (cm) | Eligible events | Feasible after fixed caps | Planned training | Planned surplus |
|---:|---:|---:|---:|---:|
| -8 | 54,299 | 33,294 | 13,707 | 40,592 |
| -3 | 23,339 | 14,475 | 13,707 | 9,632 |
| 0 | 160,951 | 92,188 | 13,707 | 147,244 |
| +3 | 25,088 | 13,707 | 13,707 | 11,381 |
| +8 | 46,836 | 24,862 | 13,707 | 33,129 |

Exactly **68,535 planned training events**, limited by the **+3 cm foil after
hole caps**. Original eligible availability is 310,513. Hole ceilings remove
131,987 events from feasible capacity; delta ceilings cause no additional loss;
strict foil equality leaves another 109,991 feasible events unused. Total
planned surplus is 241,978. The 200,000 campaign ceiling is not binding.

The disjoint frozen count partition reconciles as
499,868 = 310,513 eligible development core + 101,258 protected holdout + 88,097
noncore/unsupported development. There are 667 development events at blocked
labels, already noncore/unsupported in this parent, so the geometry mask removes
no additional accepted core events in this trial. Unsupported and holdout-core
columns are subsets, not extra disjoint categories to add to the partition.

There are **53 positive available allocation leaves below `qa_min=10`**; none
receive zero events in this campaign preview. Tests also exercise positive leaves
receiving zero. QA does not change quotas. All seven rungroups retain nonzero
contributions; totals range from 5,775 to 18,848. Distinct setting leaves are
preserved throughout allocation.

## Main remaining imbalances

Delta slice 0, [-10,-8)% (width 2 percentage points):

| Physical foil (cm) | Available | Capped slice capacity | Planned training |
|---:|---:|---:|---:|
| -8 | 4,538 | 2,794 | 2,741 |
| -3 | 1,965 | 1,495 | 1,495 |
| 0 | 14,667 | 10,338 | 2,741 |
| +3 | 2,123 | 1,430 | 1,430 |
| +8 | 5,712 | 4,001 | 2,742 |

The -3 and +3 cm slice-0 capacities bind. Equal foil totals deliberately do not
force foil-by-delta equality; these foils use more of their budgets in other
slices. Settings still contribute differently and dense populations are thinned
strongly, especially at 0 cm. Slice widths remain 2, 3, 5, 5, 5 percentage points,
without width weighting. Count equality establishes neither equal leverage nor
good conditioning or a stable optics matrix.

Historical actual SVD foil totals, in target order, are 44,310 / 25,505 / 30,000 /
25,731 / 40,467 (166,013 total). These are the saved, log-verified historical GMM
selection, not the new core population. Original historical solver ID lists were
not saved; the preserved historical diagnostics document their reconstruction.
No historical output was changed.

## Targeted verification

18 tests passed across the four documented suites:

- Seven allocation tests cover interpolated positive-count P25, fixed capacities,
  exact foil/parent/child budgets, stable integer ties and ID selection after row
  and setting reordering, 64-bit identities above 2^54, empty/tiny populations,
  zero feasible foils, sub-foil campaign ceilings, QA-only thresholds, exclusions,
  protected holdout, configuration rejection and the complete campaign preview.
- Nine existing selector tests remain passing, including legacy allocation,
  holdout independence, density behavior and ROOT sample output.
- One frozen-tag integration test verifies shared preview/event quotas, repeat
  reallocation IDs, unchanged quality and models, unchanged holdout, real cloud
  rendering, immutable tags, metadata mismatch rejection, a zero-budget export
  stop and HMS TFit export with exact membership (588 synthetic events).
- One real ROOT solver integration test admits all 18,018 supplied synthetic IDs
  through the current preallocated solver admission/design-row path, then stops
  before SVD. Each of 18 Y/delta populations contains 1,001 events, exceeding the
  old limits. The same fixture in legacy mode admits 15,000, as expected. Resource
  and memory ceilings fail rather than truncate; an invalid supplied event is
  audited and rejected; a changed build hash fails before ROOT invocation.

The legacy-mode regression used only a bounded synthetic fixture and disposable
outputs. **No campaign production SVD fit was run and no matrix was replaced.**
The preallocated fixture wrote no matrix. ROOT loaded the modified macro; its
existing `sprintf` deprecation warnings remain unrelated to these changes.

## Remaining data-dependent work

The local `min10` mirror lacks seven event-mask ROOT files, 60 frozen model files,
two archived source/default files and seven upstream exclusion audit TSVs.
Replay inputs for all seven rungroups are also absent. Exact local and originating
ifarm paths are in the index's
[INPUT_AVAILABILITY.md](../HMS_6p667GeV/05d_core_balance/equal15_review/INPUT_AVAILABILITY.md).

Consequently, **campaign event IDs, event-level reconciliation, real allocation
clouds, campaign TFit export and campaign solver admission are not verified or
produced here**. The counts preview is complete; the event implementation is
fixture-tested. Run the documented frozen-tag command on the complete saved
inputs to finish those data-dependent outputs without re-estimating density.
SHMS adapter/geometry support was not validated and is not claimed.
