# HMS 6.667 — frozen matrix comparison and revised diagnostics

Conversation recap, 15 September 2026. Runbook:
[MATRIX_COMPARISON.md](../../docs/MATRIX_COMPARISON.md).
Numerical record: [compare/MATRIX_COMPARISON.md](compare/MATRIX_COMPARISON.md),
`compare/tsv/summary.tsv`, `changes.tsv`, `residuals.tsv` and `overlap.tsv`.

## Which matrices and events were compared

| Label | Identity |
|---|---|
| Old | `config/oldfit.dat`, with explicit legacy external offsets in this comparison |
| GMM | Historical `06b_svd_fit/matrices/nps_hms_newfit_6p667_gmm_clean.dat`; original selection and ROOT solve, **not** GMM events refitted using the beam basis |
| Full core | `06e_beam_search/beam10/matrices/svd.dat` |
| EN-selected | `06e_beam_search/beam10/matrices/start.dat`; unpenalized compact-basis refit |
| Beam | `06e_beam_search/beam10/matrices/beam.dat` |

Every matrix sees identical saved focal-plane and xtar inputs, including its
valid fixed terms and constants. No fit is performed. Pools: 77,864 protected
cores, 22,532 protected noncores, 638 protected unsupported-region events, and
241,978 surplus cores; 224 blocked events are counted separately. Noncores
combine shoulders and surrounding labeled watershed regions; these are not two
separately identified tests. Unsupported means no accepted core region for the
inherited label, not an “unsupported core.” Label errors/scattering can affect
these residuals, so noncore tails are not all attributable to the matrix.

Historical GMM training overlap was reconstructed and reported; old-matrix
training membership is unknown. Consequently the all-events comparison is not
an independently held-out test of every historical matrix. Surplus cores are
development coverage, and inspected protected results cannot remain a pristine
test for future tuning.

## Interpretation established in this conversation

Protected-core RMS (target angles in mrad, positions in cm):

| Matrix | xptar | ytar | yptar | derived ztar |
|---|---:|---:|---:|---:|
| GMM | 1.5367 | 0.13726 | 1.0344 | 0.64396 |
| Full core | 1.4886 | 0.13582 | 1.0292 | 0.63741 |
| EN-selected | 1.5538 | 0.13713 | 1.0363 | 0.64389 |
| Beam | 1.4956 | 0.13602 | 1.0329 | 0.63853 |

The clear beam-over-compact-seed gain is xptar (about 3.7% pooled RMS, improvement
in all 25 foil/delta cells). Beam versus GMM gains about 6–14% in delta-0 core
xptar across the five foils. Much of that gain already appeared at the full-core
stage; selection, balance and solver effects were not isolated. There is no
broad pooled shoulder win over GMM and no global beam win over full-core SVD.
Smaller basis/better conditioning alone is not proof of superior reconstruction.

RMS contains bias; spread subtracts the mean and is not a Gaussian-fit sigma.
Use the per-foil/cell tables rather than calling the pooled ztar RMS a resolution
at every foil. Angular residuals are at the target. These comparisons and the
derived ztar diagnostic use fixed inputs/truth geometry, not a complete replay.

## Offsets and plot changes

Legacy additions were supplied as hphi=7.11802209e-4 rad to xptar and
htheta=5.84959117e-4 rad to yptar. The solved matrices have fitted constants and
receive no extra legacy correction. **The old-matrix offset convention remains
unresolved**; do not attribute gains simply to adding a constant. The historical
GMM's one malformed exponent row is excluded only through an exact row/checksum
policy; no corrected-parser refit was substituted.

Plot fixes (`2634d39`) change presentation only: central percentages per bin,
linear-horizontal tail probabilities with 10%/1%/0.1% guides, separate extremes,
and percentage RMS-change maps. `--replot` uses saved residuals and tables,
preserving numerical outputs and avoiding the slow grouped-statistics pass.
The farm refresh has not been confirmed in this audit. A fresh comparison still
has the earlier statistics bottleneck; that optimization was deferred.

Code: `compare_matrices.py`, `comparison_plots.py`. Seven targeted tests passed
when the replot change was developed; `tests/test_comparison.py` covers parser/
offset behavior, exact historical IDs, no-evaluation preflight, common subsets,
tail probabilities and replot checksum/numerical preservation. The full
`compare/residuals.npz` remains on ifarm and is required to redraw histograms.
