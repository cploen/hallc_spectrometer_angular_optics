# Frozen matrix comparison

Compare the starting replay, historical GMM fit, full-basis core SVD, compact
elastic-net-selected SVD refit, and beam result on identical events. This step
never fits coefficients or recenters residuals. It evaluates all valid saved
polynomial terms, including fixed xtar terms and fitted constants.

## Run on ifarm

Use the same campaign and balanced sample that supplied beam10. If the surplus
TFit export does not already exist, build it once:

```bash
./run_build_core_fit.sh HMS_6p667GeV equal15 surplus
```

Keep the existing fit and holdout exports. Then:

```bash
./run_compare.sh HMS_6p667GeV equal15 compare --source beam10 --check
./run_compare.sh HMS_6p667GeV equal15 compare --source beam10 --threads 8
```

The check verifies matrices, checksums, exact event membership and separation
from modern training; it does not evaluate residuals. Outputs go to
`HMS_6p667GeV/07_diagnostics/compare/`. Existing outputs cannot be overwritten;
use a different short output name for another explicitly documented comparison.
Local Mac evaluation requires the allocation masks/manifest and both evaluation
TFit exports. PNGs and beam matrices alone are insufficient.

## Refresh existing plots quickly

For the completed `compare` result, run on ifarm, where `residuals.npz` is saved:

```bash
./run_compare.sh HMS_6p667GeV equal15 compare --replot
```

This checks the saved residual-array and table checksums, then redraws the plots
in the existing `compare/plots/` directory. It does **not** reload ROOT files,
evaluate matrices, refit, or repeat the slow grouped statistics. It prints
progress for each event pool. Keep the saved `residuals.npz` on ifarm; PNGs and
summary tables alone cannot reproduce the residual histograms.

`--replot` explicitly replaces the figures, updates their manifest checksums,
and saves the plotting code under `plot_code/`. The original numerical code
snapshot, matrices, offsets, residual arrays and TSVs remain unchanged. The
plot guide is `plots/README.md`. The historical report's old plot description
is replaced with a link to that guide. This command does not fix the uncertain
legacy offset convention or the slow statistics calculation in a fresh run.

## Offsets and the historical file

`config/comparison.json` specifies external additions in this order:
`[xptar mrad, ytar cm, yptar mrad]`. Each model requires a source note. Values can
be shared by all rungroups or supplied as a rungroup-to-vector mapping. Null or
missing values fail preflight. No correction is estimated from evaluation data.

Christine supplied the legacy pre-rotation/clustering/GMM offsets on 2026-09-13:

- hphi_offset = 7.11802209e-4 rad: add 0.711802209 mrad to xptar.
- htheta_offset = 5.84959117e-4 rad: add 0.584959117 mrad to yptar.
- No external ytar correction is supplied.

These additions match the ordering in the upstream
[HCANA target reconstruction](https://github.com/JeffersonLab/hcana/blob/master/src/THcHallCSpectrometer.cxx).
Commented constant rows in the starting matrix are not active coefficients.
The GMM, full core, compact and beam matrices have active fitted constants;
the policy adds zero external correction to these solved matrices. We do not
add the old offsets again. This follows replay without offsets followed by
solving for them; it does not reproduce every historical runtime flag.

The archived GMM output has one corrupt exponent row from the documented old
parser bug (see `08_preliminary_conditioning/PLAIN_LANGUAGE_CLARIFICATIONS.md`,
section 10). The policy explicitly excludes that exact row, pins the entire
original file checksum, and explains why. All valid historical coefficients
remain unchanged. This is **not** the corrected-parser refit. No arbitrary
malformed rows are silently skipped. The original file and exclusion policy
are archived in each comparison. Historical HCANA behavior on the corrupt row
is not established by this evaluation.

## Read the outputs

Start with `MATRIX_COMPARISON.md`, then:

- `plots/*_center.png`: signed residuals near zero, common bins, **events per
  bin (%)**, no recentering. Window is ±beam P90. Each curve lists its visible
  fraction; percentages use the full pool so clipping cannot hide event loss.
- `plots/*_tails.png`: percentage of events at or beyond each absolute residual.
  Linear horizontal axis out to the largest P99.9 among the matrices; logarithmic
  percentage axis with 10%, 1%, and 0.1% guides. For example, 1% at 3 mrad means
  1% of the pool has |residual| ≥3 mrad. Lower is better.
- `plots/*_extremes.png`: separate full-range tails on logarithmic axes, including
  the rare events beyond the main tail window. Neither tail view renormalizes
  when zoomed.
- `plots/*_foil_delta.png`: all physical foils and delta slices, all five models,
  absolute RMS with common color scale per target and counts in each cell.
- `plots/*_change.png`: percent RMS change for each consecutive stage, plus beam
  versus full core and beam versus GMM. Blue/negative means smaller RMS;
  red/positive means larger RMS. Colors saturate at ±20%, but printed values
  are not clipped. Undefined changes (missing cells or zero reference RMS) show
  “—”. Counts are shown; N<10 has an asterisk and never changes allocation.
- `tsv/summary.tsv` and `changes.tsv`: bias, mean-subtracted spread, RMS, median
  absolute residual, P90 absolute residual and equal-cell MSE; changes relative
  to both the previous stage and the starting matrix.
- `tsv/residuals.tsv`: the same residual statistics by foil/delta, rungroup and
  hole. N<10 is a reporting flag only. No sparse events are discarded here.
- `tsv/overlap.tsv`: historical GMM training overlap; unknown membership stays
  unknown. Summary tables also report a common subset outside reconstructed
  GMM training membership.
- `tsv/training_conditioning.tsv`: frozen modern training-design diagnostics
  from beam10. These are not condition numbers inferred from residuals.

Angular residuals are at the target in mrad; ytar and derived ztar are in cm.
RMS includes bias; spread removes only the mean mathematically and is not a
Gaussian fit. Comparing both prevents an offset shift from masquerading as
better width. P90 is the absolute error below which 90% of events fall.

Protected cores, protected noncore regions, and surplus cores are separate.
Unsupported labels are tabulated separately; blocked positions have counts
only. Surplus is development coverage. Old-matrix training membership is
unknown, and GMM membership was reconstructed, so historical comparisons are
not an independently held-out test of every model. Noncore interpretation also
depends on label correctness.

Old to GMM measures the historical change. GMM to full core combines selection,
balancing and solver changes; it cannot isolate core selection alone. Full core
to compact to beam keeps modern training and the scaled direct-X SVD method
fixed. The compact model is `beam10/matrices/start.dat`, the unpenalized full
training refit, not elastic net's penalized coefficients.

All models use common saved focal-plane and xtar inputs and saved geometry.
This comparison is not a full iterative HCANA replay. No matrix is tuned here.
If protected results guide more tuning, subsequent use is development evidence.

## Reproduce or use another HMS campaign

Keep the complete output directory: manifest, original matrices, external-offset
policy, code snapshots, residual arrays, tables and plots. The manifest records
input/output checksums. Large residual arrays need not be pushed to inspect
plots and summary tables; preserve them on ifarm. Expose the report, manifest,
plots and TSVs using the repository's existing output-file allowlist convention.

Provide that campaign's own `config/comparison.json`, completed beam output,
frozen allocation, holdout/surplus exports and historical GMM membership. Verify
historical offsets and matrix identity for that campaign. Do not copy 6.667
values or row exclusions blindly. The TFit adapter currently accepts HMS only;
SHMS evaluation requires its geometry/reconstruction adapter first.
