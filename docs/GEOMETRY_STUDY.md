# Core-event geometry study

## Question

Does the saved replay angle follow the ray geometry implied by an assigned sieve
hole and the measured vertical vertex? Start with run 1544, zero foil, central
hole. Compare representative holes at the same run setting. No matrix or target
definitions are changed by this diagnostic.

`reacty` is the laboratory vertical reaction coordinate in cm. `xptar` is the
saved replay transport-x slope. `xptarT` is the target calculated by the existing
`hallc::targetTruth` routine for the assigned foil and hole. The Beamer
[proof](hms_geometry_proof.pdf) defines the variables and shows the contradiction.

For HMS the constructed target pair gives

```
xtarT + L*xptarT - xsT = -reacty - xmis
```

Under the same longitudinal approximation, the target correction is
`(reacty+xmis)/(L-zfoil*cos(angle))`. The study plots
`1000*(xptar-xptarT)` versus `reacty`, then subtracts that correction from the
residual and repeats the plot. At zero foil the predicted original slope is
`1000/168 = 5.95238 mrad/cm`; the corrected slope would be zero **if the replay
angle follows this reference geometry**. SHMS already includes the vertical
vertex in its target slope: the diagnostic applies no extra correction, and
its reference residual slope is zero.

A flat original residual does not disprove the target inconsistency. A matrix
trained toward the old targets could suppress the expected dependence.
Subtracting a linear correction necessarily changes the fitted slope by its
coefficient; that subtraction alone supplies no independent validation.

## Inputs and selection

- `config/rungroups_*_inputs.tsv` supplies group names, optics IDs, angles and
  original replay paths. The campaign name selects HMS or SHMS through the
  existing `spectrometer_config.py` and shared profiles.
- `config/geometry_study.json` supplies the core-file pattern, output path,
  physical foil position and selection settings. The HMS configuration uses
  the existing `05c_core_sample/min10/root/CoreSample_{rungroup}.root` files.
- The optics metadata supplies foil locations. Each group and physical foil
  gets separate outputs; angles and groups are never pooled. Delta-slice fits
  accompany the pooled fit to expose momentum dependence.
- Use only `core_keep == 1`, including the fit, holdout and surplus partitions.
  This is descriptive inspection of all cores; these inspected holdout events
  should not later be presented as an untouched validation set for a choice
  informed by these results.
- Retain the established Cherenkov > 2, calorimeter > 0.65, arm-specific delta
  acceptance, assigned delta slice, and reaction-z within 2 cm of the foil.
  There is no additional reconstructed-sieve cut or angle-residual cut.
- Default holes are the center and two index steps in either direction along
  each sieve axis. For HMS these are `(4,4), (2,4), (6,4), (4,2), (4,6)`.
  A `holes` list of `[xscol,yscol]` pairs can override this in the configuration.
- When the original replay is available, join by the saved merged-tree `entry`
  and verify agreement in all four focal-plane coordinates and delta. Otherwise
  use the replay values already saved in `CoreSample`. The terminal and manifest
  explicitly record this choice. Missing saved branches are an error.

Christine confirmed that the 2024 matrix matches NPS replay, checked with the
afterburner closure test; all zero-order offsets were removed for that replay,
and the HCANA version is unchanged. This study reads that existing reconstruction.

## Selection-sensitivity comparison

The core selection uses focal-plane labels reconciled with sieve geometry.
It can therefore retain dependence on reconstructed sieve position.

Within each hole and delta slice, compute the mean and covariance of the four
measured focal-plane coordinates among retained cores. Keep the inner 50% by
covariance-normalized distance from that mean. Singular directions are omitted.
This additional restriction uses no sieve coordinate, vertex or reconstructed
angle. It is still a restriction on physical phase space: it can change a real
vertex-angle relationship as well as selection bias. Agreement would demonstrate
some stability; disagreement cannot by itself identify which mechanism changed.
Previously rejected events remain absent in both samples.

## Plots and fits

ROOT `TH2D` histograms contain event counts, with a fitted straight line overlaid.
Fits use individual events (unweighted ordinary least squares via `TGraph`),
not histogram-bin centers or profile errors. Each plot displays the measured
slope, its statistical uncertainty and event count above the data. No theoretical
line or predicted slope appears on the plots. The terminal reports counts,
vertex range/RMS, slope/error, residual RMS, predicted slope and standardized
differences. Errors are statistical OLS errors, not selection uncertainties.

`slopes.tsv` records those quantities and the intercept. `ndel=-1` means pooled
delta slices. `geometry.root` stores the histograms, fitted functions and canvases.
Each figure also has PNG and PDF versions. Empty/underpopulated holes are reported
and skipped. The default minimum is 30 events per fit.

## Commands

From the repository root, with Python 3 and ROOT available:

```bash
# Inputs and paths only
bash run_geometry_study.sh HMS_6p667GeV --check

# First study: central foil, run 1544 and five representative holes
bash run_geometry_study.sh HMS_6p667GeV --rungroup rg01_theta12p490_foil0

# Other zero-foil groups remain separate
bash run_geometry_study.sh HMS_6p667GeV

# Extension to an outer foil at the matching 12.490 degree setting
bash run_geometry_study.sh HMS_6p667GeV --foil 8 --rungroup rg03_theta12p490_foilpm8
```

For another campaign, put a `geometry_study.json` in its existing `config/`
directory with its core-file pattern, or pass `--config PATH`. No run or replay
path is embedded in the analysis code. `OPTICS_METADATA` can override the existing
run-metadata file. Results go to the configured output directory under
`foil_<position>cm/<rungroup>/`. Rerunning a group overwrites its matching products.

## First result: run 1544, zero foil

Local run on 19 September 2026, using saved replay fields in the existing core
file because the original `/volatile/...` replay is on iFarm. Fourteen events
across the five holes failed the additional quality/foil checks. The target-pair
identity agrees with its algebraic expression to `6.1e-16 cm`; this is an internal
arithmetic check, not a measurement of physical correctness.

All slopes below are mrad/cm. Uncertainties are statistical.

| (xscol,yscol) | Core N | Original slope | Corrected slope | Inner FP-core original slope |
|---|---:|---:|---:|---:|
| (4,4) | 1584 | 8.73 +/- 0.61 | 2.78 +/- 0.61 | 0.51 +/- 0.82 |
| (2,4) | 4357 | 8.24 +/- 0.37 | 2.28 +/- 0.37 | -1.26 +/- 0.51 |
| (6,4) | 4293 | 8.97 +/- 0.35 | 3.02 +/- 0.35 | 0.48 +/- 0.47 |
| (4,2) | 2449 | 7.77 +/- 0.52 | 1.82 +/- 0.52 | -0.53 +/- 0.69 |
| (4,6) | 9322 | 7.40 +/- 0.25 | 1.45 +/- 0.25 | -1.62 +/- 0.33 |

Central-hole `reacty` RMS is 0.05333 cm. Its residual standard deviation changes
from 1.37060 to 1.29756 mrad after the full target correction, including `xmis`.
The inner FP subset has 791 events, with original slope 0.51 +/- 0.82 and
corrected slope -5.44 +/- 0.82 mrad/cm. The central-hole delta-slice original
slopes range from 2.42 +/- 2.72 to 10.81 +/- 1.53 mrad/cm.

**Interpretation:** the retained full cores have a reproducible vertex dependence
that a constant angular offset cannot remove. They do not isolate the proposed
missing term: the correction leaves nonzero slopes, and changing the FP selection
changes the slopes strongly. The algebraic target mismatch remains established
under its stated conventions. These data do not yet establish its contribution
to replay resolution or justify a target replacement/refit.

Results: [numerical table](../HMS_6p667GeV/07_diagnostics/geometry/foil_0cm/rg01_theta12p490_foil0/slopes.tsv),
[terminal output](../HMS_6p667GeV/07_diagnostics/geometry/foil_0cm/rg01_theta12p490_foil0/terminal.txt),
[central-hole original](../HMS_6p667GeV/07_diagnostics/geometry/foil_0cm/rg01_theta12p490_foil0/x4_y4_core_all_delta_original.pdf),
[central-hole corrected](../HMS_6p667GeV/07_diagnostics/geometry/foil_0cm/rg01_theta12p490_foil0/x4_y4_core_all_delta_corrected.pdf).

## Validation and next step

### Follow-up: densest half of the central-hole core

This is a sieve-density selection, distinct from the inner focal-plane subset.
Within each delta slice of the retained central-hole core, rank `core_score`
from highest to lowest and retain the top half, including all ties at the
threshold. `core_score` is the background-subtracted reconstructed-sieve density
relative to the peak. The selections retain 93, 161, 246, 189 and 119 events in
delta slices 0 through 4, respectively, for 808 events total.

```bash
bash run_geometry_study.sh HMS_6p667GeV --rungroup rg01_theta12p490_foil0 --central-dense-half
```

The command explicitly selects only the central hole. Results go beneath
`07_diagnostics/geometry/foil_0cm/central_dense_half/rg01_theta12p490_foil0/`.
It also regenerates the full-core central-hole comparison with the same cuts.

| Selection | Events | Original slope (mrad/cm) | Corrected slope (mrad/cm) |
|---|---:|---:|---:|
| Full retained central-hole core | 1584 | 8.7313 +/- 0.6077 | 2.7789 +/- 0.6077 |
| Densest half, ties included | 808 | 6.7440 +/- 0.4818 | 0.7916 +/- 0.4818 |

The densest-half slope differs from 5.95238 mrad/cm by 1.64 statistical standard
errors. Its residual standard deviation changes from 0.82537 to 0.74152 mrad
after correction. This supports agreement with the expected dependence under
this selection. It does not independently establish the cause, because the
selection uses reconstructed sieve coordinates and the slope is selection
dependent. Densest-half results have not been extended to the other holes.

For Charles, share the original and corrected densest-half PDFs, the original
full-core PDF for comparison, and `docs/hms_geometry_proof.pdf`. All are committed
with their PNG versions; the numerical tables, terminal output and configuration
are retained too. ROOT histogram files and delta-slice figures can be regenerated
with the command above and are not checked in.

Suggested accompanying statement: the target equations have an algebraic
mismatch. The densest-half central-hole result is consistent with the predicted
vertex slope within statistical uncertainty, but the difference from the
full-core result shows that selection matters. No target change or matrix refit
has been made.

`python3 -m unittest discover -s tests -p test_geometry_study.py` checks synthetic
reference rays for HMS and SHMS at zero and +8 cm. It verifies known slopes,
the corrected zero slope, arm selection, and rejection of noncore/other-foil
events. The real zero-foil run completed and representative plots were inspected.
No real outer-foil or SHMS data study has been run.

The next useful physics check is to compare these retained cores with the
original focal-plane-assigned central-hole candidates before sieve-based core
selection. That would change the current core-only scope and needs a separate
decision. No such selection expansion was made here.
