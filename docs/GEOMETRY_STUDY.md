# Core-event geometry study

## Observable and prediction

Plot the unchanged saved replay residual `1000*(xptar-xptarT)` against `reacty`
(in cm). `xptarT` comes from the existing `hallc::targetTruth` routine for the
assigned foil and hole. The factor 1000 expresses the residual in mrad.
No predicted slope or vertex term is subtracted from any event.

The HMS target definitions satisfy

```
xtarT + L*xptarT - xsT = -reacty - xmis
```

Under the same longitudinal approximation, a replay angle following the physical
reference ray would produce a residual slope `1000/(L-zfoil*cos(angle))`.
At zero foil this is 5.95238 mrad/cm. SHMS targets already include the vertex term,
so the reference residual slope is zero. The Beamer derivation is in
[docs/hms_geometry_proof.pdf](hms_geometry_proof.pdf).

The plots contain data histograms, an unconstrained fitted line, ROOT's native
fit-statistics box, and a subtitle giving the geometry prediction. The fitted
line comes from unweighted ordinary least squares on individual events using
`TGraph`; histogram binning does not determine the fitted coefficients. The
statistics box shows entries, intercept, slope and their errors, plus ROOT's
chi-square/ndf. With no individual measurement errors supplied, this chi-square
uses unit weights and is not a calibrated goodness-of-fit probability. Reported
parameter errors are statistical OLS errors and exclude selection systematics.

## Inputs and selection

The campaign name selects HMS or SHMS via the existing spectrometer profiles.
`config/rungroups_*_inputs.tsv` supplies run groups, optics IDs, angles and replay
paths. `config/geometry_study.json` supplies the core-file pattern and output
location. Foil positions come from the existing run metadata, with an optional
`OPTICS_METADATA` environment override. Each group and foil stays separate.
Momentum-slice fits accompany the pooled fit.

Use only `core_keep == 1`, including fit, holdout and surplus events, with the
established Cherenkov > 2, calorimeter > 0.65, arm-specific delta acceptance,
assigned delta interval and reaction-z within 2 cm of the selected foil.
These are descriptive results; inspected holdout events should not be called
untouched validation for subsequent choices informed by these results.

When accessible, the original replay is joined by entry and checked against
saved focal-plane coordinates and delta. Otherwise use the replay fields saved
in CoreSample and record this in the terminal output and manifest. Run 1544 was
analyzed locally from these saved fields because its configured replay path is
on iFarm. Christine confirmed the 2024 matrix matches NPS replay, all zero-order
offsets were removed for that replay, and the HCANA version is unchanged.

Default holes are the center and two indices in each direction along both sieve
axes. HMS indices are `(4,4), (2,4), (6,4), (4,2), (4,6)`. Configuration can override
these with a `holes` list. Underpopulated selections are skipped (minimum 30).

Two distinct selection comparisons are available:

- **Inner FP core:** retain the inner 50% by covariance-normalized distance in
  four focal-plane coordinates, independently within each hole/delta slice.
  This restriction can change physical vertex-angle correlations as well as bias.
- **Densest half:** central hole only, highest `core_score` half in each delta
  slice, including ties. The score is reconstructed-sieve density above background
  relative to the peak. This is not a focal-plane distance selection.

Both begin with the existing sieve-reconciled cores. Neither restores previously
rejected events or supplies an independent test of the initial core selection.

## Reproduction

From the repository root with Python 3 and ROOT:

```bash
bash run_geometry_study.sh HMS_6p667GeV --check
bash run_geometry_study.sh HMS_6p667GeV --rungroup rg01_theta12p490_foil0
bash run_geometry_study.sh HMS_6p667GeV --rungroup rg01_theta12p490_foil0 --central-dense-half
```

Other zero-foil groups can be run separately by omitting `--rungroup`. To extend
to an outer foil at the matching 12.490 degree setting:

```bash
bash run_geometry_study.sh HMS_6p667GeV --foil 8 --rungroup rg03_theta12p490_foilpm8
```

Another campaign uses its own `config/geometry_study.json`, or `--config PATH`.
No replay path is hardcoded in the analysis. Rerunning overwrites matching products.
ROOT files contain histograms, fit functions and canvases; PNG/PDF exports and
`slopes.tsv`, `terminal.txt`, `manifest.json` accompany them. Pooled plots and
summaries are checked in. ROOT files and momentum-slice plots are regenerable.

## Results: run 1544, zero foil

| Hole | Core events | Fitted slope (mrad/cm) | Inner FP-core slope (mrad/cm) |
|---|---:|---:|---:|
| (4,4) | 1584 | 8.73 +/- 0.61 | 0.51 +/- 0.82 |
| (2,4) | 4357 | 8.24 +/- 0.37 | -1.26 +/- 0.51 |
| (6,4) | 4293 | 8.97 +/- 0.35 | 0.48 +/- 0.47 |
| (4,2) | 2449 | 7.77 +/- 0.52 | -0.53 +/- 0.69 |
| (4,6) | 9322 | 7.40 +/- 0.25 | -1.62 +/- 0.33 |

The densest-half central-hole sample contains 808 events and gives
**6.7440 +/- 0.4818 mrad/cm**, compared with 5.95238 mrad/cm predicted.
The difference is 1.64 statistical standard errors. Its reacty RMS is 0.05412 cm
and residual standard deviation is 0.82537 mrad.

The highest-score half retains 93, 161, 246, 189, 119 events in delta slices 0--4,
with thresholds 0.8249194, 0.8547786, 0.8580519, 0.8201271, 0.810554, respectively.
Ties at each threshold are retained. This selection has been applied only to the
central hole. The five-hole full-core study has 22,005 events after cuts.

The dense central sample agrees with the prediction within statistical
uncertainty. A constant angular offset cannot remove its fitted slope. However,
the substantial selection dependence prevents identifying the cause uniquely.
The algebraic mismatch is established under the stated conventions; its effect
on fitted/replayed resolution is not established. No targets or matrices were changed.

## Plots for Charles

Under `HMS_6p667GeV/07_diagnostics/geometry/foil_0cm/central_dense_half/rg01_theta12p490_foil0/`:

- `x4_y4_core_all_delta.pdf`: full central-hole core, N=1584.
- `x4_y4_dense_half_all_delta.pdf`: densest half, N=808.

PNG versions are alongside the PDFs. These are the two data-only central-hole
comparisons to share with the geometry proof. All superseded plot variants were
removed from the current branch. There are no plots with a predicted term subtracted.

## Checks and next step

Synthetic ROOT tests verify known slopes, HMS/SHMS selection, zero and +8 cm foil
isolation, noncore rejection and dense-half selection with ties:

```bash
python3 -m unittest discover -s tests -p test_geometry_study.py
```

The real zero-foil run completed and representative renderings were inspected.
Real outer-foil and SHMS data studies have not been run. The next discriminating
physics comparison would use the original focal-plane-labeled central-hole events
before sieve-based core selection, which requires a separate scope decision.
