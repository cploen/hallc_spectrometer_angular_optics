# Reproducing the preliminary 6.667 GeV angular-optics study

This is the runbook for the work performed during the preliminary angular
optics conditioning discussion on 26--27 August 2026.  It records the data
provenance, the parser ambiguity we found, the order of calculations, and the
commands corresponding to the committed scripts.

The delta analysis was used only as the organizational and plotting model.  No
delta fit or delta study was rerun.

## Repository state and scope

The two repositories were checked out locally:

- `hallc_spectrometer_angular_optics`, branch `prelim-6p667-data`
- `hms-delta-fit-diagnostics`, branch `main`, commit `60a631d`

The angular branch's input-data commit is `5cb66ee` (`Add minimal 6.667 GeV
angular fit reproduction inputs`).  That commit had already been pushed from
ifarm before this analysis began.  The ifarm source directory was
`~/nps_link/ML_HMS_OPTICS_DEV`.

The preliminary study was intentionally limited to the existing 6.667 GeV fit
ntuples.  It did not regenerate the GMM selections or fit ntuples, and it did
not rerun any delta study.

## Software environment used locally

The recorded local environment was:

```text
macOS arm64
ROOT / PyROOT 6.36.000
Python 3.13.5
NumPy 2.3.1
```

The scripts require ROOT with PyROOT and NumPy.  Run all commands below from
the root of `hallc_spectrometer_angular_optics`.

## Inputs

The direct inputs are:

```text
DATfiles/list_of_optics_run.dat
HMS_6p667GeV/config/oldfit.dat
HMS_6p667GeV/config/rungroups_6p667_inputs.tsv
HMS_6p667GeV/06a_fit_ntuple/root/Optics_666701_-1_fit_tree_gmm.root
HMS_6p667GeV/06a_fit_ntuple/root/Optics_666702_-1_fit_tree_gmm.root
HMS_6p667GeV/06a_fit_ntuple/root/Optics_666703_-1_fit_tree_gmm.root
HMS_6p667GeV/06a_fit_ntuple/root/Optics_666704_-1_fit_tree_gmm.root
HMS_6p667GeV/06a_fit_ntuple/root/Optics_666705_-1_fit_tree_gmm.root
HMS_6p667GeV/06a_fit_ntuple/root/Optics_666706_-1_fit_tree_gmm.root
HMS_6p667GeV/06a_fit_ntuple/root/Optics_666707_-1_fit_tree_gmm.root
HMS_6p667GeV/06b_svd_fit/matrices/nps_hms_newfit_6p667_gmm_clean.dat
```

`pre_gmm_veto.tsv`, the original fit logs, and the build summaries are retained
with the branch as provenance, but the conditioning scripts operate on the fit
ntuples rather than rerunning the preceding GMM workflow.

SHA-256 checksums for the numerical inputs used here are:

```text
0199c6fed3faa0c567b6e3a0fd016662acdc3c5c46cd72c7d1ea9a93b4a769d4  DATfiles/list_of_optics_run.dat
3bf27d1e057d76c525074c15b89b7f19d070f3e04a519af2d98f0d4e9c32e8c1  HMS_6p667GeV/config/oldfit.dat
aea3dcae3e47735d7ea57ea0c7079a6d6f6ab6708fc58057bf7270d11444c500  HMS_6p667GeV/config/rungroups_6p667_inputs.tsv
ca3820ce1c8df08d5b4fb4c6a189652478562c1fd751fc9a785f42459f8e3044  HMS_6p667GeV/06b_svd_fit/matrices/nps_hms_newfit_6p667_gmm_clean.dat
655e27c897060c8d751bfcce64290d8192757eccf74277aab8d7473dbf7ed5ed  HMS_6p667GeV/06a_fit_ntuple/root/Optics_666701_-1_fit_tree_gmm.root
72e7e904256ef053433584a33d0e86cf55651ec49b71899a6a4b7bf578689cb3  HMS_6p667GeV/06a_fit_ntuple/root/Optics_666702_-1_fit_tree_gmm.root
545c61a62a1fd12d7d0487ac6f1dce5d92491ac187946ae4c4b7a62e74b4134f  HMS_6p667GeV/06a_fit_ntuple/root/Optics_666703_-1_fit_tree_gmm.root
3f4be0f90cd116c1b1606a6932d99277a0127b421d25107f58319f5dbf199d4d  HMS_6p667GeV/06a_fit_ntuple/root/Optics_666704_-1_fit_tree_gmm.root
3e3c3ae5e15de566637e57b02c2b4ee3dc072312ebd9f5ffffa0f5e8a5103b5c  HMS_6p667GeV/06a_fit_ntuple/root/Optics_666705_-1_fit_tree_gmm.root
c1fb88cd18c5a856b9c9e9dfcec23cd205b40f425a85fb83a1846ba50eb95824  HMS_6p667GeV/06a_fit_ntuple/root/Optics_666706_-1_fit_tree_gmm.root
90e1fb240041f7bb418b1bf7eeab61af594e609ecf35e6a590c8fe50e3682538  HMS_6p667GeV/06a_fit_ntuple/root/Optics_666707_-1_fit_tree_gmm.root
```

## Step 1: reproduce the existing angular fit

The first local rerun used the existing fit ntuples, seed matrix, setting list,
ROOT `TDecompSVD`, and the original parser.  The equivalent command is:

```bash
root -l -b -q \
  'fit_opt_matrix_gmm.C+("local_repro",-1,200000,-1,"HMS_6p667GeV/06a_fit_ntuple/root","HMS_6p667GeV/06b_svd_fit_local_repro","HMS_6p667GeV/config/oldfit.dat","HMS_6p667GeV/config/rungroups_6p667_inputs.tsv","DATfiles/list_of_optics_run.dat")'
```

This selected exactly the same 166,013 events and reproduced the prior
per-setting and per-foil counts.  It did not produce a trustworthy byte-for-byte
matrix reproduction because the seed parser accepted a blank line as though it
were a coefficient row.

The preserved uncorrected matrix and QA files are under:

```text
HMS_6p667GeV/08_preliminary_conditioning/reproduction_check/uncorrected_local_rerun/
```

Its matrix checksum is:

```text
4bd0370e13cdf61b912407affb7e69c023f1a1fd4c2259546874abdb8d867801  nps_hms_newfit_local_repro.dat
```

## Step 2: make the seed-matrix parser strict

The return value of `sscanf` in `fit_opt_matrix_gmm.C` was not checked.  A blank
line after the first matrix separator could therefore append uninitialized
coefficients and exponents.  Which retained `xtar` row appeared depended on the
machine's memory contents.

The approved correction accepts a row only when all four coefficients and five
exponents are parsed (`parsed == 9`); blank or malformed lines are skipped.  No
other fitting rule was changed.

The corrected rerun command is:

```bash
root -l -b -q \
  'fit_opt_matrix_gmm.C+("local_clean_parser",-1,200000,-1,"HMS_6p667GeV/06a_fit_ntuple/root","HMS_6p667GeV/06b_svd_fit_local_clean","HMS_6p667GeV/config/oldfit.dat","HMS_6p667GeV/config/rungroups_6p667_inputs.tsv","DATfiles/list_of_optics_run.dat")'
```

The corrected parser finds 461 valid seed rows: 209 `xtar`-independent rows to
refit and 252 explicit-`xtar` rows retained unchanged.  The fit code adds the
constant term, giving 210 adjusted coefficients.  The prior ifarm output had
253 retained `xtar` rows because of the malformed-row ambiguity.

The corrected local matrix and QA files are archived under:

```text
HMS_6p667GeV/08_preliminary_conditioning/reproduction_check/strict_parser_local_rerun/
```

Its matrix checksum is:

```text
a04cddbdd5323f42b593934aad10f8d903e4932e159964844a33fe91a363eb02  nps_hms_newfit_local_clean_parser.dat
```

The corrected coefficients are not byte-identical to the prior ifarm matrix.
Their relative coefficient differences are about 2.3--3.0%, while predictions
on the common selected sample differ by only 0.00455 mrad RMS in `xptar`,
0.000781 cm RMS in `ytar`, and 0.00264 mrad RMS in `yptar`.  This ambiguity was
retained and reported rather than silently choosing one matrix.

## Step 3: define the design matrix and samples

The design matrix columns follow the `xtar`-independent terms in `oldfit.dat`,
preceded by a constant column.  The focal-plane variables are
`xfp/100`, `xpfp`, `yfp/100`, and `ypfp`.  Explicit-`xtar` terms are evaluated
with their seed coefficients and subtracted from the three angular-fit targets;
they are not refitted in this workflow.

An event is eligible when it:

- has `-15 < delta < 30`;
- is within 2.5 cm of a configured foil center;
- falls in one of the configured delta intervals; and
- is within 0.5 cm of one of the nine horizontal sieve-row centers.

Training events are selected in the original tree order, matching
`fit_opt_matrix_gmm.C`, with these caps:

- at most 1,000 events per foil/delta/sieve-row cell;
- at most 15,000 events per foil within each setting; and
- at most 200,000 events globally.

This selects 166,013 training events.  The 278,684 held-out events pass the same
eligibility cuts but were not selected because of those caps.  They come from
the same runs as the training events; this is a within-dataset validation, not
an independent replay or independent-run test.

The foil totals used for the per-foil plots are:

| nominal foil | training | held out |
|---:|---:|---:|
| -8 cm | 44,310 | 36,160 |
| -3 cm | 25,505 | 8,218 |
| 0 cm | 30,000 | 199,492 |
| +3 cm | 25,731 | 9,857 |
| +8 cm | 40,467 | 24,957 |

## Step 4: initial conditioning and reproduction comparison

For a safe rerun that does not overwrite the committed results, define a new
output directory:

```bash
CAMPAIGN=HMS_6p667GeV
META=DATfiles/list_of_optics_run.dat
OUT=HMS_6p667GeV/08_preliminary_conditioning_reproduced
mkdir -p "$OUT"
```

Then run:

```bash
python3 diagnostics/conditioning/preliminary_angular_conditioning.py \
  --campaign "$CAMPAIGN" \
  --metadata "$META" \
  --output "$OUT" \
  --file-id -1 \
  --nfit-max 200000 \
  --matrix-a HMS_6p667GeV/06b_svd_fit/matrices/nps_hms_newfit_6p667_gmm_clean.dat \
  --matrix-b HMS_6p667GeV/08_preliminary_conditioning/reproduction_check/strict_parser_local_rerun/matrices/nps_hms_newfit_local_clean_parser.dat
```

This produces the full-design singular-value table, raw and unit-column
condition numbers, leave-one-setting-out conditioning, and the saved-versus-
corrected matrix prediction comparison.

## Step 5: run every term count with the current XTX method

The full low-N ladder was run; no step-of-ten shortcut was needed:

```bash
python3 diagnostics/conditioning/angular_term_ladder.py \
  --campaign "$CAMPAIGN" \
  --metadata "$META" \
  --output "$OUT/term_ladder" \
  --file-id -1 \
  --nfit-max 200000 \
  --step 1
```

For each prefix `N=1..210`, this forms `G_N = X_N^T X_N` and solves with ROOT
`TDecompSVD`, matching the current production method.  It also evaluates the
unit-column Gram matrix for conditioning comparison.

ROOT's solve-specific cutoff is

```text
tau_N = ROOT tolerance * largest singular value at that N.
```

The plotted truncation margin is `log10(sigma_min / tau_N)`.  Zero therefore
means that the weakest direction has reached that particular solve's cutoff;
it is not a fixed threshold line reused for every `N`.

## Step 6: compare four numerical solve routes

```bash
python3 diagnostics/conditioning/angular_solver_comparison.py \
  --campaign "$CAMPAIGN" \
  --metadata "$META" \
  --saved-matrix HMS_6p667GeV/06b_svd_fit/matrices/nps_hms_newfit_6p667_gmm_clean.dat \
  --gram-ladder "$OUT/term_ladder/angular_term_ladder.tsv" \
  --output "$OUT/solver_comparison" \
  --file-id -1 \
  --nfit-max 200000
```

The four numerical routes are:

1. current unscaled `X^T X`, solved with ROOT `TDecompSVD`;
2. direct SVD/least-squares solution of unscaled `X`;
3. direct SVD/least-squares solution after scaling every column of `X` to unit
   norm, followed by conversion back to the original coefficients; and
4. ROOT SVD of the normal equations formed from the unit-column matrix.

The comparison records condition numbers, retained numerical rank, training
residuals, held-out residuals, coefficient differences, prediction differences,
and the preliminary geometry-derived vertex stress metric.

## Step 7: make training and held-out residual ladders for each foil

```bash
python3 diagnostics/conditioning/angular_per_foil_residual_ladder.py \
  --campaign "$CAMPAIGN" \
  --metadata "$META" \
  --current-coefficients "$OUT/term_ladder/angular_term_ladder_coefficients.npz" \
  --output "$OUT/per_foil_residual_ladder" \
  --reference-n 126 \
  --file-id -1 \
  --nfit-max 200000
```

This independently refits the rescaled direct-X solution at every `N` and
compares it with the saved current-XTX coefficient ladder.  For every foil and
each of `xptar`, `ytar`, and `yptar`, the plots show training and held-out RMS.

The black dashed line is the minimum held-out RMS for `N >= 35`.  The red
dotted marker is the first run of five consecutive `N` values more than 5%
above that minimum.  It is a descriptive marker, not a recommended truncation
or production cutoff.  The lower panel is relative to the complete fifth-order
solution at `N=126`.

## Step 8: review and presentation outputs

The numerical results were summarized in:

```text
HMS_6p667GeV/08_preliminary_conditioning/PRELIMINARY_FINDINGS.md
HMS_6p667GeV/08_preliminary_conditioning/COLLABORATION_SUMMARY.md
HMS_6p667GeV/08_preliminary_conditioning/ANALYSIS_BUNDLE_README.md
HMS_6p667GeV/08_preliminary_conditioning/HMS_6p667_angular_stability_preliminary.pptx
```

The slide deck uses committed plots and contains plain-language speaker notes.
It reports conditioning and residual behavior without recommending that a
particular matrix be committed.

## Checks that should agree on another device

At minimum, a successful rerun should reproduce these bookkeeping and
conditioning landmarks, allowing for ordinary floating-point text formatting:

- 461 valid seed rows;
- 209 nonconstant `xtar`-independent rows plus the added constant, for 210 fit
  coefficients;
- 252 explicit-`xtar` rows retained unchanged;
- 166,013 training events and 278,684 held-out events;
- first current-XTX truncated direction at `N=54`;
- retained current-XTX rank `111/210` at `N=210`;
- retained rank remains 111 from `N=176` through `N=210`;
- full-order current-XTX condition number approximately `5.732e26`;
- full-order rescaled-direct-X condition number approximately `3.468e5`;
- the direct unscaled-X, rescaled-direct-X, and rescaled-XTX predictions agree
  to floating-point precision on the training problem; and
- the strongest high-N held-out degradation is concentrated in the +8 cm
  subset, with the final 12.495-degree +/-8 setting a major contributor.

Do not interpret the rank statement as saying that every coefficient added
after a particular `N` is zero.  An SVD direction is a mixture of named
coefficients, and some later terms add independent information after the first
truncation.  The supported statement is that truncation begins at `N=54`,
remains active thereafter, and the number of retained independent directions
stops increasing at 111 by `N=176`.

## Retained limitations

- The prior ifarm matrix contains the parser-dependent extra retained-`xtar`
  row; it was preserved rather than overwritten.
- The held-out set shares runs and selection cuts with the training set.
- Per-foil curves aggregate settings and can hide a setting-specific effect.
- The vertex comparison is a preliminary transport/geometry stress test, not a
  full replay validation.
- No production matrix choice was made in this study.
