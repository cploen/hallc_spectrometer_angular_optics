# Separate x-target predictor demonstration

This implements the separate predictor experiment described in
`DIRECT_XTAR_HANDOFF_20261111.md`. It does not use the superseded
conditional-coefficient runner. All existing transport coefficients stay fixed.
No real-data fit has yet been executed with this code at publication time.

## Numerical contract and checks

`elastic_net.tune`, `beam_search.Folds` and `beam_search.fit_all` accept explicit
keyword-only `targets` and `scales`; old three-target defaults remain unchanged.
The new path has one response of shape `(n,1)`, never duplicated angular targets.
Native positions are cm; beam scores use the explicit `10 mm/cm` reporting
scale. Elastic-net CV tables retain native cm-squared scores. Exported native
coefficients are never multiplied by an angular mrad scale.

Within each fold, fitting rows alone determine scaling:
`Z=(X-mean_X)/scale_X`, `t=(y-mean_y)/scale_y`. Selected-term unpenalized fits use
`lstsq(Z_selected,t,rcond)`. Beam candidates use orthogonal QR compression of
`[Z,t]`, then direct least squares on the selected R columns. No normal-equation
SVD/refit is introduced. EN's existing Gram precomputation remains a penalized
coordinate-descent optimization, not a normal-equation SVD solve.

Native correction slopes are `beta*scale_y/scale_X`; their separate intercept is
`mean_y-mean_X @ slopes`. `xtar_predictor.py` explicitly subtracts the measured
anchor before fitting and adds it back for predictions. Its JSON format records
features, expressions, units, anchor, provenance and training validity. An
anchored model rejects prediction calls without its measured anchor. Exports
are new files only, separate from HCANA matrices. Feature descriptions alone do
not certify absence of truth leakage; the adapter and tests check construction.

Validation at publication: 26 targeted tests passed: 11 existing elastic tests,
4 existing beam tests, 3 existing conditional-extension regression tests on
synthetic fixtures (no stopped production run), 3 new single-target numerical
checks, 3 new anchored-pipeline checks, and 2 new diagnostic checks. Coverage
includes direct scaled-X vs QR/SVD equivalence, mm scoring, rank rejection,
export identity, truth-feature exclusion, mispointing signs, real-loader
synthetic integration, unchanged matrix bytes, exact reserved IDs, and extended
bin edges/RMS signs/blocked counts. A synthetic plot smoke check was rendered
and inspected; it is not an experimental result.

## Frozen reference evidence

Both ifarm `oldfit.dat` candidates match
`b2036a935db3f17189fcf489d1017da030a8bb58de7d5a091384ef3d3b60c92c`.
Their header says HMS tune 2018, created 03/02/18. These facts do not establish
the complete historical chain of custody; the original-2018 provenance remains
qualified. The reference bytes are never overwritten.

The exact completed-refit export is
`f715879845f777f2da2aa62b7e66756907929bdf4b7b40809ac55370156f8801`.
It is the entire-00000-row-removed zero-offset replay experiment, not a fit
that constrained every intercept during optimization. Its remaining
coefficients stay unchanged here.

The accepted training allocation has 118,850 rows, including 95,080 core and
23,770 inner shoulder. This 80/20 population mixture is not a train/test split.
Saved folds have 39,624, 39,576 and 39,650 rows. No GMM fit, relabeling or
reallocation occurs. A legacy `_gmm` filename is not a method declaration.

`reactycalc` was already exported from `H.react.y` in repository commit
`d60afd9` (July 15, 2026); it was not introduced for this experiment.

## Proof-of-principle runner (scope clarified by the user)

Historical replay settings are unavailable. This is a demonstration and a
possible starting point for lower-momentum NPS settings. Recovering that
historical configuration is **not a gate**. Its unknown settings are not inferred.

`run_xtar_predictor.py` now implements that demonstration with a declared
raw-coordinate comparison convention. It evaluates each frozen matrix once
on saved focal-plane measurements, with positions converted cm to m, optional
focus/offset/slope corrections set to zero for this comparison, and no external
target offsets. Matrix constants are evaluated exactly as stored. Oldfit and
the completed refit share the saved `xtar` input; the candidate substitutes the
separate predictor once. There is no vertex update or iterative HCANA replay.
No historical replay equivalence is claimed. Labels retain the unresolved
original-2018 provenance qualification instead of silently certifying it.

The compact starting basis has 20 columns: all degree-1 and degree-2 monomials
of `(xfp_m, xpfp_rad, yfp_m, ypfp_rad, beam_y_m)`, each multiplied by
`cos(angle_rad)`. Beam y is `reactycalc/100`. The fixed anchor is
`-reactycalc-x_mispointing` in cm, using the same capped-angle x-pointing formula
as the geometric target. The single learned response is `xtarT-anchor` in cm.
A free intercept belongs to this separate correction model. No foil, hole,
true slope, true target position, reconstructed z or saved xtar is a predictor
feature. Truth labels define responses and diagnostic grouping only. Feature
selection still needs real-data scrutiny; the model is not guaranteed to
identify a unique inverse or improve experimental reconstruction.

`--fit` opts into fitting; default invocation is read-only input/design/fold
validation. `--evaluate` then compares the three models on exact reserved IDs,
using the requested extended delta edges and four variables. Core and
supported noncore are separate; support is determined from training physical
foil/original-delta/hole coverage. Blocked/unsupported/outside-range counts are
explicit. All angle residuals are mrad and both position residuals are mm.
Negative RMS differences indicate improvement. Original delta coordinates and
membership never change with the candidate. Eight wide figures, the complete
cell CSV and pooled summary are written, along with local event-level evidence.
Event-level outputs are retained on ifarm; they are not uploaded by this code.

The user's selected validation design uses the original saved folds, each
containing data from all three accepted angles. The runner verifies this
property and writes `angle_folds.tsv`. It does not perform a leave-angle-out
stress check or change event allocations. Only 583 accepted training electrons
are at about 30 degrees. The user explains that good electron statistics are
limited by high pi0 contamination there and that the transport matrix remains
useful for positrons. Sparse electron coverage is not evidence that transport
is unusable for positrons; this electron-trained predictor has not separately
been validated on positron data.

Output is exclusively a new `06h_xtar_predictor/<name>` directory. Existing
output names are refused. The stopped `06g_xtar_refit` experiment and all
existing matrix files are untouched. The runner records hashes, copies frozen
reference matrices and code into its new evidence directory, and rechecks
inputs after completion. It does not install anything into replay.

## HCANA source audit evidence

Read-only ifarm preflight verified 118,850 exact training rows, all 44 input
hashes, saved geometric equations, and the three saved folds. It reported a
missing optional `config/comparison.json`; this does not invalidate membership.

The fetched JeffersonLab/hcana master sources were compared with the installed
source hashes supplied by the preflight. `THcExtTarCor.cxx` is byte-identical
(SHA256 `13e999f2bb9a76fca7b5e2948797e124698cd4050073db7c6fc84f39cfe4555a`), as is
`THcReactionPoint.cxx` (SHA256
`57842f0f7096f0e7b4e4abda16602072d641d07bb6b8b88a7ad2bf20144884b4`).
`THcHallCSpectrometer.cxx` differs by whole-file hash, but the supplied installed
polynomial/input/output excerpts and upstream equations agree on the checked
coordinate chain. This does not establish equality of every version detail.

Checked formulas:

- ReactionPoint sets `v=(beam_org.x,beam_org.y,ztarg)`. Thus reacty is the beam
  position coordinate in this implementation, rather than its reconstructed z.
- ExtTarCor uses `x_tg=-vertex[1]-pointing_off[0]` in cm, followed by
  `xtar_new=x_tg-xptar*ztarg*cos(theta)`. The subtraction signs agree with the
  saved geometric target check. Golden-track target position is stored in cm.
- Matrix input 5 is `xtar/100`. Focal-plane positions are cm -> m; slopes rad.
  Column order is xptar, ytar, yptar, delta; delta is not an xtar output.
- xptar adds `fPhiOffset`; yptar adds `fThetaOffset` (names are easy to confuse).
  This demonstration applies neither external offset.
- Upstream `InitializeReconstruction` zeroes the focus/detector/angular terms;
  its focal-plane header parser currently contains commented Fortran lines,
  not active assignments. The formula-based HMS x mispointing is
  `0.1*(2.37-.086*d+.0012*d*d)` with `d=min(abs(angle_deg),50)`. These observations
  support the declared convention, not a certification of historical settings.

Sources:
https://github.com/JeffersonLab/hcana/blob/master/src/THcHallCSpectrometer.cxx
https://github.com/JeffersonLab/hcana/blob/master/src/THcExtTarCor.cxx
https://github.com/JeffersonLab/hcana/blob/master/src/THcReactionPoint.cxx

No analyzer defect has been demonstrated. A difference between saved xtar and
geometric xtarT alone is not evidence of a defect. No messages have been sent
to Mark, Charles, Josh or anyone else.
