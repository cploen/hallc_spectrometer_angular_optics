# Separate x-target predictor: implementation checkpoint

This is the separate predictor experiment described in
`DIRECT_XTAR_HANDOFF_20261111.md`. The superseded conditional-coefficient runner
is not used. This checkpoint implements the numerical core; it is **not a
completed real-data fit**, a production feature definition, or an installed
replay change.

## Numerical contract

`elastic_net.tune`, `beam_search.Folds` and `beam_search.fit_all` now accept
explicit keyword-only `targets` and `scales`. Their existing three-target
interfaces and default units remain unchanged. A single response has shape
`(n, 1)`; a mismatched shape or reporting-scale count is rejected.

`xtar_predictor.fit_training_design` accepts an already-audited measured-input
design, a single native `xtarT` response in **cm**, existing training folds and
scoring-cell labels. It uses the actual single response, never duplicated
angular targets. It requires beam width **30**. There is a free intercept in
this separate predictor; no transport-function intercept is altered.

Within each training fold, scaling comes only from that fold's fitting rows.
For centered fits, `Z=(X-mean_X)/scale_X` and `t=(y-mean_y)/scale_y`.
The unpenalized solve is `lstsq(Z_selected, t, rcond)`; beam candidates use the
established orthogonal QR compression of `[Z,t]` and then direct least squares
on selected columns of R. No normal-equation solve is introduced. EN retains
its existing Gram precomputation for penalized coordinate descent; this is not
a normal-equation SVD/refit.

Beam CV scores use pooled equal-cell MSE with `10 mm/cm` scaling. Elastic-net
CV tables remain in native cm squared, as the established tuner reports native
response MSE. The reporting scale is never applied to exported coefficients.
Final native slopes are `beta*scale_y/scale_X`; the intercept is
`mean_y - mean_X @ slopes`.

The new exporter creates a distinct JSON file exclusively and requires ordered
feature names, expressions, source descriptions and units, plus provenance and
validity information. Its output is cm. Its algebraic export check is separate
from the end-to-end test of predictions on event rows. It does not write any
transport file. Supplying a feature description does not certify that the
feature is measured or free of truth leakage: the production adapter and audit
must establish that before fitting.

## Evidence already available

- ifarm pasted hashes: both `oldfit.dat` candidates equal
  `b2036a935db3f17189fcf489d1017da030a8bb58de7d5a091384ef3d3b60c92c`.
  Equality is established; historical original-2018 provenance is not yet
  certified.
- Completed zero-row-removed refit export matches the handoff hash
  `f715879845f777f2da2aa62b7e66756907929bdf4b7b40809ac55370156f8801`.
- ifarm manifest records 118,850 fit events: 95,080 core and 23,770 inner
  shoulder; 23,770 per physical foil. The 80/20 mixture is distinct from the
  reserved allocation.
- ifarm run-1338 tree contains 21,058 entries and `xtarT`, `reactxcalc`,
  `reactycalc`, but no exported raster/BPM branches. A filename containing
  `_gmm` does not establish the method used. No GMM fitting is performed by
  this implementation or the read-only preflight.
- Repository commit `d60afd9` (July 15, 2026) already binds `H.react.y` and
  exports it as `reactycalc`; the variable was not introduced by this experiment.
  Whether it is an appropriate measured predictor input remains an installed
  HCANA dependency question.

## Read-only real-data checkpoint

`xtar_direct_preflight.py` uses the existing ifarm checkout. It prints installed
HCANA excerpts and hashes, matrix headers/hashes, relevant configuration,
exporter/allocation excerpts, Python versions, and exact training membership,
geometry, saved-fold and coverage validation through the existing loader.
It reads training data only and does not call a selector, fitter, replay or
matrix writer. It disables Python bytecode writes. It reports missing evidence
or validation errors explicitly. Excerpts report truncation; absent text is
not evidence of absent behavior. It can take time because it verifies frozen
ROOT input hashes.

Do not interpret a successful preflight as complete HCANA closure. Historical
2018 replay closure is still unresolved and is not a prerequisite reinstated
by this work.

## Still required before a production fit/comparison

1. Review installed transformations, angle signs, beam/vertex dependencies,
   mispointing and external configuration. Obtain complete sections if an
   excerpt is truncated. Establish original-matrix historical provenance.
2. Define a compact measured-input basis, units, intercept convention and
   validity. No foil/hole/truth-derived input features. Preserve exact saved
   folds, memberships, masks and exclusions.
3. Supply the verified data adapter, then run training-only selection and SVD.
   Check real-event JSON prediction identity and all source/input hashes.
4. Evaluate original 2018, completed refit and completed refit plus predictor
   on identical reserved IDs. Preserve original delta coordinates; use edges
   `-15,-12,-10,-8,-5,0,5,8,10,12,14,15` for all four variables. Report signed
   biases, RMS including bias, signed RMS differences, counts, core and
   supported noncore, unsupported/blocked counts, sparse/blank cells.
5. Treat reused holdouts as development evidence. Add a run/setting stress
   check; retain the offline-iteration assumptions and closure limitation.
   Prepare a factual review note; no messages to other people are authorized.
6. Review before any replay installation. Matched baseline/candidate HCANA
   replay is eventual validation; no installation is included here.
