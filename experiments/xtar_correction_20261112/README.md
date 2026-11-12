# Conservative additive xtar correction — 12 November 2026

The approved experiment preserves the complete 2018 refit transport matrix and
learns a separate correction to its iterative xtar coordinate. No replay is
modified. `oldfit.dat` stays original 2018. Both matrices retain all their own
stored constant terms, with zero external additions. Width (standard deviation)
is evaluated per foil/delta cell; signed means are shown separately.

## Human command

From this experiment's source directory on ifarm:

```bash
A=/home/cploen/nps_link/ML_HMS_OPTICS_JOINT_SURFACE
C="$A/HMS_LowMom_LowDelta_2018_Accepted"
python3 -u xtar_correction.py --checkout "$A" --campaign "$C" \
  --out /home/cploen/nps_link/xtar_correction_2018_20261112/results \
  --threads 4
```

Output must not already exist. Existing fits, event allocations, and results
are never overwritten. Source code is snapshotted in the new output.

## Frozen protocol

- Existing exact training allocation: 80% core / 20% inner shoulder, not a
  training/test split. Original exclusions, masks and IDs are unchanged.
- xtarT is independently checked against vertex-to-hole geometry:
  `xv=-reactycalc-xmis`, `zv=zfoil*cos(angle)`,
  `xptarT=(xsT-xv)/(168-zv)`, `xtarT=xv-xptarT*zv` (cm and radians).
- Baseline uses the established fixed-saved-vertex afterburner, five updates
  maximum, stopping at a 2 mrad xptar change. Zero correction is checked against
  that implementation on every training event. This does not claim historical
  replay closure or recompute the reaction vertex.
- Correction: `xtar_proposal=anchor-xptar*saved_react_z*cos(angle)+g`.
  `g` is recomputed from fixed measured inputs, never accumulated across updates.
- Basis: 20 degree-one/two monomials of focal-plane coordinates and beam y,
  multiplied by cos(angle). No foil/hole/run identity or truth inputs, and
  **no free correction constant**. Matrix positions are converted cm to m.
- Three folds keep complete runs together. EN screens the basis. Beam width
  **30**, at most **six nonzero terms**, ridge strength 0.03 retained through
  final solve. Ridge uses augmented scaled X with direct SVD; no normal-equation
  solve. EN uses the actual penalized coefficients, not an unpenalized refit.
- Beam candidates use amplitudes 0.1, 0.25, 0.5, 1; zero always remains a candidate.
  Reject corrections exceeding 0.2 cm on validation/reserved events rather than
  clipping them. Reject rank deficiency, scaled condition above 10,000, or new
  iteration-limit failures.
- No preferential delta weights. Train EN/ridge on correction-to-xtarT error.
  Select using event-count-pooled within-cell xtar residual variance, requiring
  no increase in xtar MSE. Prefer fewer terms within 0.5% variance of best;
  require at least 0.1% width improvement to retain a correction.
- **2% maximum downstream SD increase** overall, separately for core and
  shoulder, and in each foil/delta cell of each population with N>=100. Apply
  in each validation fold and pooled validation, and again after final solve
  on the reserved development sample. Absolute mean bias may not increase by
  more than 2% of the corresponding baseline SD. Thus pure alignment changes
  are neither credited as width gains nor allowed to hide large displacements.
- Delta edges: -15,-12,-10,-8,-5,0,5,8,10,12,14,15 percent. Original saved delta
  defines all bins. Cells with N<100 are displayed and explicitly uncertified.
- Freeze EN and beam coefficients before opening reserved data. Beam is the
  predeclared final proposal. If its reserved checks or xtar improvement fail,
  the accepted result is **zero correction**; do not retune using that sample.
  Previously reviewed reserved events remain development evidence, not a fresh
  independent test. Original 2018 comparison does not certify historic offsets.

## Review

`plots/xtar_refit_summary.png`: four-panel basis/width/conditioning/reserved
summary analogous to the delta-refit figure. Other pages show core and shoulder
SD changes, absolute widths, and separate signed mean changes for all four
variables. Comparisons include refit/original, EN/refit, beam/refit, beam/EN,
and final/original. Exact values/counts are in `reserved_metrics.tsv`.

`matrices/final_transport.dat` is byte-identical to the completed 2018 refit.
`final_correction.json` is the separate xtar model; the fourth transport column
remains **delta**. The manifest declares acceptance or fallback, hashes, settings,
and limitations. Nothing is installed into hcana or replay by this command.
