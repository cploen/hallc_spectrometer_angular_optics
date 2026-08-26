# Preliminary 6.667 GeV angular-optics conditioning and stability

Date: 2026-08-26

## Bottom line

The current 210-term angular fit is strongly ill-conditioned in its production form, which solves the unscaled normal-equation matrix `X^T X`. ROOT retains only 111 of 210 modes at `N=210`. The first discarded mode appears at `N=54`.

The low-`N` ladder separates three useful regions:

- `N <= 39`: all modes retained with a meaningful margin. `N=35` is the last complete polynomial degree; `N=39` is an incomplete quartic prefix but gives a material xptar improvement.
- `N=40..53`: all modes are still retained, but the weakest mode is only about 0.13--0.17 decades (roughly 1.36--1.48 times) above that solve's own truncation threshold. This is not a comfortable numerical margin.
- `N >= 54`: the production Gram solve truncates modes. At `N=210`, 99 modes are discarded.

Column scaling changes the numerical picture substantially: the unit-column Gram matrix retains all 210 modes, with its weakest mode 4.57 decades above its own threshold. This points to basis scaling/preconditioning, rather than merely adding or removing data settings, as the first numerical issue to investigate.

These are preliminary training-sample results, not yet a recommendation to replace the production matrix.

## Reproduction and retained ambiguity

The clean local rerun selected exactly the same 166,013 events and reproduced every per-setting and per-foil event count from the prior ifarm run.

The seed-matrix parser had accepted a blank line as an uninitialized row because the return value of `sscanf` was not checked. The prior ifarm matrix therefore contains one machine-dependent extra xtar row. The corrected parser rejects malformed rows, yielding 461 valid seed rows and 252 retained xtar rows rather than 253. The prior output and the first uncorrected local rerun were preserved; neither was overwritten.

After correction, valid coefficient vectors differ from the prior ifarm output by 2.3--3.0% in relative L2 norm, but their predictions on the common selected sample differ only by:

| quantity | RMS prediction difference | maximum difference |
|---|---:|---:|
| xptar | 0.00455 mrad | 0.156 mrad |
| ytar | 0.000781 cm | 0.0336 cm |
| yptar | 0.00264 mrad | 0.204 mrad |

Changing only how `X^T X` and the right-hand side are accumulated, while retaining ROOT's solver, changes coefficient vectors by about 2.3--4.6% but predictions by only 0.0024--0.0050 mrad for the angular quantities and 0.0015 cm for ytar. This is direct evidence that coefficients are numerically unstable even when in-sample predictions remain close.

## Complete term ladder

All 210 prefixes were evaluated; no step-10 approximation was needed. Every rung uses ROOT `TDecompSVD`, matching the production solver and its retention rule

`sigma_i > tolerance * sigma_max(N)`.

The plotted threshold distance is solve-dependent:

`log10(sigma_min(N)) - log10(tolerance * sigma_max(N))`.

Thus zero means equality with that rung's own threshold. The horizontal zero reference is not a global truncation threshold.

| N | basis endpoint | rank | production margin (decades) | scaled margin (decades) | xptar RMS (mrad) | yptar RMS (mrad) | ytar RMS (cm) |
|---:|---|---:|---:|---:|---:|---:|---:|
| 5 | complete linear | 5/5 | +11.047 | +14.101 | 8.316 | 3.184 | 1.168 |
| 15 | complete quadratic | 15/15 | +6.485 | +12.180 | 3.122 | 1.290 | 0.2849 |
| 35 | complete cubic | 35/35 | +1.977 | +10.252 | 2.359 | 1.154 | 0.1714 |
| 39 | partial quartic through `13000` | 39/39 | +1.976 | +9.357 | 2.154 | 1.153 | 0.1702 |
| 40 | partial quartic through `04000` | 40/40 | +0.170 | +8.730 | 2.154 | 1.153 | 0.1701 |
| 53 | last untruncated production prefix | 53/53 | +0.133 | +8.554 | 2.083 | 1.144 | 0.1614 |
| 54 | first truncated prefix | 53/54 | -0.876 | +8.542 | 2.084 | 1.144 | 0.1614 |
| 70 | complete quartic | 63/70 | -2.154 | +8.349 | 2.061 | 1.139 | 0.1536 |
| 126 | complete quintic | 90/126 | -6.807 | +6.448 | 1.992 | 1.128 | 0.1492 |
| 210 | complete sextic | 111/210 | -11.105 | +4.573 | 1.978 | 1.126 | 0.1478 |

The largest xptar prediction change after `N=35` occurs when term `22000` enters at `N=38` (0.899 mrad RMS). This explains why the numerically healthy `N=39` prefix improves xptar substantially, but it remains an order-dependent, incomplete-degree basis and should not be treated as a final model choice without reviewing the term ordering and physics symmetry.

Using `N=210` as the comparison target does not make it ground truth: its production solve already discards 99 modes. Residual improvements beyond the lower-order models are real on the fit sample, but must be weighed against this truncation and checked out of sample.

## Conditioning diagnostics

For the direct design matrix `X` at `N=210`:

- raw condition number: `3.48e13`;
- unit-column condition number: `3.47e5`;
- raw column-norm range: factor `3.42e10`;
- maximum absolute unit-column correlation: `0.997866`.

Leave-one-setting-out unit-column condition numbers range only from about 332,000 to 365,000. No single run group is responsible for the ill-conditioning.

## Preliminary interpretation

1. Use `N=35` as the conservative complete-degree reference model for immediate discussion.
2. Also show `N=39` as a numerically healthy exploratory prefix because it captures a large xptar improvement; label the incomplete-quartic caveat explicitly.
3. Do not call `N=53` robust merely because ROOT retains all modes: it is only 0.133 decades above its own threshold.
4. Treat `N=54` as the onset of production-solver truncation for this ordering and selected sample.
5. Next, compare scaled/preconditioned solves and validate candidate term sets by held-out setting or foil. That is beyond this preliminary one-hour pass.

## Files

- Full ladder table: `term_ladder/angular_term_ladder.tsv`
- Coefficient archive: `term_ladder/angular_term_ladder_coefficients.npz`
- Solve-specific threshold-distance plots: `term_ladder/truncation_margin_vs_N.png` and `term_ladder/truncation_margin_lowN.png`
- Conditioning plots: `term_ladder/condition_vs_N.png` and `term_ladder/condition_lowN.png`
- Stability and fit plots: `term_ladder/prediction_step_vs_N.png` and `term_ladder/residual_rms_vs_N.png`
- Design-matrix/leave-one-setting-out report: `preliminary_angular_conditioning.txt`

## Known limitation not silently changed

The production macro's post-fit xptar QA section appears to add the xtar correction twice internally and then fills the displayed histogram without that correction. This does not affect the written matrix, and it was not silently changed for this preliminary analysis. The ladder residuals were computed directly from the corrected targets and solved coefficients instead.
