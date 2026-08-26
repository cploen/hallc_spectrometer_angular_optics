# Preliminary numerical stability study of the 6.667 GeV HMS angular optics

## Bottom line

The present HMS angular-optics matrix gives reasonable reconstruction for most of the calibration events, but the way it is calculated is not numerically robust. At high order, the calibration data cannot cleanly distinguish many combinations of matrix coefficients. ROOT therefore removes those combinations using a cutoff determined by floating-point precision, not by an optics or physics criterion.

This accidental removal is doing two things at once:

1. It makes individual high-order coefficients sensitive to software version, rounding, and accumulation order.
2. It also acts as an accidental regularizer: it suppresses some high-order behavior that would otherwise produce large excursions near the edge of the sampled acceptance.

Simply rescaling the full 210-coefficient problem and retaining every coefficient combination is therefore **not** yet a safe production replacement. It improves the numerical formulation and the calibration-sample residuals, but the complete sixth-order solution develops rare, very large angular excursions on events excluded by the fit-sample caps.

The most promising preliminary path is to use a rescaled direct solve together with an explicit, validated choice of model complexity. The complete fifth-order basis at `N=126` is a useful first candidate to compare against the current production-style `N=210` solution. This is not yet a final matrix recommendation.

## What the matrix is doing

The transport matrix converts focal-plane measurements into target quantities needed for HMS reconstruction, including

- the in-plane target angle, `xptar`;
- the out-of-plane target angle, `yptar`;
- the transverse target position, `ytar`.

The combination of `ytar`, `yptar`, the HMS central angle, and beam position is used to reconstruct the event vertex along the beam direction. Therefore, a matrix can look acceptable in one coefficient or one angular residual while still biasing or broadening the reconstructed vertex.

## Why there are 210 fitted combinations but more than 461 matrix entries

The clean seed matrix contains 461 valid rows:

- 209 rows have no explicit `xtar` dependence and are refitted;
- 252 rows contain explicit `xtar` dependence and are carried forward unchanged in this fit.

The fit also introduces a constant row, so it adjusts `209 + 1 = 210` coefficients for each reconstructed quantity. The clean output matrix therefore contains `210 + 252 = 462` rows.

When this report says “210 directions” or “210 modes,” it means 210 independent combinations that could be formed from the 210 adjustable coefficients. It does not mean that the entire transport matrix has only 210 rows.

## What ROOT is doing in the present solve

Some polynomial columns in the fit are much larger than others, and some produce nearly the same event-by-event pattern. This is analogous to trying to determine many adjustment knobs when turning one knob up and another down gives almost the same detector response.

SVD rewrites the 210 coefficient knobs as 210 independent combinations, ordered from easiest to hardest for the data to distinguish. ROOT compares each combination with a numerical cutoff. Combinations below that cutoff are set to zero rather than fitted.

At `N=54`, the present unscaled `X^T X` calculation suppresses its first coefficient combination. At `N=210`, it retains 111 combinations and suppresses 99.

This does **not** mean that ROOT deletes 99 named polynomial rows from the matrix. Each suppressed direction is generally a mixture of many rows.

The phrase “not an unconstrained 210-parameter determination” means that the output file contains 210 adjusted coefficient values, but the data/solver combination did not independently determine all 210 possible combinations. ROOT forced 99 weak combinations to zero internally. In practical language, the solution behaves more like a 111-direction fit, with the exact 111 directions chosen by the numerical cutoff.

## What rescaling changes

Rescaling temporarily divides each column of the fit matrix by its overall size. This is like expressing all adjustment knobs in comparable units before solving. Afterward, the coefficients are converted back to the standard HMS transport-matrix convention.

Rescaling does not change the polynomial model or the predicted quantity in exact arithmetic. It changes how accurately a finite-precision computer can solve the problem.

Solving `X` directly is also preferable to first forming `X^T X`. Forming `X^T X` magnifies the difference between strong and weak directions: approximately speaking, it squares the condition number.

## Condition number in plain language

The condition number is a warning indicator for how strongly rounding, small input changes, or data perturbations can be amplified in the fitted coefficients. A value near one is easy to solve. Larger values are worse. It is a worst-case numerical diagnostic, not a direct prediction of the physics error.

For the full `N=210` angular problem:

| numerical formulation | condition number | combinations retained |
|---|---:|---:|
| current unscaled `X^T X` | `5.7e26` | 111/210 |
| direct solve of unscaled `X` | `3.5e13` | 210/210 |
| direct solve of rescaled `X` | `3.5e5` | 210/210 |
| solve of rescaled `X^T X` | `1.2e11` | 210/210 |

The direct unscaled-`X`, direct rescaled-`X`, and rescaled-`X^T X` solutions agree essentially exactly for this dataset. Thus, rescaling is not inventing a different least-squares answer. It allows the computer to recover the full least-squares answer without losing weak directions. The current unscaled-`X^T X` result differs because its weak directions have already fallen below ROOT's cutoff.

## Tangible differences on the selected calibration sample

| solution | xptar residual RMS | ytar residual RMS | yptar residual RMS |
|---|---:|---:|---:|
| current clean `X^T X` matrix | 1.978 mrad | 0.1478 cm | 1.126 mrad |
| rescaled direct `X` | 1.906 mrad | 0.1449 cm | 1.117 mrad |

On these same selected events, the current and rescaled solutions differ in their predictions by approximately

- 0.529 mrad RMS in `xptar`;
- 0.0291 cm RMS in `ytar`;
- 0.141 mrad RMS in `yptar`.

The rescaled solution fits the calibration sample somewhat better. That alone does not establish that it will behave better on production data.

## Behavior on events excluded by the fit-sample caps

There are 278,684 otherwise eligible events in the fit ntuples that were not used because of the per-cell and per-foil caps. These are useful as an immediate stress sample, although they are not statistically independent runs and are strongly dominated by one setting.

For 99% of these events, the absolute difference between the current and full rescaled solutions is below approximately

- 1.52 mrad in `xptar`;
- 0.079 cm in `ytar`;
- 0.382 mrad in `yptar`.

However, the full rescaled `N=210` solution has rare, extremely large angular excursions, concentrated in the final `theta=12.495 deg, foil=+/-8 cm` setting. Those tails make its overall held-out angular residuals substantially worse:

| solution | xptar residual RMS | ytar residual RMS | yptar residual RMS |
|---|---:|---:|---:|
| current clean `X^T X` matrix | 1.874 mrad | 0.1455 cm | 1.069 mrad |
| full rescaled direct `X` | 3.049 mrad | 0.5141 cm | 3.838 mrad |

This is consistent with high-order overfitting: the full rescaled solution follows the selected calibration sample more closely but is not controlled at sparsely represented edges of the sample.

## Consequences for vertex reconstruction

The errors in `ytar` and `yptar` are correlated, so they partly cancel when they are combined to reconstruct the beam-direction vertex. On the cap-excluded stress sample:

| solution | vertex-z RMS | median absolute error | 99th-percentile absolute error |
|---|---:|---:|---:|
| current clean matrix | 0.727 cm | 0.439 cm | 1.715 cm |
| full rescaled direct `X` | 0.692 cm | 0.436 cm | 1.710 cm |

The median difference in reconstructed vertex between the two matrices is about 0.070 cm, and 99% differ by less than about 0.367 cm. Both matrices still contain rare extreme vertex events that require inspection.

For this preliminary calculation, the event beam intercept was inferred from the known foil position and generated target coordinates stored in the fit ntuple. This isolates the transport-matrix contribution but is not a substitute for replay-level `H.react.z` validation with the production beam information.

This means the full rescaled solution is not simply “bad”: its vertex reconstruction is comparable and slightly better by RMS in this preliminary stress test. Nevertheless, its separate angular excursions make it unsafe to adopt without additional controls, because those angles are also physics outputs and enter other reconstruction quantities.

## What the term-count scan says

The present and rescaled methods are essentially equivalent through the low-order region where no directions are suppressed. The methods begin to differ when the current `X^T X` calculation starts dropping weak combinations.

Important complete-degree points are:

- `N=35`: complete cubic basis; numerically comfortable, but vertex resolution is not yet competitive.
- `N=70`: complete quartic basis; improved, but still weaker vertex performance.
- `N=126`: complete quintic basis; the rescaled direct solution gives about 0.654 cm vertex RMS on the stress sample without the severe sixth-order angular blow-up.
- `N=210`: complete sixth-order basis; the full rescaled solution develops rare angular excursions, while the current solution suppresses 99 directions.

Partial sixth-order prefixes around `N=140-160` also give approximately 0.65 cm vertex RMS in this preliminary sample, but they depend on the arbitrary ordering of terms and are less natural model definitions than the complete `N=126` basis.

## Decision before committing production data

The current evidence does not support committing all production analysis to either of these without one more validation step:

1. **Do not replace the current matrix with the full rescaled N=210 matrix.** Its improved conditioning and training residual do not compensate for its rare high-order angular excursions.
2. **Do not interpret the current N=210 coefficient list as uniquely measured.** Its apparent stability comes partly from ROOT suppressing 99 weak combinations at a machine-precision cutoff.
3. **Use a direct, rescaled solve for future candidate generation.** It is the cleanest numerical formulation and makes any regularization or term removal explicit rather than accidental.
4. **Compare at least the current clean N=210 baseline and a rescaled complete N=126 candidate in replay-level validation.** The N=126 candidate is preliminary, not selected for production yet.
5. **Require per-setting, per-foil, delta-slice, and acceptance-edge checks.** The final choice should be based on reconstructed foil positions and widths, angular residuals, tail behavior, and stability under leaving out a setting—not only on the calibration-sample RMS.
6. **Choose any cutoff or regularization deliberately.** It should be documented as an analysis choice, not inherited from ROOT's floating-point default.

The defensible collaboration statement today is:

> The present angular matrix reconstructs typical events reasonably, but its high-order coefficients are not numerically unique. A properly rescaled direct solve is much better conditioned, yet retaining every sixth-order direction overfits sparse regions. We should validate a controlled lower-order or explicitly regularized rescaled solution—beginning with the complete fifth-order basis—before assigning a production matrix.

## Supporting plots and tables

- `solver_comparison/condition_by_solve_method_vs_N.png`
- `solver_comparison/discarded_directions_by_solve_method_vs_N.png`
- `solver_comparison/candidate_residuals_fit.png`
- `solver_comparison/candidate_residuals_heldout.png`
- `solver_comparison/candidate_vertex_residuals_heldout.png`
- `solver_comparison/angular_solver_comparison.txt`
- `solver_comparison/angular_solver_comparison.tsv`
- `solver_comparison/angular_solver_comparison_heldout.tsv`
- `solver_comparison/angular_candidate_N_comparison.tsv`
