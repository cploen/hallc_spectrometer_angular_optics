# Elastic-net convergence study

Saved fold 0 of enet; tolerance 1e-06, unchanged. Each case starts at zero and continues the same objective between checkpoints.

| Target | L1 | Alpha | Iterations | Converged | Gap / limit | Validation MSE / SVD |
|---|---:|---:|---:|---|---:|---:|
| xptar | 0.1 | 0.0001 | 74423 | True | 1 | 1.292 |
| xptar | 0.1 | 3e-05 | 354669 | True | 1 | 1.145 |
| xptar | 0.5 | 0.0001 | 34547 | True | 1 | 1.413 |
| xptar | 0.5 | 3e-05 | 96960 | True | 1 | 1.214 |
| ytar | 0.1 | 0.0001 | 489158 | True | 1 | 1.139 |
| ytar | 0.1 | 3e-05 | 1081997 | True | 1 | 1.054 |
| ytar | 0.5 | 0.0001 | 120233 | True | 1 | 1.14 |
| ytar | 0.5 | 3e-05 | 421450 | True | 1 | 1.056 |
| yptar | 0.1 | 0.0001 | 155579 | True | 1 | 1.027 |
| yptar | 0.1 | 3e-05 | 508963 | True | 1 | 1.011 |
| yptar | 0.5 | 0.0001 | 32968 | True | 1 | 1.046 |
| yptar | 0.5 | 3e-05 | 75189 | True | 1 | 1.017 |

See convergence.png and checkpoints.tsv. Crosses remain unfinished even if their validation scores look good. One-fold scores guide the next experiment; they do not select a final matrix. No protected events were fitted or evaluated, and no matrix was exported. The SVD reference solves scaled X directly.

If convergence is reached and scores stabilize, extend the useful settings to all saved folds. If not, inspect objective, gap, KKT residual, prediction changes and runtime before raising the budget again or changing solver. Do not loosen tolerance to obtain a passing status. Starts differ from the original descending-alpha warm path, so this is a controlled continuation study rather than an exact replay of its 20,000-iteration fit.

## Unpenalized selected-term SVD comparison

| Target | Alpha | L1 | Terms including constant | Penalized MSE / full SVD | Refit MSE / full SVD |
|---|---:|---:|---:|---:|---:|
| xptar | 0.0001 | 0.1 | 132 | 1.292 | 1.04 |
| xptar | 3e-05 | 0.1 | 151 | 1.145 | 1.011 |
| xptar | 0.0001 | 0.5 | 80 | 1.413 | 1.325 |
| xptar | 3e-05 | 0.5 | 95 | 1.214 | 1.116 |
| ytar | 0.0001 | 0.1 | 166 | 1.139 | 1.009 |
| ytar | 3e-05 | 0.1 | 178 | 1.054 | 1.004 |
| ytar | 0.0001 | 0.5 | 112 | 1.14 | 1.033 |
| ytar | 3e-05 | 0.5 | 134 | 1.056 | 1.03 |
| yptar | 0.0001 | 0.1 | 129 | 1.027 | 1.004 |
| yptar | 3e-05 | 0.1 | 145 | 1.011 | 1.001 |
| yptar | 0.0001 | 0.5 | 75 | 1.046 | 1.019 |
| yptar | 3e-05 | 0.5 | 98 | 1.017 | 1.011 |

Below 1 means lower validation MSE than full-basis SVD. Only converged elastic-net cases are refitted. Term selection, centering, scaling and both SVD solves use only the fitting portion of the fold. The constant stays free of penalties; fixed xtar terms and delta coefficients are retained in seed.dat. A score near or below 1 is a candidate for the three-fold study, not a final matrix choice.

comparison.tsv records terms, rank and scores; validation.tsv records N, bias, RMS and P95 by physical foil/delta, setting and setting/hole. ztar is included only where all three target fits converged for the same penalty settings. terms.tsv and coefficients.npz preserve the fitting-fold selection and native coefficients. These are development fits, not full-training replay matrices.

See plots/comparison.png and each case foil/delta and hole map. Error ratios use common scales; black outlines flag sparse validation holes, and blanks mean no validation events. Counts are not training allocation counts. Sparse populations may be absent from this fold. Complete all folds before choosing the basis; keep the protected pools closed.
