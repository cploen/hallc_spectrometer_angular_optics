# Elastic-net convergence study

Saved fold 0 of enet; tolerance 1e-06, unchanged. Each case starts at zero and continues the same objective between checkpoints.

| Target | L1 | Alpha | Iterations | Converged | Gap / limit | Validation MSE / SVD |
|---|---:|---:|---:|---|---:|---:|
| xptar | 0.1 | 0.0001 | 74423 | True | 1 | 1.292 |
| xptar | 0.1 | 3e-05 | 300000 | False | 1.99 | 1.145 |
| xptar | 0.5 | 0.0001 | 34547 | True | 1 | 1.413 |
| xptar | 0.5 | 3e-05 | 96960 | True | 1 | 1.214 |
| ytar | 0.1 | 0.0001 | 300000 | False | 26.7 | 1.139 |
| ytar | 0.1 | 3e-05 | 300000 | False | 761 | 1.054 |
| ytar | 0.5 | 0.0001 | 120233 | True | 1 | 1.14 |
| ytar | 0.5 | 3e-05 | 300000 | False | 6.22 | 1.056 |
| yptar | 0.1 | 0.0001 | 155579 | True | 1 | 1.027 |
| yptar | 0.1 | 3e-05 | 300000 | False | 22.2 | 1.011 |
| yptar | 0.5 | 0.0001 | 32968 | True | 1 | 1.046 |
| yptar | 0.5 | 3e-05 | 75189 | True | 1 | 1.017 |

See convergence.png and checkpoints.tsv. Crosses remain unfinished even if their validation scores look good. One-fold scores guide the next experiment; they do not select a final matrix. No protected events were fitted or evaluated, and no matrix was exported. The SVD reference solves scaled X directly.

If convergence is reached and scores stabilize, extend the useful settings to all saved folds. If not, inspect objective, gap, KKT residual, prediction changes and runtime before raising the budget again or changing solver. Do not loosen tolerance to obtain a passing status. Starts differ from the original descending-alpha warm path, so this is a controlled continuation study rather than an exact replay of its 20,000-iteration fit.
