# Equal-foil weighting with smooth Huber

Four prespecified fits were completed on the same development and reserved events as the unweighted trial. All cores: 310,513 training events. Cores plus supported shoulders: 396,320. Evaluation: 77,864 reserved cores and 22,532 reserved shoulders. No event was discarded to balance a foil.

## Weighting definition

Each event receives `w_i = N / (5 N_foil)`, so every physical foil has exactly 20% of the total base weight. Minimize `sum(w_i rho(r_i / sigma)) / sum(w_i)`; the weight is outside the loss, so it does not change the residual threshold. Huber uses `c = 1.5` and the previous balanced-training MAD scales: 1.559451, 0.1406577, 1.083016 (mrad, cm, mrad). All 210 polynomial terms, fitted constants, inherited xtar/delta coefficients, and geometry remain fixed in definition. Weights are constructed from training counts only. Delta slices and holes are not separately balanced in this first follow-up. The original balanced sample also applied slice/hole caps; equal-foil weights alone do not recreate that allocation.

## Reserved cores: RMS

| Fit | xptar (mrad) | ytar (cm) | yptar (mrad) |
|---|---:|---:|---:|
| Original balanced SVD | 1.488619 | 0.135820 | 1.029166 |
| Unweighted cores + shoulders, Huber | 1.465615 | 0.134538 | 1.026959 |
| Equal-foil cores, squared | 1.483745 | 0.135024 | 1.027723 |
| Equal-foil cores, Huber | 1.484112 | 0.134997 | 1.027491 |
| Equal-foil cores + shoulders, squared | 1.573990 | 0.135489 | 1.044300 |
| Equal-foil cores + shoulders, Huber | 1.518633 | 0.135417 | 1.037378 |

## Reserved shoulders: RMS

| Fit | xptar (mrad) | ytar (cm) | yptar (mrad) |
|---|---:|---:|---:|
| Original balanced SVD | 4.318663 | 0.167218 | 1.675688 |
| Unweighted cores + shoulders, Huber | 4.214493 | 0.164664 | 1.640240 |
| Equal-foil cores, squared | 4.316905 | 0.166354 | 1.651778 |
| Equal-foil cores, Huber | 4.317120 | 0.166670 | 1.651974 |
| Equal-foil cores + shoulders, squared | 4.134985 | 0.164640 | 1.615737 |
| Equal-foil cores + shoulders, Huber | 4.205898 | 0.164985 | 1.621754 |

## Changes relative to the original balanced SVD

Positive means smaller error. Equal-cell RMS is the square root of the mean of 25 foil/delta-cell MSEs. Within-hole width is the square root of the event-weighted mean within-hole variance, after subtracting each run/foil/delta/hole mean; only common groups with at least 10 reserved events are used. It is a descriptive width, not fitted Gaussian resolution.

### core

| Fit | Target | Pooled RMS gain (%) | Equal-cell RMS gain (%) | Improved cells | Within-hole width gain (%) |
|---|---|---:|---:|---:|---:|
| Equal-foil cores, squared | xptar | 0.327 | 0.336 | 18/25 | 0.162 |
| Equal-foil cores, squared | ytar | 0.586 | 0.473 | 17/25 | 0.593 |
| Equal-foil cores, squared | yptar | 0.140 | 0.345 | 17/25 | 0.735 |
| Equal-foil cores, Huber | xptar | 0.303 | 0.277 | 18/25 | 0.222 |
| Equal-foil cores, Huber | ytar | 0.606 | 0.428 | 18/25 | 0.676 |
| Equal-foil cores, Huber | yptar | 0.163 | 0.345 | 18/25 | 0.770 |
| Equal-foil cores + shoulders, squared | xptar | -5.735 | -4.661 | 1/25 | 1.244 |
| Equal-foil cores + shoulders, squared | ytar | 0.243 | 0.136 | 17/25 | 0.784 |
| Equal-foil cores + shoulders, squared | yptar | -1.470 | -0.444 | 14/25 | 1.412 |
| Equal-foil cores + shoulders, Huber | xptar | -2.016 | -1.192 | 8/25 | 1.078 |
| Equal-foil cores + shoulders, Huber | ytar | 0.296 | 0.186 | 15/25 | 0.965 |
| Equal-foil cores + shoulders, Huber | yptar | -0.798 | -0.122 | 15/25 | 1.320 |

### shoulder

| Fit | Target | Pooled RMS gain (%) | Equal-cell RMS gain (%) | Improved cells | Within-hole width gain (%) |
|---|---|---:|---:|---:|---:|
| Equal-foil cores, squared | xptar | 0.041 | 0.081 | 15/25 | 0.184 |
| Equal-foil cores, squared | ytar | 0.517 | 0.634 | 18/25 | 0.487 |
| Equal-foil cores, squared | yptar | 1.427 | 1.565 | 21/25 | 2.311 |
| Equal-foil cores, Huber | xptar | 0.036 | 0.120 | 16/25 | 0.209 |
| Equal-foil cores, Huber | ytar | 0.328 | 0.365 | 16/25 | 0.637 |
| Equal-foil cores, Huber | yptar | 1.415 | 1.523 | 21/25 | 2.398 |
| Equal-foil cores + shoulders, squared | xptar | 4.253 | 4.968 | 20/25 | 3.487 |
| Equal-foil cores + shoulders, squared | ytar | 1.542 | 1.553 | 22/25 | 1.937 |
| Equal-foil cores + shoulders, squared | yptar | 3.578 | 3.725 | 21/25 | 5.319 |
| Equal-foil cores + shoulders, Huber | xptar | 2.611 | 3.068 | 23/25 | 2.104 |
| Equal-foil cores + shoulders, Huber | ytar | 1.335 | 1.322 | 22/25 | 2.091 |
| Equal-foil cores + shoulders, Huber | yptar | 3.219 | 3.201 | 22/25 | 4.757 |

## Residual weights after foil balancing

For smooth Huber the IRLS factor is `h_i = 1/sqrt(1+(r_i/(c sigma))²)`. The following shares sum `w_i h_i`; they describe weights in the estimating equation, not full coefficient influence (which also depends on the polynomial design). No second foil normalization is applied after Huber weighting.

### all_core

| Foil (cm) | Base share | xptar IRLS share | ytar IRLS share | yptar IRLS share |
|---|---:|---:|---:|---:|
| -8 | 20.00% | 19.97% | 19.76% | 19.99% |
| -3 | 20.00% | 19.81% | 20.06% | 20.03% |
| +0 | 20.00% | 20.04% | 20.23% | 20.17% |
| +3 | 20.00% | 20.01% | 20.14% | 20.01% |
| +8 | 20.00% | 20.17% | 19.82% | 19.79% |

### all_supported

| Foil (cm) | Base share | xptar IRLS share | ytar IRLS share | yptar IRLS share |
|---|---:|---:|---:|---:|
| -8 | 20.00% | 19.71% | 19.67% | 20.01% |
| -3 | 20.00% | 19.84% | 20.11% | 20.03% |
| +0 | 20.00% | 20.05% | 20.25% | 20.17% |
| +3 | 20.00% | 20.13% | 20.16% | 19.99% |
| +8 | 20.00% | 20.27% | 19.81% | 19.80% |

### Shoulder contributions in the full supported sample

| Target | Raw event share | Base-weight share | IRLS-weight share | Mean residual factor |
|---|---:|---:|---:|---:|
| xptar | 21.65% | 21.56% | 16.90% | 0.648 |
| ytar | 21.65% | 21.56% | 20.81% | 0.837 |
| yptar | 21.65% | 21.56% | 19.64% | 0.777 |

## Files and limits

[Comparison](plots/comparison.png) · [Foil/delta changes](plots/foil_delta.png) · [Effect of adding weights](plots/weighting_effect.png) · [Foil weights](plots/foil_weights.png).

`tsv/summary.tsv` contains pooled, equal-cell and within-hole statistics for all nine reference/candidate fits. `tsv/metrics.tsv` includes biases, spreads and tails; `tsv/holes.tsv` preserves sparse groups; `tsv/weight_contributions.tsv` contains totals, foil shares, core/shoulder breakdowns and Kish effective counts. `training_weights.npz` stores exact IDs and base weights. `tsv/paired_intervals.tsv` provides conditional paired cluster-bootstrap intervals for pooled RMS changes against the original balanced fit; it does not include training or systematic uncertainty.

All four candidates were specified before viewing their reserved scores. Source inputs and residual identities were verified, all optimizations converged, and exported matrices reproduce predictions. This follows earlier examination of the same reserved sample and is not fresh independent validation. Shoulder labels, target geometry and saved xtar remain inherited. No replay matrix was installed; these are fixed-input reconstruction comparisons.
