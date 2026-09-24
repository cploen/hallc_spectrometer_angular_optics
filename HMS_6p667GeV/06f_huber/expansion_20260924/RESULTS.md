# Smooth Huber and training-sample expansion

Real HMS 6.667 GeV frozen event data. Candidate matrices are exported but not installed in replay.

## Experimental definition

The regression loss is the smooth (pseudo-)Huber loss, `rho(u) = c² (sqrt(1+(u/c)²)-1)`, with `u = residual / sigma` and `c = 1.5`. `sigma` is 1.4826 times the residual MAD from the balanced-training squared-error fit, frozen for every later fit. Constants are fitted; xtar-dependent and delta coefficients are inherited from the archived seed. No foil, slice or hole population weights are applied. The regression loss follows [SciPy pseudo_huber](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pseudo_huber.html); the classification loss called modified_huber is not used.

The main experiment retains all 210 terms (including the constant). Half increments use a reproducible event-ID hash and retain the entire preceding sample. Additional cores mostly increase statistics within existing accepted regions; shoulders also broaden the occupied sample. The expansion stays inside the original accepted labels and delta range; it cannot recover events absent from the frozen trees.

| Training sample | Events |
|---|---:|
| balanced | 68,535 |
| half_surplus | 189,464 |
| all_core | 310,513 |
| half_shoulders | 353,341 |
| all_supported | 396,320 |

Residual scales (xptar mrad, ytar cm, yptar mrad): 1.55945, 0.140658, 1.08302.

## Common reserved-sample RMS

### Reserved cores

| Training | Loss | xptar (mrad) | ytar (cm) | yptar (mrad) |
|---|---|---:|---:|---:|
| balanced | squared | 1.488619 | 0.135820 | 1.029166 |
| balanced | huber | 1.489359 | 0.135679 | 1.028882 |
| half_surplus | squared | 1.449631 | 0.134368 | 1.023541 |
| half_surplus | huber | 1.449688 | 0.134355 | 1.023376 |
| all_core | squared | 1.448471 | 0.134301 | 1.023150 |
| all_core | huber | 1.448502 | 0.134349 | 1.023146 |
| half_shoulders | squared | 1.474560 | 0.134424 | 1.027414 |
| half_shoulders | huber | 1.453178 | 0.134400 | 1.024795 |
| all_supported | squared | 1.505663 | 0.134579 | 1.031097 |
| all_supported | huber | 1.465615 | 0.134538 | 1.026959 |

### Reserved shoulders

| Training | Loss | xptar (mrad) | ytar (cm) | yptar (mrad) |
|---|---|---:|---:|---:|
| balanced | squared | 4.318663 | 0.167218 | 1.675688 |
| balanced | huber | 4.318268 | 0.167305 | 1.675119 |
| half_surplus | squared | 4.332004 | 0.165888 | 1.675000 |
| half_surplus | huber | 4.331379 | 0.166340 | 1.676470 |
| all_core | squared | 4.342222 | 0.165936 | 1.673868 |
| all_core | huber | 4.341001 | 0.166519 | 1.675342 |
| half_shoulders | squared | 4.194653 | 0.164775 | 1.637049 |
| half_shoulders | huber | 4.272289 | 0.165247 | 1.656051 |
| all_supported | squared | 4.136230 | 0.164282 | 1.623674 |
| all_supported | huber | 4.214493 | 0.164664 | 1.640240 |

## Expanded Huber versus balanced squared error

Positive gain means lower RMS. The intervals below are paired 500-replicate bootstrap intervals over run/foil/delta/hole clusters, stratified by physical foil/delta. They describe evaluation-sample variation conditional on these fitted models, not training uncertainty or systematic calibration error.

| Pool | Target | RMS reduction (%) | Conditional 95% interval | Cell-macro RMS reduction (%) | Improved cells | Worst cell change (%) |
|---|---|---:|---|---:|---:|---:|
| core | xptar | 1.545 | [1.049, 2.065] | -2.043 | 7/25 | -8.825 |
| core | ytar | 0.944 | [0.657, 1.216] | -0.348 | 12/25 | -8.648 |
| core | yptar | 0.214 | [-0.111, 0.495] | -0.329 | 12/25 | -6.124 |
| shoulder | xptar | 2.412 | [1.777, 3.080] | 2.158 | 19/25 | -1.953 |
| shoulder | ytar | 1.528 | [1.143, 1.976] | 0.545 | 14/25 | -6.289 |
| shoulder | yptar | 2.115 | [1.073, 3.301] | 1.769 | 17/25 | -6.365 |

## Elastic net and beam-search controls

The existing elastic-net-selected and beam10 term supports are reused and refitted on all cores and on cores plus shoulders with both losses. The elastic-net-selected basis here has unpenalized coefficients; it is not a new penalized elastic-net fit. No new beam search or hyperparameter optimization is claimed. See `plots/bases.png` and the complete `tsv/metrics.tsv`.

## Checks and limitations

Frozen input hashes and balanced event IDs were verified against the archived beam run. The locally rebuilt balanced SVD predictions reproduce the saved SVD matrix; maximum absolute differences in the reported physical units are [5.121875146230082e-10, 1.5086864049718152e-10, 8.017128627635373e-10]. All Huber optimizations pass an explicit gradient check; candidate matrix files reproduce in-memory predictions.

All variants are fixed before examining the common reserved residuals. The reserved sample was already examined during earlier work, so this is diagnostic evidence, not a new independent final validation. Density labels themselves were estimated from development events. Unsupported and blocked labels are excluded. Shoulder residuals remain conditional on their inherited labels. Original geometry and saved xtar are held fixed; full iterative replay and independent physical calibration are untested. RMS includes bias and is not a Gaussian resolution. Bias, standard deviation, central 68% half-width, 95th absolute residual, and sparse-hole statistics are retained in the TSV files.

Files: `plots/expansion.png`, `plots/foil_delta.png`, `plots/bases.png`, `tsv/coverage.tsv`, `tsv/metrics.tsv`, `tsv/holes.tsv`, `tsv/paired_intervals.tsv`, `tsv/convergence.tsv`, `tsv/membership.tsv`, candidate `matrices/*.dat`, and the input/output manifest.
