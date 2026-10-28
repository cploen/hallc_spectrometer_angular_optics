# Frozen matrix comparison

**Matrix-file limitation (gmm):** Known uninitialized fixed-xtar row from the historical parser bug, documented in 08_preliminary_conditioning/PLAIN_LANGUAGE_CLARIFICATIONS.md section 10. Exclude this invalid exponent code only; preserve every valid historical fitted coefficient. This is not the corrected-parser refit and does not establish how HCANA treated this corrupt row. Exact excluded row(s) and the original file checksum are recorded in manifest.json.

Order: old replay plus external offsets → GMM fit plus external offsets → full-basis core SVD → elastic-net-selected SVD refit → beam.

## protected_core

| Matrix | xptar RMS (mrad) | ytar RMS (cm) | yptar RMS (mrad) | ztar RMS (cm) |
|---|---:|---:|---:|---:|
| Old replay + offsets | 2.5092 | 0.14229 | 1.242 | 0.6669 |
| GMM SVD | 1.6338 | 0.13746 | 1.0358 | 0.64484 |
| Full core SVD | 1.4528 | 0.13592 | 1.0304 | 0.63777 |
| EN-selected SVD | 1.5124 | 0.13739 | 1.039 | 0.64507 |
| Beam | 1.459 | 0.13631 | 1.033 | 0.63967 |

## protected_noncore

| Matrix | xptar RMS (mrad) | ytar RMS (cm) | yptar RMS (mrad) | ztar RMS (cm) |
|---|---:|---:|---:|---:|
| Old replay + offsets | 4.4704 | 0.17252 | 1.7799 | 0.80751 |
| GMM SVD | 3.9833 | 0.16527 | 1.5817 | 0.77512 |
| Full core SVD | 3.8724 | 0.16503 | 1.5887 | 0.77474 |
| EN-selected SVD | 3.9108 | 0.16614 | 1.5906 | 0.78021 |
| Beam | 3.8766 | 0.16551 | 1.5891 | 0.7774 |

## surplus_core

| Matrix | xptar RMS (mrad) | ytar RMS (cm) | yptar RMS (mrad) | ztar RMS (cm) |
|---|---:|---:|---:|---:|
| Old replay + offsets | 2.3641 | 0.13611 | 1.2195 | 0.6487 |
| GMM SVD | 1.5359 | 0.13397 | 1.0239 | 0.63869 |
| Full core SVD | 1.4366 | 0.1331 | 1.0212 | 0.63474 |
| EN-selected SVD | 1.489 | 0.13431 | 1.0271 | 0.64085 |
| Beam | 1.4418 | 0.13333 | 1.022 | 0.63592 |

## What these comparisons establish

Every matrix is evaluated on the same saved focal-plane coordinates and xtar input, with its entire angular polynomial and explicit external corrections. No coefficient or offset is fitted on these evaluation events. Matrices and offset sources are frozen in this output. This is a common-input reconstruction comparison, not a full iterative HCANA replay with new event selection.

Old → GMM quantifies the historical change. GMM → core combines event selection, balancing and solver changes; it does not isolate tighter cores alone. Core → EN-selected SVD → beam holds the training events, fixed terms and direct-X solver constant and measures the basis-change effect. EN-selected SVD is the all-training refit of the chosen compact seed basis, not the penalized elastic-net prediction. Its identity is the beam source matrices/start.dat.

External offsets must come from historical settings, not from centering these residuals. tsv/offsets.tsv lists additions in mrad/cm. New core, EN and beam fits already contain freely fitted intercepts, so no historical correction is added to those by default. Bias, spread and RMS are reported separately: RMS² = bias² + spread². A lower RMS caused only by removing a constant offset is not an improvement in mean-subtracted width. Spread is not automatically an experimental fitted Gaussian resolution. Local bias/spread must be examined as well as pooled values, where opposing cell biases can cancel.

Protected cores and protected noncore events are separate. Surplus cores are additional development coverage, not an independent final test. Unsupported labels are tabulated separately; blocked positions have counts only. Shoulder residuals depend on the validity of inherited labels.

tsv/overlap.tsv identifies events used by the historical GMM fit, based on reconstructed membership rather than an original solver ID log. tsv/summary.tsv and changes.tsv also report the common subset outside that reconstructed GMM training membership. Unknown membership is not counted as unseen. Training membership for the old replay matrix remains unknown. Thus the all-events historical comparison is descriptive; new-fit held-out performance must not be confused with an independently held-out test of every historical matrix.

Read plots/README.md for the central residual, tail, extreme-tail and foil/delta figures. Central bins and tail probabilities are percentages of the entire pool. Percent-change heat maps complement absolute RMS maps. Full bias, spread, median and P90 down to setting/hole are in tsv/residuals.tsv; qa_low flags N<10 and never changes the sample.

tsv/changes.tsv reports both previous-stage and starting-matrix differences; negative width changes mean improvement. Mean cell MSE gives equal weight to populated foil/delta cells; other pooled statistics weight events equally. tsv/training_conditioning.tsv preserves the source training-design condition numbers for the modern fits only. These are not condition numbers of the held-out residuals, and historical effective rank cannot be inferred just by counting nonzero polynomial coefficients.

Once these protected results guide additional tuning, subsequent performance on them is development evidence. This report makes no automatic claim of improvement.

Historical membership method: Sequential TFit order; 1000 per foil/delta/Y column, 15000 per rungroup/foil, 200000 global
