# Frozen matrix comparison

**Matrix-file limitation (gmm):** Known uninitialized fixed-xtar row from the historical parser bug, documented in 08_preliminary_conditioning/PLAIN_LANGUAGE_CLARIFICATIONS.md section 10. Exclude this invalid exponent code only; preserve every valid historical fitted coefficient. This is not the corrected-parser refit and does not establish how HCANA treated this corrupt row. Exact excluded row(s) and the original file checksum are recorded in manifest.json.

Order: old replay plus external offsets → GMM fit plus external offsets → full-basis core SVD → elastic-net-selected SVD refit → beam.

## protected_core

| Matrix | xptar RMS (mrad) | ytar RMS (cm) | yptar RMS (mrad) | ztar RMS (cm) |
|---|---:|---:|---:|---:|
| Old replay + offsets | 1.5709 | 0.14281 | 1.0487 | 0.66875 |
| GMM fit + offsets | 1.5367 | 0.13726 | 1.0344 | 0.64396 |
| Full core SVD | 1.4886 | 0.13582 | 1.0292 | 0.63741 |
| EN-selected SVD | 1.5538 | 0.13713 | 1.0363 | 0.64389 |
| Beam | 1.4956 | 0.13602 | 1.0329 | 0.63853 |

## protected_noncore

| Matrix | xptar RMS (mrad) | ytar RMS (cm) | yptar RMS (mrad) | ztar RMS (cm) |
|---|---:|---:|---:|---:|
| Old replay + offsets | 4.3981 | 0.17488 | 1.6882 | 0.81838 |
| GMM fit + offsets | 4.3157 | 0.16644 | 1.6595 | 0.78139 |
| Full core SVD | 4.3187 | 0.16722 | 1.6757 | 0.78729 |
| EN-selected SVD | 4.3346 | 0.1674 | 1.6756 | 0.7865 |
| Beam | 4.3222 | 0.16692 | 1.6736 | 0.78439 |

## surplus_core

| Matrix | xptar RMS (mrad) | ytar RMS (cm) | yptar RMS (mrad) | ztar RMS (cm) |
|---|---:|---:|---:|---:|
| Old replay + offsets | 1.4816 | 0.13771 | 1.0299 | 0.65581 |
| GMM fit + offsets | 1.498 | 0.13388 | 1.0234 | 0.63822 |
| Full core SVD | 1.4716 | 0.13319 | 1.0213 | 0.6352 |
| EN-selected SVD | 1.5293 | 0.13419 | 1.0257 | 0.64022 |
| Beam | 1.4779 | 0.13323 | 1.0227 | 0.63545 |

## What these comparisons establish

Every matrix is evaluated on the same saved focal-plane coordinates and xtar input, with its entire angular polynomial and explicit external corrections. No coefficient or offset is fitted on these evaluation events. Matrices and offset sources are frozen in this output. This is a common-input reconstruction comparison, not a full iterative HCANA replay with new event selection.

Old → GMM quantifies the historical change. GMM → core combines event selection, balancing and solver changes; it does not isolate tighter cores alone. Core → EN-selected SVD → beam holds the training events, fixed terms and direct-X solver constant and measures the basis-change effect. EN-selected SVD is the all-training refit of the chosen compact seed basis, not the penalized elastic-net prediction. Its identity is the beam source matrices/start.dat.

External offsets must come from historical settings, not from centering these residuals. tsv/offsets.tsv lists additions in mrad/cm. New core, EN and beam fits already contain freely fitted intercepts, so no historical correction is added to those by default. Bias, spread and RMS are reported separately: RMS² = bias² + spread². A lower RMS caused only by removing a constant offset is not an improvement in mean-subtracted width. Spread is not automatically an experimental fitted Gaussian resolution. Local bias/spread must be examined as well as pooled values, where opposing cell biases can cancel.

Protected cores and protected noncore events are separate. Surplus cores are additional development coverage, not an independent final test. Unsupported labels are tabulated separately; blocked positions have counts only. Shoulder residuals depend on the validity of inherited labels.

tsv/overlap.tsv identifies events used by the historical GMM fit, based on reconstructed membership rather than an original solver ID log. tsv/summary.tsv and changes.tsv also report the common subset outside that reconstructed GMM training membership. Unknown membership is not counted as unseen. Training membership for the old replay matrix remains unknown. Thus the all-events historical comparison is descriptive; new-fit held-out performance must not be confused with an independently held-out test of every historical matrix.

Read [plots/README.md](plots/README.md) for the revised figures. Plotting changes do not change the saved residuals or numerical tables.

tsv/changes.tsv reports both previous-stage and starting-matrix differences; negative width changes mean improvement. Mean cell MSE gives equal weight to populated foil/delta cells; other pooled statistics weight events equally. tsv/training_conditioning.tsv preserves the source training-design condition numbers for the modern fits only. These are not condition numbers of the held-out residuals, and historical effective rank cannot be inferred just by counting nonzero polynomial coefficients.

Once these protected results guide additional tuning, subsequent performance on them is development evidence. This report makes no automatic claim of improvement.

Historical membership method: Sequential TFit order; 1000 per foil/delta/Y column, 15000 per rungroup/foil, 200000 global
