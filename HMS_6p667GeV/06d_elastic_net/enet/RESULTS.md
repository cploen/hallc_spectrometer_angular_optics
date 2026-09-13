# Elastic-net angular fit: HMS_6p667GeV

68,535 balanced training cores; 3 training-only validation folds.

| Target | Alpha | L1 ratio | Active slopes | CV MSE / scaled SVD | Alpha at grid edge |
|---|---:|---:|---:|---:|---|
| xptar | 0.0003 | 0.5 | 59 | 1.699 | False |
| ytar | 0.001 | 0.5 | 62 | 1.860 | False |
| yptar | 0.0003 | 0.5 | 55 | 1.108 | False |

Scores average squared error equally over populated physical-foil/delta cells within each fold. The one-standard-error rule favors fewer terms; its error bar is a fold-variability heuristic, not a confidence interval.

Compare [validation curves](plots/cv.png), [term stability](plots/terms.png), [training strata](tsv/training.tsv), and [term candidates](tsv/terms.tsv). The scaled SVD baseline uses identical events, basis, fixed terms and training-only transformations; it is not the historical ROOT normal-equation fit.

Beam search is a possible development step, not a demonstrated improvement. Use fold frequency and standardized coefficients as candidate rankings, not proof of physical importance. Keep the constant and fixed xtar/delta terms. Use these saved training folds, recomputing scaling in each fold. Do not rank terms or tune search settings on protected core/noncore results. Prefer completing any beam-search study before the one-time protected evaluation. If protected diagnostics subsequently guide changes, they become development evidence and are no longer an untouched final test.

Do not advance a matrix based on pooled RMS alone: inspect edge-hole counts/residuals and all foil/delta cells. Grid-edge selections, unstable term membership, or little CV advantage over SVD call for further training-only study. Candidate matrices have not been installed in replay. Reconstructed quantities here keep the saved replay xtar input fixed; full replay/iteration validation remains necessary.
