# Beam-search development results

| Target | Seed terms | Beam terms | Seed MSE | Beam MSE | Full-basis MSE | Stop |
|---|---:|---:|---:|---:|---:|---|
| xptar | 95 | 110 | 3.08344 | 2.80502 | 2.77347 | step limit |
| ytar | 132 | 142 | 0.023472 | 0.0230283 | 0.0228738 | step limit |
| yptar | 75 | 78 | 1.20312 | 1.19302 | 1.18081 | no material gain |

MSE is averaged equally over populated physical foil/delta cells. Its units are mrad² for angles and cm² for ytar. Median and P90 describe event-pooled absolute residuals. A hole with six events does not receive the same weight as an entire well-populated cell, and sparse-hole flags never veto a candidate. No events are omitted from scoring.

plots/search.png shows best score reached, terms, and the worst fold condition number; stars mark the final choice, which can differ from the lowest-MSE candidate because of the configured sparsity allowance. The full-basis reference always uses these same events and fold scaling, not the historical ROOT normal-equation solve on a different sample. The constant is counted as a term but excluded from the centered slope-block condition number.

plots/residuals.png shows signed residual distributions. plots/foil_delta.png shows absolute RMS in mrad or cm, with a shared color scale for the three models within each target. RMS includes bias; spread in tsv/residuals.tsv is standard deviation after subtracting the mean, not a fitted Gaussian detector resolution. ztar uses the saved-xtar reconstruction diagnostic and is not a new replay.

tsv/residuals.tsv gives bias, spread, RMS, median, P90 and N down to individual setting/hole. qa_low marks fewer than ten events and is informational. tsv/summary.tsv and fold_scores.tsv give overall and individual-fold scores; conditioning.tsv gives ranks and conditioning by fold, final_fit.tsv the final full-training solves. residuals.npz preserves exact event IDs and physical-unit residuals for all three models, ordered xptar, ytar, yptar, ztar.

**These are adaptive development scores.** Candidate term choices use the same three folds that score the search, and the original elastic-net seeds also depended on these training events. Each prediction uses coefficients fitted without that event, but the basis selection is not independent of the scored events. There are no confidence intervals or claims of independent generalization. Protected core and noncore pools remain closed.

Each target starts from all distinct saved bases for its chosen elastic-net setting. The Seed reference is the best of those fixed bases under this shared-fold scoring, not the earlier pooled result that used a different seed in each fold. No averaging or intersection of the seed masks is used.

matrices/beam.dat is the chosen basis refitted on all balanced training events. matrices/start.dat and svd.dat provide corresponding seed and full-basis refits. seeds/fold*.dat preserve the actual input seed matrices; seed.dat preserves the fixed transport terms. Delta coefficients and xtar-dependent terms are unchanged. No matrix has been installed in replay or evaluated on protected events.
