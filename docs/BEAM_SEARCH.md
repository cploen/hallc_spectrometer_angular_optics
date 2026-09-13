# Beam search from elastic-net seed bases

This step searches combinations of angular polynomial terms. **Every candidate
is refitted by unpenalized SVD on scaled X.** Elastic-net penalties are used only
to obtain the starting term sets; they do not enter the beam-search fits.

## Run on ifarm

After pulling the code, run from the repository root:

```bash
./run_beam.sh HMS_6p667GeV equal15 beam --check
./run_beam.sh HMS_6p667GeV equal15 beam --beam 4 --threads 8
```

The default source is `06d_elastic_net/pooled/`. It must be a completed pooled
study with the same balanced training membership. The first command verifies
source products, event IDs, and seed selections without fitting or writing
outputs. It needs the **fit** TFit exports, not the protected holdout exports.
The second writes `HMS_6p667GeV/06e_beam_search/beam/`. Outputs are immutable:
use a new output name for another study, for example:

```bash
./run_beam.sh HMS_6p667GeV equal15 wide --beam 8 --threads 8 --steps 30
```

`--source NAME` chooses another pooled elastic-net study. `--beam` changes the
number of paths retained each round; `--threads` changes candidate concurrency
and defaults to eight. `--steps` limits rounds per target. Other small policy
choices live in `config/beam.json`; command-line values override that file.

For another HMS momentum campaign, substitute its directory, balanced sample
tag and pooled source. Seed penalties can be changed in that campaign's config.
The search engine has no hard-coded foil locations, hole counts or delta
boundaries. The existing event/geometry adapter currently supports HMS only;
SHMS needs that adapter extended before use, rather than silently applying HMS
geometry.

## Seed handoff

For HMS 6.667 the agreed starting settings are:

| Target | Color in the pooled study | Alpha | L1 fraction |
|---|---|---:|---:|
| xptar | Red | 0.00003 | 0.5 |
| ytar | Red | 0.00003 | 0.5 |
| yptar | Green | 0.0001 | 0.5 |

The pooled study already provides `terms.tsv`, `coefficients.npz`, `seed.dat`
and `folds.tsv`. Beam search reads them directly, with checksum verification.
No manual copying or conversion is needed. It starts from **all distinct bases
selected in the three folds** for each chosen setting, rather than intersecting,
averaging or choosing one fold's selection. Search is independent per target.

`seeds/fold0.dat`, `fold1.dat`, etc. preserve actual native-unit seed matrices
using those selected coefficients; `seeds/terms.tsv` preserves their masks.
These seed coefficients were fitted on the corresponding two-thirds sample.
They identify the handoff, but every candidate is freshly refitted in every
scoring fold. The reported Seed reference is the best of these initial fixed
bases under the new shared-fold scoring. It can differ slightly from the old
pooled elastic-net score, which used a different basis in each fold.

## Search and stopping rules

1. Keep the constant and all fixed xtar-dependent transport terms. The search
   can add or remove any of the original fitted nonconstant columns. Previously
   excluded elastic-net terms may return. Delta is not fitted. No new polynomial
   powers beyond the original candidate dictionary are invented in this step.
2. For each basis, fit coefficients on two folds and predict the third. Repeat
   for every saved fold. Assemble one coefficient-excluded prediction per event.
3. Score the pooled residuals by averaging MSE equally across populated physical
   foil/delta cells. Within each cell, events contribute equally. This preserves
   broad acceptance coverage without giving a six-event hole the weight of an
   entire well-populated cell. No extra event balancing or resampling occurs.
4. Evaluate all one-term additions and removals from the current frontier. Keep
   the best `beam` new candidates and retain a separate record of the best basis
   ever seen. A removal followed by an addition allows replacement; the frontier
   can temporarily worsen, so replacements need not improve at their first move.
   Identical bases are evaluated once per target, with deterministic ordering.
5. Stop after `steps` rounds (default 20), exhaustion of unvisited neighbours,
   or `patience` rounds (default 3) without a cumulative relative improvement of
   `gain` (default 0.001, or 0.1%) over the last improvement anchor.
6. Among all evaluated, full-rank candidates within `slack` (default 0.005,
   or 0.5%) of the best MSE, choose the smallest. Break equal-size ties by MSE,
   then condition number and term ordering. Set slack to zero to select strictly
   by MSE. Conditioning is a diagnostic and tie-breaker, not a hidden penalty or
   cutoff. A candidate must retain full rank at the source's saved SVD cutoff.

The choices above are explicit, adjustable heuristics, not statistical confidence
limits. A finite beam can miss useful combinations, especially combinations that
require several temporarily unhelpful moves. Increase beam width or patience
only if the first results give a reason. A step-limit stop means the budget was
used, not that a plateau has been proven. There is no requirement to beat the
full-basis reference before a candidate can be retained.

## SVD and performance

Centering and column scaling are learned separately on each fold's fitting
events. To avoid repeatedly decomposing tens of thousands of rows, an orthogonal
QR factorization compresses `[X, y]` once per fold. Each candidate is solved by
SVD on the selected columns of that R factor. Orthogonal compression preserves
the least-squares problem and singular values of X; **X-transpose-X is not used
for the solve**. Validation residual norms are similarly compressed, with
weights that reproduce the exact pooled mean foil/delta MSE. Actual event-level
residuals are reconstructed for the seed, chosen and full-basis diagnostics.

Candidates run concurrently on eight worker threads by default, with one BLAS
thread per worker to avoid multiplying the thread count. One-time preparation
and final full-training solves can use the configured thread count. Per-round
elapsed time and evaluation counts are recorded in `tsv/search.tsv`. Eight
threads are not a promise of eightfold speedup; farm runtime depends on hardware
and how many terms and rounds are evaluated.

Final coefficients are fitted on **all** balanced training cores using the
selected target-specific masks and direct scaled-X SVD. `matrices/beam.dat`,
`start.dat` and `svd.dat` are candidate transport matrices. Unselected angular
coefficients are zero; rows can remain for other targets or retained delta
coefficients. Matrix round-trip predictions are checked before completion.
Nothing is installed in replay, and no protected pool is read.

## Read the results

- `RESULTS.md`: concise scores, selected term counts and stop reasons.
- `plots/search.png`: best MSE reached, term count, and worst fold condition
  number through the search. Stars identify the chosen smaller basis, which can
  differ from the absolute best-MSE basis because of slack.
- `plots/residuals.png`: signed residual distributions for seed, beam and full
  basis; includes xptar, ytar, yptar and the derived ztar diagnostic.
- `plots/foil_delta.png`: absolute RMS for every foil/delta cell and each model,
  with units and event counts. Color scales are shared between models for each
  target. This makes a change such as 1.9 to 2.1 mrad visible directly.
- `tsv/summary.tsv`, `fold_scores.tsv`: mean cell MSE, median absolute residual,
  P90, RMS, bias and spread, pooled and per fold. MSE uses squared physical units.
- `tsv/residuals.tsv`, `coverage.tsv`: corresponding local diagnostics and
  counts down to setting/hole. `qa_low` flags fewer than ten events for reporting
  only. These sparse estimates do not veto widespread improvement.
- `tsv/conditioning.tsv`, `final_fit.tsv`: term counts, ranks and condition
  numbers for fold fits and final full-training fits. The separately handled
  constant is counted in terms but excluded from slope-block conditioning.
- `tsv/terms.tsv`, `candidates.tsv`, `search.tsv`, `search.json`: selected masks,
  evaluated candidates and search history. Five-digit term codes are strings.
- `model.npz`, `residuals.npz`, `matrices/`, `seeds/`, `seed.dat`, `tsv/folds.tsv`,
  `manifest.json`, `code/`: coefficients, exact residual/event IDs, source seed
  snapshots, effective settings and checksums for reproduction.

RMS includes systematic offsets. `spread` is the standard deviation after
subtracting the mean; neither is automatically the experimental Gaussian
resolution requirement. The fixed-input ztar diagnostic is not a full replay.

The full-basis reference is the rescaled direct-X solve on these same balanced
cores and saved folds. It isolates the effect of changing the basis. It is not
a claim about improvement over the historical ROOT solve with different events.

**Search scores are adaptive development scores.** The basis is chosen using
these folds, although each event's predicted coefficients were fitted without
that event. Seeds also originated from these training events. The plotted gains
therefore are not independent estimates of generalization. Once the search
configuration and final basis are settled, compare frozen candidates on the
protected cores and noncore events before drawing the final conclusion. Keep
those tests separate from search and retain the original full-basis reference.

To expose a completed run for review, stage the entire
`HMS_6p667GeV/06e_beam_search/beam/` directory, including its manifest and tables;
no additional ROOT files are needed for interpreting the outputs.
