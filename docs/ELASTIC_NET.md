# Elastic-net angular optics: balanced HMS samples

This step fits xptar, ytar and yptar from the **same balanced core event IDs**
used for the preallocated SVD handoff. It does not repeat density selection,
change quotas, borrow protected events, or run the historical ROOT event caps.
HMS 6.667 is the first campaign to run. Other HMS campaigns use their own saved
metadata, seed matrix and balance tag; no momentum or five-foil layout is hardcoded.
SHMS needs an adapter/geometry validation before enabling it. The regression
module itself has no spectrometer assumptions.

## Commands on ifarm

Use the repository root and the ROOT/HCANA environment used for core exports.
`PYTHON` selects the Python interpreter, as in the other steps.

```bash
python3 -m pip install --user -r requirements_core.txt
./run_build_core_fit.sh HMS_6p667GeV equal15 fit
./run_build_core_fit.sh HMS_6p667GeV equal15 holdout
./run_elastic.sh HMS_6p667GeV equal15 enet --check
./run_elastic.sh HMS_6p667GeV equal15 enet
```

The first two commands require the **complete event-level** `equal15` tag and
its replay ROOT inputs on ifarm. An allocation plot directory is not an
allocation manifest. Existing immutable exports are never overwritten; use
existing exports if their checks pass. Older holdout builds without ROOT-file
checksums need a new export using this version of `build_core_fit.py`.
Do not delete or overwrite a prior tag simply to make a command proceed.

The fit uses training events only. Review the training validation curves and
term stability first. If beam search is to follow, do its tuning on these same
training folds before opening protected evaluation. When ready for the final
comparison, run:

```bash
./run_elastic.sh HMS_6p667GeV equal15 enet --evaluate
```

This evaluates the saved fit; it does not refit. Evaluation is saved once and
refuses overwrite. Changing result names cannot make a previously inspected
holdout statistically untouched again.

Outputs go to `HMS_6p667GeV/06d_elastic_net/enet/`. Substitute the campaign,
input tag and short output name for another experiment. No matrix is installed
in replay and no production SVD macro is changed by this command.

## Fit choices

- **Basis and units:** use the non-xtar terms in the frozen training build's
  `oldfit.dat`, plus one constant. HMS 6.667 supplies 210 fitted terms including
  that constant. Positions enter the polynomial in meters and slopes in radians,
  matching the existing fit. Subtract the seed's xtar-dependent contribution
  using the saved replay `xtar` for both fitting and prediction. Preserve all
  fixed xtar coefficients and the complete delta coefficient column at export.
  This avoids mixing replay xtar in fitting with truth xtar in evaluation.
- **Scaling:** center and standardize polynomial columns and each target using
  the training portion of each fold. Exactly constant columns have zero slope;
  their indices are reported. An unpenalized intercept is transformed back into
  the constant transport term. Refit scaling on all training events only after
  choosing regularization. Write native coefficients with a numerical prediction
  round-trip check. No scaler uses protected data.
- **Elastic net:** fit each of the three targets separately with scikit-learn's
  coordinate-descent solver. Its objective is mean squared error divided by two,
  plus alpha times a mixture of L1 and L2 penalties. L1 can remove terms; L2
  stabilizes correlated terms. See the [solver reference](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html).
- **Validation:** three deterministic event folds within each rungroup/local
  foil/delta/hole training population. Hash ordering preserves 64-bit entry IDs
  and makes splits independent of input order. Round-robin spreads tiny cells
  where possible; a singleton cannot occur in every fold. These are within-data
  splits, not tests of independence between runs or time-correlated events.
  Density regions and allocation are already frozen; their estimation is not
  repeated inside the folds.
- **Score:** average validation MSE equally across populated physical-foil/delta
  cells within each fold, then across folds, separately per target. This keeps
  a populous slice from dominating tuning. Hole balance comes from the existing
  allocation; there is no additional inverse-hole weighting. The diagnostics
  separately expose sparse-hole errors. No folding procedure can give a missing
  population a performance estimate.
- **Choice:** among candidates within one standard error of the lowest mean CV
  score, choose the one with the fewest average nonzero slopes, then larger alpha
  and L1 ratio to break ties. The standard error across three folds is a
  practical selection heuristic, not a formal confidence interval. A candidate
  must converge in every fold; the final refit must converge too. Grid-edge
  choices are flagged. No automatic grid expansion or protected-data tuning.
- **Comparison:** fit a centered/scaled direct least-squares SVD baseline on
  identical events, terms, fixed contributions and transformations. Its singular
  values, numerical rank and explicit cutoff are saved. This comparison isolates
  regularization; it does **not** reproduce the historical unscaled ROOT
  normal-equation solve or its older sample.

## Campaign configuration

Optional `CAMPAIGN/config/elastic.json` overrides these defaults. The effective
configuration and package versions are saved with every result; unknown keys
and invalid ranges fail rather than being ignored.

| Key | Default | Meaning |
|---|---|---|
| `seed` | 667 | Deterministic fold assignment |
| `folds` | 3 | Training-only cross-validation folds |
| `alphas` | 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001, 0.00003, 0.00001 | Penalty strengths on standardized columns/targets |
| `l1` | 0.1, 0.5, 0.9 | L1 share of the penalty |
| `tol` | 0.000001 | Coordinate-descent convergence tolerance |
| `max_iter` | 20000 | Iteration limit per fit |
| `rcond` | 0.000000000001 | Scaled SVD relative singular-value cutoff |
| `qa_min` | 10 | Sparse diagnostic warning only; never changes membership |

For example, change only `alphas` to widen a training-only study. Use a new
short result name, such as `enet2`, to preserve the first result. A full run
uses a modest grid and runs serially; it does not launch an unbounded search.

## Diagnostic portfolio

Each figure has matching numerical tables. Counts and residuals describe all
supplied events: there is no plotting subsample or hidden residual-range cut.

| Output | Question answered |
|---|---|
| `RESULTS.md` | What was selected, and how did training validation compare with SVD? |
| `plots/cv.png` | Does the score support the chosen alpha/L1 ratio, or hit a grid edge? |
| `plots/terms.png`, `tsv/terms.tsv` | Which terms survive and how consistently across folds? |
| `tsv/training.tsv` | Did the training error hide a foil, delta or individual setting/hole problem? |
| `evaluation/plots/residuals.png` | Bias and tails for protected cores, supported noncores and unsupported labels, separately |
| `evaluation/plots/foil_delta.png` | Elastic/SVD RMS ratio for every physical foil and delta, including ztar |
| `evaluation/plots/slice_N_holes.png` | Log-colored numbered event counts and angular RMS ratios for core/noncore holes across all foils at that delta |
| `evaluation/plots/slice_N_core_clouds.png` and `slice_N_noncore_clouds.png` | Identical events under saved replay, SVD and elastic-net reconstruction, with common axes across foils |
| `evaluation/tsv/metrics.tsv` | N, bias, RMS and P95 absolute error, overall and by foil/delta, setting, and setting/hole |
| `evaluation/tsv/ids.tsv` | Exact protected event identities and pool classifications, including excluded blocked labels |

The angular hole-map score combines xptar and yptar squared errors in mrad;
ytar and ztar remain separate in the tables and foil/delta plot. RMS-ratio
colors saturate outside 0.5–1.5, with extension indicators; printed numbers
retain the actual ratios. Blank positions have no events, not zero error.
Black hole outlines flag fewer than `qa_min` events. Foil/delta tables preserve
setting-level detail even when the figures aggregate multiple settings.

All protected events remain accounted for. `quality=2` is the core test pool,
`quality=1` the supported noncore test pool, and `quality=0` a separate unsupported
stress test. Blocked labels receive counts but no target residual scores.
Noncore results measure extrapolation away from the selected cores within this
campaign, assuming inherited labels are correct. Scattering or mistaken labels
can produce irreducible tails; these residuals are not all matrix error.

Sieve clouds use the predicted angles/ytar and saved replay xtar; ztar uses the
HMS target identity with beam position inferred from exported truth geometry.
These are fixed-input reconstruction diagnostics, **not full replay**. A candidate
still needs replay/iteration and independent physics checks before adoption.

## Beam-search handoff

The reusable handoff is `tsv/folds.tsv`, `tsv/terms.tsv`, `seed.dat`, `model.npz`,
`manifest.json` and the unchanged balanced training build. The term field is a
five-digit exponent string: preserve its leading zeros when reading the TSV.
Use fold selection frequency and standardized coefficient magnitude to propose
small candidate subsets. Correlated terms can substitute for one another;
never treat elastic-net zeros as proof that a physical term is unnecessary.
The constant, fixed xtar terms and delta coefficients are mandatory.

Before committing to a beam search, check for a reproducible training-validation
advantage or useful reduction in terms without a material score loss. Widespread
unstable membership, alpha at the search boundary, or poor convergence calls for
resolving the elastic-net fit first. Search and compare subsets on training folds
with fold-specific scaling and the same score. Holdouts must not pick the subset,
beam width, stopping rule or penalty. After final protected evaluation, report
both pooled performance and sparse/edge-cell changes; do not hide a local
regression behind a better global RMS. There is no automatic claim of a better
matrix or beam-search readiness.

## Reproducibility and verification status

Input hashes tie the fit to the allocation, geometry, seed and exact TFit
exports. `model.npz` stores scalers and native coefficients; transport files are
written separately for elastic net and SVD. Source snapshots and package versions
are saved. The checks fail on changed inputs, missing identities, blocked/noncore
training labels, unequal physical-foil totals, invalid geometry, or nonfinite
values. They never silently remove events or fall back to historical GMM trees.

Local verification covers a five-foil synthetic ROOT campaign, both delta
slices, every pool, 64-bit identities above 2^54, deterministic folds/refits,
unchanged fits when holdout exports are absent, convergence rejection, coefficient
round-trip checks, fixed/delta preservation, input tampering, and rendered
coverage/cloud diagnostics. **Synthetic verification is not a 6.667 physics result.**

The September 12 local preflight finds these absent from the GitHub mirror:

- `HMS_6p667GeV/05c_core_sample/equal15/manifest.json`
- `HMS_6p667GeV/06c_core_ntuple/equal15/fit/build.json` and the associated export
- `HMS_6p667GeV/06c_core_ntuple/equal15/holdout/build.json` and the associated export

The mirrored equal15 diagnostic figures do not supply those inputs. Run the
commands on the complete ifarm data, or expose the full balanced tag and both
fit/holdout export directories to run the identical study locally. No real
campaign elastic-net result or superiority claim has been produced here yet.

## Bounded convergence study after the first 6.667 result

The first run selected the lowest converged penalties on their paths. Many
weaker penalties did not converge, so the original `alpha_edge=false` field
must not be read as evidence that the search found an interior optimum. The
original CV plot omits failed candidates; use `tsv/cv_summary.tsv` to see them.

Before another full grid, run this smaller study from the repository root:

```bash
./run_elastic.sh HMS_6p667GeV equal15 conv --convergence
```

It reads the saved `enet` result and exactly reuses fold 0 from its `folds.tsv`.
It fits all three targets at alpha 0.0001 and 0.00003, with L1 fractions 0.1 and
0.5: twelve cases on one training/validation split. Each case starts at zero
and continues at the **same** alpha/L1 through cumulative iteration budgets
of 20,000, 100,000 and 300,000, stopping when converged. These budgets are totals,
not extra iterations at each checkpoint. The original tolerance remains 1e-6;
the script reads it from the saved fit and does not offer a tolerance override.
The scaled direct-X SVD reference is recomputed on that same fold.

This is deliberately a controlled continuation study, not an exact repetition
of the original descending-alpha warm-start path. More iterations are the only
change within a given case. Starts can differ from the original first run.
Optional `config/elastic_conv.json` sets the fold, alphas, L1 fractions and
budgets for other campaigns/studies. `--source enet2` uses a different saved
fit if needed; use a new short output name to preserve the first study.

Outputs are in `06d_elastic_net/conv/`:

- `RESULTS.md`: final status and validation error relative to SVD for each case.
- `convergence.png`: convergence gap and validation error versus iterations;
  circles mean converged and crosses mean unfinished.
- `checkpoints.tsv`: actual cumulative iterations, cumulative solver seconds,
  objective, gap relative to its stopping limit, maximum coefficient optimality
  violation (`kkt_max`), active terms, validation MSE, and prediction changes.
  Prediction-change RMS uses mrad/cm/mrad for xptar/ytar/yptar and compares
  consecutive checkpoints on the same validation events; the first is blank.
- `manifest.json` and `code/`: input hashes, effective policy, versions and code.

The maximum optimality violation is an additional diagnostic, not a replacement
stopping rule. Solver time excludes input loading, SVD and plot generation.
Unfinished cases remain visible even when their validation error looks good.
The table is saved after every checkpoint. An interrupted study is not resumed
or overwritten; its partial table can still inform our next decision.

If useful cases converge, repeat those settings across all saved folds before
choosing a matrix. If cases remain unfinished, use the gap/objective/prediction
changes and runtime to decide whether more iterations or a solver change is
justified. Do not loosen tolerance merely to obtain a passing status. This
study neither exports a matrix nor evaluates protected samples.

## Complete convergence and refit the selected terms without penalties

Run this next, using the same exports and original `enet` folds:

```bash
./run_elastic.sh HMS_6p667GeV equal15 refit --refit
```

The study uses the same twelve target/penalty cases on fold 0. It repeats the
short initial passes because the previous `conv` output saved diagnostics but
not resumable coefficient states. Cases stop when converged. Cases that remain
unfinished can continue through cumulative totals of 1,000,000 and 3,000,000
iterations at the original tolerance. It does not relax the stopping rule.
The policy is in `config/elastic_refit.json`; existing `enet` and `conv` results
remain unchanged.

For each converged elastic-net fit, the nonzero coefficients define a proposed
term set. The selected columns of centered/scaled **X are then solved directly
by SVD without penalties**, using the same fitting events. The constant is
always included and is unpenalized. The seed's xtar-dependent angular terms
remain fixed, and delta is not fitted. A case that still fails convergence has
checkpoint records but does not receive a selected-term refit.

Three predictions are compared on identical validation events:

1. `full_svd`: the full original polynomial basis, solved by scaled direct SVD.
2. `enet`: the converged penalized elastic-net prediction.
3. `refit_svd`: the elastic-net-selected basis, solved by scaled direct SVD
   without penalties.

Selection, centering, scaling and all coefficient estimation use **only the
fitting portion of the fold**. Validation events score the resulting predictions;
they do not supply the selected terms or fitted coefficients. There is no
protected-sample evaluation, no change to allocation, and no final matrix is
installed or exported for replay.

Outputs go to `06d_elastic_net/refit/`:

- `RESULTS.md`: convergence status and the three-way comparison.
- `plots/comparison.png`: penalized and refitted validation MSE relative to
  full-basis SVD, with selected term counts including the constant.
- `comparison.tsv`: exact scores, selected term counts, refit rank and
  conditioning. A count of selected terms does not guarantee full numerical
  rank; the same saved SVD cutoff is used and the rank is reported.
- `validation.tsv`: validation N, bias, RMS and P95, overall and by physical
  foil/delta, setting, and setting/hole. It records all three predictions.
- `plots/*_foil_delta.png`: elastic-net and refitted RMS ratios across all
  physical foils and delta slices. ztar appears only for penalty combinations
  where all three targets converged; it combines their predictions using the
  same fixed-input HMS reconstruction diagnostic as before.
- `plots/*_holes.png`: selected-term SVD/full-SVD RMS ratios at sieve holes,
  with all foils and delta slices on each target/case page. Settings are pooled
  by event count for these maps but remain separate in the table. Black outlines
  mark fewer than `qa_min` validation events, and blank positions have no
  validation estimate. The numbers printed in circles are **RMS ratios**, not
  counts; counts are in `validation.tsv`. Colors saturate outside 0.5–1.5 and
  printed numbers preserve the actual ratio.
- `terms.tsv`, `coefficients.npz`, `seed.dat`, `folds.tsv`: term selections,
  native coefficients, fixed terms, and exact saved fold assignments for the
  next study. Coefficients here describe fits on two thirds of the training
  events, not a final fit to the full balanced sample.
- `checkpoints.tsv`, `convergence.png`, `manifest.json`, `code/`: continuation
  diagnostics and reproducibility records as in the earlier study.

The constant is included in the term counts in this comparison, whereas the
original elastic-net report's "active slopes" count excluded it. Five-digit
term codes must be read as strings to preserve leading zeros.

If the unpenalized refit recovers the accuracy that the penalized fit lost,
elastic net is providing useful candidate bases. If the refit remains worse,
the removed terms or their alternatives may still be needed. Inspect local
residuals as well as overall scores before choosing settings for the full
three-fold study. A promising one-fold score is development evidence, not a
final selection or proof of improvement. Beam search comes after that check,
and it must be able to reintroduce terms excluded by elastic net.
