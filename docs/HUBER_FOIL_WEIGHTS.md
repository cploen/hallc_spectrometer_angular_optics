# Equal-foil weights with all cores and supported shoulders

This follow-up adds fixed population weights to the smooth-Huber experiment.
It retains every eligible development event and performs four full-basis fits:
all cores and all cores plus supported shoulders, each with weighted squared
error and weighted smooth Huber. It does not retune the threshold or the basis.

```bash
python3 -m unittest discover -s tests -p 'test_huber*.py' -v
python3 fit_huber_balanced.py HMS_6p667GeV equalfoil_20260924 \
  --source expansion_20260924 --threads 4
```

Output is immutable under `HMS_6p667GeV/06f_huber/<name>/`. Use a new name for
another run. The completed preceding expansion study and original frozen trees
are required; source hashes and event identities are checked before fitting.

## Weight definition

For training event i in physical foil f, set

\[
 w_i=\frac{N}{5N_f},\qquad
 \mathcal L=\frac{\sum_i w_i\rho((X_i\beta-y_i)/\sigma)}{\sum_i w_i}.
\]

Every foil has 20% of the total base weight. Each event remains in the fit;
events from densely populated foils receive smaller individual weights. Foil
counts are taken only from the relevant training pool. Supported shoulders are
included in those counts when training on the broader pool.

Huber uses the same threshold `c = 1.5` and balanced-training MAD scales as the
preceding study. The weights multiply the loss **outside** the Huber function.
Multiplying residuals by sqrt(w) inside Huber would change the physical
transition threshold from foil to foil and is not the intended objective.

Weighted centering/scaling and QR/SVD preconditioning improve numerical
conditioning without changing the objective. The intercept is freely fitted.
Weighted squared error is checked against a direct sqrt(w)-weighted SVD solve;
weighted Huber is independently checked against unweighted fits to duplicated
observations representing integer weights. Existing unweighted tests also pass.

No separate delta-slice or sieve-hole weights are applied in this first test.
The old `equal15` allocation imposed both foil equality and slice/hole caps,
so these population weights do not exactly reproduce its distribution.

## What the contribution diagnostic means

At the converged Huber fit, the residual factor is

\[
 h_i=\left[1+\left(\frac{r_i}{c\sigma}\right)^2\right]^{-1/2}.
\]

The table reports each foil's share of `sum(w_i h_i)`, separately for each
target. Base shares are exactly equal; residual-weighted shares need not be.
There is no second normalization after applying Huber, since that would change
the declared objective. These are IRLS/estimating-equation weight shares, not a
complete measure of coefficient influence or leverage.

`tsv/weight_contributions.tsv` also separates cores from shoulders and gives
their mean residual factor and Kish effective count. `training_weights.npz`
stores exact training IDs and all base weights. `tsv/coverage.tsv` shows how
those weights are distributed over foil/delta cells.

## Evaluation

The same 77,864 reserved cores and 22,532 reserved shoulders are used throughout.
All four candidate fits finish before their reserved predictions are scored.
The output contains the five saved unweighted references plus four new fits.

Read `RESULTS.md` and `plots/comparison.png`. `plots/foil_delta.png` compares the
weighted broad-sample Huber fit with the original balanced SVD.
`plots/weighting_effect.png` compares it with the unweighted Huber fit on
identical training events, isolating the foil-weighting change.

`tsv/summary.tsv` includes pooled RMS and spread, equal-cell RMS, and within-hole
spread. Within-hole spread removes each run/foil/delta/hole group's own mean,
then takes the square root of the event-weighted average variance across common
groups containing at least 10 reserved events. It is not a Gaussian resolution
fit. Sparse groups remain recorded in `tsv/holes.tsv`.

This reserved sample has been examined before, so results remain diagnostic
rather than fresh independent validation. Existing labels, target geometry,
saved xtar and the full 210-term basis are unchanged. Candidate matrices are
exported and prediction-checked; none is installed in replay.
