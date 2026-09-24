# Smooth Huber with a larger training sample

Run the actual HMS 6.667 GeV study from the repository root:

```bash
python3 -m pip install -r requirements_core.txt
python3 fit_huber.py HMS_6p667GeV --check
python3 fit_huber.py HMS_6p667GeV expansion --threads 4
```

Outputs are immutable under `HMS_6p667GeV/06f_huber/<name>/`. A failed run retains
its `status: running` manifest and partial files for inspection; use a new name
after correcting the cause. No existing fit output or replay configuration is
changed. `--check` verifies inputs, labels, balanced IDs and pool separation.

## Model

For each angular reconstruction target, minimize

\[
 \frac1N\sum_i c^2\left[\sqrt{1+\left(\frac{X_i\beta-y_i}{c\sigma}\right)^2}-1\right].
\]

Here `c = 1.5`, and `sigma = 1.4826022 MAD(residual)` is estimated separately
for each target from the original balanced training SVD fit. The same scales
and threshold apply to all later fits. There is no foil, momentum-slice, or
hole-population weighting. Consequently, robust loss is an **alternative to the
balancing procedure in this experiment**, not a guarantee of equal foil influence.

This uses the continuous regression
[pseudo-Huber loss](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pseudo_huber.html).
The classification loss named `modified_huber` has different semantics and is
not appropriate for these continuous targets.

The polynomial design is centered/scaled using its own training sample, then
preconditioned with QR/SVD. Optimization in the orthogonal coordinates changes
neither the basis nor the objective. Coefficients are transformed back to the
native optical variables. The intercept is fitted explicitly; frozen xtar and
delta coefficients remain in every exported matrix. Fits require optimizer
success and a maximum gradient below `2e-6`. Rank, conditioning, gradient,
iterations, residual scale, and effective residual weights are recorded.

## Experiment and existing elastic-net/beam work

Five nested sample definitions are fixed before evaluation:

1. Exact archived equal15 balanced training IDs.
2. Those IDs plus approximately half the surplus cores, selected by event hash.
3. All development cores, with no balancing or caps.
4. All cores plus approximately half the supported shoulder events.
5. All supported development events (core and shoulder).

Each receives a squared-error and a smooth Huber fit with the full 210-term
basis. This isolates sample expansion and residual loss. Reusing the completed
elastic-net-selected and beam10 term sets supplies additional controls: both
losses are refitted on all cores and all supported events, with target-specific
supports frozen. This does **not** run new elastic-net penalty tuning or a new
beam search. The archived EN-selected matrix is an unpenalized SVD refit, not
the penalized elastic-net prediction.

No model or threshold is selected from this experiment's reserved scores.
Adding cores mainly increases the event count in already supported regions;
adding shoulders broadens the within-hole distribution. Coverage counts and
hole-level residuals document the distinction. This study does not expand past
the input trees' original cuts, unidentified holes, or unsupported labels.

## Inputs and identity checks

The seven frozen `05c_core_sample/min10/root/CoreSample_*.root` trees contain
the saved focal-plane coordinates, xtar, beam coordinates, and labels needed
to reproduce the HMS TFit target equations locally. Original replay files are
not necessary for this fixed-input experiment. The adapter mirrors the existing
HMS `targetTruth` equations; it makes no geometry correction.

All tree hashes are checked against the frozen manifest. Metadata, physical
foil/delta labels, run identity, hole bounds, unique integer IDs and finiteness
are checked. `equal15/tsv/selected_ids.tsv` must exactly match the checksummed
archived beam training IDs. The seed, archived SVD, EN-selected and beam matrices
are checked against the beam manifest. Before any expanded fit, the local
balanced SVD must reproduce archived SVD predictions to `1e-5` in the reported
physical units. Every exported matrix is reread and compared with predictions.

`sample == 2` remains reserved throughout. Blocked holes and unsupported labels
never train or define the primary residuals. This sample has previously been
examined in the campaign, so its results are diagnostic and cannot be advertised
as fresh independent test performance. All candidate variants are specified
before scoring it. Shoulder scores assume the inherited labels remain valid.

## Reading the result

Start with `RESULTS.md` and `plots/expansion.png`. `plots/foil_delta.png` checks
whether gains persist across physical foils and momentum slices;
`plots/bases.png` checks reuse of compact bases. Positive percent gain means
lower RMS relative to the balanced squared-error full basis.

`tsv/metrics.tsv` separates bias, standard deviation, RMS, central 68% half-width,
and 95th absolute residual. The report also includes equal-cell-weighted RMS
and the number of improved cells. `tsv/holes.tsv` flags groups below ten events.
`tsv/paired_intervals.tsv` contains conditional 95% intervals from 500 paired
bootstrap resamples of run/foil/delta/hole clusters within physical foil/delta
strata. These intervals omit training uncertainty and systematic errors.

This is a common-input optical reconstruction comparison, with saved replay
xtar and existing target conventions. It is not a new iterative HCANA replay,
an independent test of geometry, or a measurement of Gaussian resolution.

Targeted solver and sample-isolation checks:

```bash
python3 -m unittest discover -s tests -p test_huber.py -v
```
