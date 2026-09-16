# HMS 6.667 — elastic-net studies and beam-seed choice

Conversation recap, 15 September 2026. The complete runbook is
[ELASTIC_NET.md](../../docs/ELASTIC_NET.md); each saved study has its own
`RESULTS.md`, manifest, tables and code snapshot.

## What changed during the study

| Saved name | Question and adjustment |
|---|---|
| `enet` | Initial 3-fold grid, 20,000 iterations. Several weaker penalties did not converge; automated selection among converged candidates was therefore not the final basis decision. |
| `conv` | Controlled continuation on saved fold 0, for alpha 1e-4/3e-5 and L1 fraction 0.1/0.5. Test whether more iterations resolve the solver issue without changing tolerance. |
| `refit` | Extend cumulative budgets to 1 million and 3 million at tolerance 1e-6; all 12 target/penalty cases converged. Refit each selected basis by unpenalized scaled direct-X SVD. |
| `pooled` | Repeat those cases across all three saved folds; all 12 cases completed every fold. Pool out-of-fold residuals to use the whole training sample in the diagnostics. |

`conv`, `refit`, and `pooled` read `enet`'s frozen folds. They do not successively
resume coefficients from the preceding output folder. Single-case continuation
holds the objective fixed; budgets are cumulative iteration totals, not extra
iterations. The hardest reported fold-0 case took 1,081,997 iterations. A
gap/tolerance ratio measures the solver gap against its stopping limit; use the
saved `converged` flag rather than rounded console values. Convergence does not
establish that the chosen basis is physically adequate.

## Statistical and numerical choices

- The 68,535 balanced training cores are partitioned into three folds. Each fit
  uses approximately two thirds, and predicts the remaining third. The pooled
  study gives every training event one validation prediction. Protected cores
  were not borrowed to enlarge a fold. Pooling improves coverage of the diagnostic
  but cannot manufacture information in sparse holes.
- Scaling, term selection and fitted coefficients use fitting events only.
  There are 210 adjustable terms including the unpenalized constant; the seed's
  xtar-dependent contribution remains fixed and uses saved replay xtar.
- Compare the full and selected bases on the **same rows, scaling convention and
  direct-X SVD cutoff (1e-12)**. Condition numbers describe the centered/scaled
  nonconstant design block, not X-transpose-X or a matrix of residuals. They
  depend on the sample and should not be compared indiscriminately with the
  earlier historical conditioning study.
- Elastic net supplies candidate terms. The intended final coefficients come
  from **unpenalized SVD**, so penalized prediction error is not the deciding
  result. There was no requirement that a compact starting basis beat full SVD.
- Keep mean foil/delta MSE, pooled median absolute residual, P90 and conditioning
  together. Median alone can hide tails; MSE alone can overemphasize them. P90 is
  the absolute residual below which 90% of events fall. Angles are target mrad;
  ytar is cm; MSE uses squared units. A plotted MSE of 2.49 mrad² corresponds to
  RMS sqrt(2.49) ≈1.58 mrad, not 2.49 mrad resolution.

The four candidate settings establish a coarse trend, not a precisely located
plateau. The user deliberately chose compact starts and allowed beam search to
add terms back: xptar and ytar **red** (alpha 3e-5, L1 0.5), yptar **green**
(alpha 1e-4, L1 0.5). Use those numerical settings, not colors, for reproduction.
Actual bases differ between folds; an average term count is not a seed matrix.

## Handoff and tests

`pooled/{terms.tsv,coefficients.npz,seed.dat,folds.tsv,manifest.json}` supplies
the next step. `pooled/plots/plateau.png` and `foil_delta.png` summarize the review;
`coverage.tsv` retains sparse-population flags. Final seeds/coefficients after
full-training refits belong to the beam step, not these two-thirds fits.

Code: `fit_elastic.py`, `elastic_net.py`, `elastic_convergence.py`,
`elastic_refit.py`, `elastic_pool.py`, `elastic_diagnostics.py`.
Tests: `tests/test_elastic*.py` exercise fold/scaling separation, nonconvergence
rejection, fixed-objective continuation, selected-term SVD/units, no validation
leakage into refits, and complete versus incomplete pooled cases.

The preliminary conditioning study under `08_preliminary_conditioning/` used a
different, GMM-based sample and held-out definition. Its full-basis residual
growth motivated this work; the stronger modern full-core baseline does not
contradict that result. This sequence did not isolate selection, balancing and
solver changes in a controlled factorial comparison. The historical ROOT cutoff
discarded singular directions; that is not equivalent to deleting the same
number of named polynomial functions from the basis.
