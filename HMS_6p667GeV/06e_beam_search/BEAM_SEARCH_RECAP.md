# HMS 6.667 — compact-seed beam search

Conversation recap, 15 September 2026. Full algorithm and commands:
[BEAM_SEARCH.md](../../docs/BEAM_SEARCH.md).

Start each target from all distinct fold-selected bases at its agreed elastic-net
setting. Do not intersect or average the masks. Every candidate is refitted by
**unpenalized scaled-X SVD** on the same development folds. QR compression
preserves that least-squares problem; the solver does not form X-transpose-X.

Search can add or remove nonconstant terms from the original candidate dictionary.
Keep a fitted constant, fixed xtar-dependent terms and the original delta column.
No higher powers beyond the original dictionary are introduced. Score pooled
MSE equally across populated physical-foil/delta cells, with equal event weight
inside each cell. Sparse holes remain visible but cannot independently veto
a broad improvement. This is not inverse-hole weighting.

| Saved run | Width | Maximum rounds | Threads | Final terms: xptar / ytar / yptar |
|---|---:|---:|---:|---|
| `beam` | 4 | 20 | 8 | 110 / 142 / 78 |
| `beam10` | 10 | 40 | 8 | 111 / 143 / 78 |

The first xptar/ytar searches reached the step limit; the wider run stopped for
no material gain on all targets. It became the frozen comparison candidate.
Wider search costs computation and increases adaptive selection on the same
folds; it is not a guarantee of better independent prediction.

The unchanged policy uses gain 0.1%, patience 3 and slack 0.5%: among full-rank
candidates within 0.5% of the best score, choose the smallest, then break ties
by score and conditioning. Lower condition number is not a hidden search penalty.
The tiny 3–5-term/two-delta plots discussed during implementation were test
fixtures, not the 6.667 campaign results.

## Result and handoff

`beam10/matrices/start.dat`, `svd.dat`, and `beam.dat` are respectively the chosen
compact seed, full basis, and beam basis, each refitted on **all 68,535 training
cores**. `start.dat` is not the penalized elastic-net matrix. Final fitted-design
condition numbers from `beam10/tsv/final_fit.tsv` are:

| Basis | xptar | ytar | yptar |
|---|---:|---:|---:|
| Full 210-term basis | 324,016 | 324,016 | 324,016 |
| Compact seed (95 / 132 / 75 terms) | 3,306 | 9,067 | 1,471 |
| Beam (111 / 143 / 78 terms) | 17,714 | 20,903 | 1,790 |

Beam recovered predictive accuracy relative to the compact starts while giving
up some of their conditioning advantage. It did not beat full-core SVD globally.
Search/fold scores are adaptive development evidence, not an independent test;
the subsequent frozen comparison opened the protected pools. No replay matrix
was installed by this step.

Code: `fit_beam.py`, `beam_search.py`, `beam_diagnostics.py`.
`tests/test_beam.py` checks QR versus direct-X coefficients/scores, add/remove
search, thread determinism, rank/slack selection, seed/matrix round trips and
absence of protected-data use. The large `tsv/candidates.tsv` is the exhaustive
candidate log; `RESULTS.md`, `tsv/search.tsv`, summaries and plots suffice for
routine review, but retain the complete run for provenance.
