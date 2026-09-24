# Results: retaining events with equal-foil weights

**Equal-foil weighting lets us retain all 310,513 development cores with small,
broadly distributed improvements. Including all supported shoulders still
trades core angular accuracy for better shoulder reconstruction.**

[Full numerical report](../HMS_6p667GeV/06f_huber/equalfoil_20260924/RESULTS.md) ·
[Comparison figure](../HMS_6p667GeV/06f_huber/equalfoil_20260924/plots/comparison.png) ·
[Implementation and reproduction](HUBER_FOIL_WEIGHTS.md).

The four planned fits were run: all cores or cores plus supported shoulders,
each with weighted squared error and weighted smooth Huber. Every foil starts
with 20% of the total weight. Basis, geometry, residual scales and Huber
threshold remain as in the previous test; evaluation uses the same 100,396
reserved events. No matrix was installed in replay.

## Statistical balancing works as intended

All 310,513 available development cores remain in the core-only fits, compared
with 68,535 in the original balanced allocation. Adding 85,807 supported
shoulders gives 396,320 training events. No balancing subsample is drawn.

After Huber's residual reduction, each foil carries between **19.67% and 20.27%**
of the total estimating-equation weight across the broader fit's three targets.
The central foil no longer dominates through its event count.

Shoulders comprise 21.65% of events and 21.56% of total base weight. Their
contributions after Huber weighting are:

| Target | Shoulder share after Huber | Mean shoulder residual factor |
|---|---:|---:|
| xptar | 16.90% | 0.648 |
| ytar | 20.81% | 0.837 |
| yptar | 19.64% | 0.777 |

These are sums of base weight times the Huber residual factor. They measure
weights in the fitting equation, not the full influence of events on individual
coefficients, which also depends on their focal-plane coordinates.

## Reconstruction with all cores

Relative to the original balanced SVD, equal-foil Huber gives:

| Reserved core metric | xptar | ytar | yptar |
|---|---:|---:|---:|
| Pooled RMS reduction | 0.303% | 0.606% | 0.163% |
| Equal foil/delta-cell RMS reduction | 0.277% | 0.428% | 0.345% |
| Improved foil/delta cells | 18/25 | 18/25 | 18/25 |
| Within-hole width reduction | 0.222% | 0.676% | 0.770% |

The gains are modest. Conditional 95% cluster-bootstrap intervals for pooled
core RMS gains are [-0.046, 0.588]% for xptar, [0.420, 0.789]% for ytar, and
[-0.049, 0.356]% for yptar. The tiny angular gains are not resolved from zero by
this diagnostic. Weighted squared error performs similarly, so Huber adds little
in core-only training.

This improves the population tradeoff from the unweighted all-core study:
equal-cell core xptar RMS previously worsened by about 1.7%; it now improves
slightly. The large earlier pooled gain was not preserved, because it favored
the high-statistics central foil.

## Reconstruction when shoulders are included

With equal-foil weights and Huber, the broader sample changes RMS as follows:

| Reserved population | xptar | ytar | yptar |
|---|---:|---:|---:|
| Cores | **2.016% worse** | 0.296% better | **0.798% worse** |
| Shoulders | 2.611% better | 1.335% better | 3.219% better |

Shoulder RMS improves in 23/25, 22/25 and 22/25 foil/delta cells respectively.
For cores, equal-cell xptar RMS worsens by 1.192%, with only 8/25 cells improved.
Weighted squared error on the broader sample worsens core xptar RMS by 5.735%;
Huber reduces that degradation but does not eliminate it.

The width/centroid decomposition explains why a visually narrower hole need
not reconstruct more accurately. For core holes containing at least ten reserved
events, all-supported Huber changes the event-weighted xptar within-hole width
from **1.3042 to 1.2901 mrad**, while the RMS of their mean residuals increases
from **0.6746 to 0.7604 mrad**. The local distributions narrow, but their
centers move farther from the existing target labels on average. This is a
comparison against inherited labels, not an independent diagnosis of geometry.

## What this establishes

Retaining all cores with foil weights is a useful candidate, with small gains
and approximately equal foil contributions. The broad shoulder-inclusive fit
is not an across-the-board replacement. Residual weighting reduces shoulder
influence but does not isolate their useful information without a core tradeoff.

Nine targeted tests passed, including weighted-SVD parity, analytic gradients,
integer-weight equivalence to duplicated observations, scale invariance and
foil accounting. All full-data optimizations converged and matrix exports
reproduce predictions. These remain diagnostics on previously examined reserved
events, not fresh independent validation or iterative replay.

## Follow-up status after review

The equal-foil results above use the full 210-term basis. The preceding
unweighted expansion study also refitted the saved, independently selected
elastic-net and beam-search bases for each target. Those reduced bases have
not yet been combined with equal-foil weights in this study. Numerical
preconditioning does not remove the underlying full-basis coefficient
sensitivity.

Shoulder-inclusive fits are set aside: some shoulder events may be electrons
scattered by the sieve, for which the nominal hole label is not a suitable
fit target. Huber residual weighting cannot identify that population.
The proposed next comparison is core-only, using the existing per-target
beam-search bases with equal-foil weights, once with squared loss and once
with smooth Huber. It has not been run; no new basis search is needed for
that comparison.

Equal-foil weighting retains all eligible cores but does not balance holes
within a foil. Controlling total base weight per foil/delta/hole group was
discussed as a way to retain more events without allowing high-count holes
to dominate; it is not implemented in these results.

## GitHub artifact availability

Both dated Huber studies include reports, compact summary tables, plots,
exported matrices, manifests, and source snapshots in Git. The sieve gallery
includes all 36 PNG comparisons and binned histograms. Those sieve comparisons
are from the initial unweighted expansion study, not the equal-foil follow-up.

Raw input trees and large event-level outputs remain local: `residuals.npz`,
`training_weights.npz`, `tsv/holes.tsv`, and `tsv/membership.tsv` are not
committed. Manifests describe the complete original runs, including local
files; a manifest entry does not imply that file is distributed in Git.
Reproduction requires the frozen local inputs described in the method guides.
