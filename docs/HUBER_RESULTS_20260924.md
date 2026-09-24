# Huber expansion results — 24 September 2026

**More data improves some reconstruction errors, but smooth Huber alone does
not replace foil balancing successfully across the acceptance.** The largest
pooled gains favor the heavily populated central foil. Peripheral foil errors
often increase.

The implementation and experiment are in [HUBER_EXPANSION.md](HUBER_EXPANSION.md).
The [complete numerical report](../HMS_6p667GeV/06f_huber/expansion_20260924/RESULTS.md)
contains all sample sizes, model comparisons, conditional uncertainty intervals,
and provenance. No candidate matrix was installed in replay.

## What was run

Five nested training samples range from 68,535 balanced cores to 396,320
supported development events. All variants use the same 77,864 reserved cores
and 22,532 reserved shoulders. Unsupported and blocked labels are excluded.
The main comparison holds all 210 polynomial terms and the existing geometry
fixed. Smooth Huber uses a threshold of 1.5 training-residual MAD scales, fixed
before evaluation, and no population weights.

Existing elastic-net-selected and beam10 term sets were also refitted on all
cores and all supported events using both losses. These are reuse of completed
selection results, not new elastic-net tuning or a new beam search.

## Main numerical comparison

Each table entry gives RMS on **reserved cores / reserved shoulders**. Lower
is better. Angular units are mrad; ytar is cm.

| Training and fit | Training N | xptar RMS | ytar RMS | yptar RMS |
|---|---:|---:|---:|---:|
| Balanced, squared error | 68,535 | 1.4886 / 4.3187 | 0.13582 / 0.16722 | 1.0292 / 1.6757 |
| All cores, squared error | 310,513 | 1.4485 / 4.3422 | 0.13430 / 0.16594 | 1.0232 / 1.6739 |
| All cores, smooth Huber | 310,513 | 1.4485 / 4.3410 | 0.13435 / 0.16652 | 1.0231 / 1.6753 |
| Cores + shoulders, squared error | 396,320 | 1.5057 / 4.1362 | 0.13458 / 0.16428 | 1.0311 / 1.6237 |
| Cores + shoulders, smooth Huber | 396,320 | 1.4656 / 4.2145 | 0.13454 / 0.16466 | 1.0270 / 1.6402 |

Adding surplus cores supplies most of the pooled core improvement. At the
half-surplus step (189,464 events), xptar RMS is already 1.4496 mrad, close to
the all-core result of 1.4485 mrad. Huber adds almost nothing for core-only
training. At the broadest sample, Huber reduces the core degradation caused by
shoulders, while squared error obtains a larger shoulder improvement.

Relative to balanced squared error, all-supported Huber improves pooled core
RMS by 1.55%, 0.94%, and 0.21% in xptar, ytar, and yptar, respectively. Shoulder
RMS improves by 2.41%, 1.53%, and 2.12%. The conditional bootstrap interval for
the small core yptar change includes zero.

RMS includes bias. For core xptar, the all-supported Huber fit changes the
mean-subtracted standard deviation from **1.4620 to 1.4651 mrad**, a 0.21%
worsening, despite the lower RMS. Its shoulder standard deviation decreases
from 4.3179 to 4.2037 mrad, a 2.65% improvement. Thus the pooled core RMS gain
must not be called a corresponding resolution improvement.

## Why the pooled gains are insufficient

For each foil/delta cell, first compute its MSE, then average the 25 cell MSEs
equally and take the square root. Against the balanced squared-error baseline:

| All-supported Huber | Core equal-cell RMS change | Improved core cells | Shoulder equal-cell RMS change |
|---|---:|---:|---:|
| xptar | **2.04% worse** | 7/25 | 2.16% better |
| ytar | 0.35% worse | 12/25 | 0.55% better |
| yptar | 0.33% worse | 12/25 | 1.77% better |

Core xptar worsens in every delta slice of both ±3 cm foils and the +8 cm foil.
The largest increase is 8.82% at +3 cm, slice 3. All five central-foil slices
improve. Even all-core squared error, which has the best pooled core xptar RMS,
worsens the equal-cell core xptar RMS by 1.74%. This is a change in which
populations influence the fit, not an overall improvement across acceptance.

The reused compact bases do not remove this tradeoff. For all-supported Huber,
reserved core/shoulder xptar RMS is 1.5252/4.2850 mrad with the EN-selected basis
and 1.4728/4.2705 mrad with the beam basis, compared with 1.4656/4.2145 mrad for
the full basis. These are descriptive results of frozen supports, not proof
that a new search could never find a better basis.

## Figures and validation

- [Training-sample expansion](../HMS_6p667GeV/06f_huber/expansion_20260924/plots/expansion.png)
- [Foil/delta changes](../HMS_6p667GeV/06f_huber/expansion_20260924/plots/foil_delta.png)
- [Elastic-net-selected and beam basis controls](../HMS_6p667GeV/06f_huber/expansion_20260924/plots/bases.png)

Five targeted tests pass, including analytic-gradient checks, agreement with
an independent SciPy robust solver, ill-conditioned and rank-deficient designs,
outlier resistance, and sample isolation. Input hashes and exact training IDs
match the frozen campaign. The rebuilt balanced SVD reproduces archived
predictions within 8.1e-10 in the reported physical units. Every robust fit
converged, and every exported matrix reproduces its fitted predictions.

These are fixed-input reconstruction diagnostics using existing target
equations and saved xtar. The reserved sample was inspected in earlier work;
this is not fresh independent validation or a full iterative replay. The
results do not justify replacing the current balanced fit. If the next study
uses all development events, preserving explicit foil/slice population weights
while varying the residual loss would separate the benefit of extra statistics
from the population shift observed here.
