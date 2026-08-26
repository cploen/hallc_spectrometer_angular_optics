# Per-foil angular residual ladders

Each primary plot shows an independently refitted rescaled direct-`X` solution at every `N=1..210`, separated into selected training events and eligible events excluded by the fit caps. The displayed range begins at `N=35` so that high-order behavior is visible.

- Black dashed line: minimum held-out RMS within the displayed `N>=35` range.
- Red dotted line: first run of five consecutive `N` values more than 5% above that minimum. This is a descriptive marker, not a production cutoff.
- Lower panel: held-out RMS change relative to the complete fifth-order `N=126` solution.

## Rescaled direct-X plots

| nominal foil | xptar | ytar | yptar |
|---:|---|---|---|
| -8 cm | [PNG](rescaled_direct_X/foil_m8cm/xptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_m8cm/xptar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_m8cm/ytar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_m8cm/ytar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_m8cm/yptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_m8cm/yptar_residual_ladder.pdf) |
| -3 cm | [PNG](rescaled_direct_X/foil_m3cm/xptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_m3cm/xptar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_m3cm/ytar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_m3cm/ytar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_m3cm/yptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_m3cm/yptar_residual_ladder.pdf) |
| 0 cm | [PNG](rescaled_direct_X/foil_0cm/xptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_0cm/xptar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_0cm/ytar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_0cm/ytar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_0cm/yptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_0cm/yptar_residual_ladder.pdf) |
| +3 cm | [PNG](rescaled_direct_X/foil_p3cm/xptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_p3cm/xptar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_p3cm/ytar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_p3cm/ytar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_p3cm/yptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_p3cm/yptar_residual_ladder.pdf) |
| +8 cm | [PNG](rescaled_direct_X/foil_p8cm/xptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_p8cm/xptar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_p8cm/ytar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_p8cm/ytar_residual_ladder.pdf) | [PNG](rescaled_direct_X/foil_p8cm/yptar_residual_ladder.png) / [PDF](rescaled_direct_X/foil_p8cm/yptar_residual_ladder.pdf) |

## Current unscaled-XTX plots

These use ROOT `TDecompSVD` on the unscaled normal-equation matrix for every independent prefix solve, matching the current numerical method and its default truncation rule.

| nominal foil | xptar | ytar | yptar |
|---:|---|---|---|
| -8 cm | [PNG](current_XTX/foil_m8cm/xptar_residual_ladder.png) / [PDF](current_XTX/foil_m8cm/xptar_residual_ladder.pdf) | [PNG](current_XTX/foil_m8cm/ytar_residual_ladder.png) / [PDF](current_XTX/foil_m8cm/ytar_residual_ladder.pdf) | [PNG](current_XTX/foil_m8cm/yptar_residual_ladder.png) / [PDF](current_XTX/foil_m8cm/yptar_residual_ladder.pdf) |
| -3 cm | [PNG](current_XTX/foil_m3cm/xptar_residual_ladder.png) / [PDF](current_XTX/foil_m3cm/xptar_residual_ladder.pdf) | [PNG](current_XTX/foil_m3cm/ytar_residual_ladder.png) / [PDF](current_XTX/foil_m3cm/ytar_residual_ladder.pdf) | [PNG](current_XTX/foil_m3cm/yptar_residual_ladder.png) / [PDF](current_XTX/foil_m3cm/yptar_residual_ladder.pdf) |
| 0 cm | [PNG](current_XTX/foil_0cm/xptar_residual_ladder.png) / [PDF](current_XTX/foil_0cm/xptar_residual_ladder.pdf) | [PNG](current_XTX/foil_0cm/ytar_residual_ladder.png) / [PDF](current_XTX/foil_0cm/ytar_residual_ladder.pdf) | [PNG](current_XTX/foil_0cm/yptar_residual_ladder.png) / [PDF](current_XTX/foil_0cm/yptar_residual_ladder.pdf) |
| +3 cm | [PNG](current_XTX/foil_p3cm/xptar_residual_ladder.png) / [PDF](current_XTX/foil_p3cm/xptar_residual_ladder.pdf) | [PNG](current_XTX/foil_p3cm/ytar_residual_ladder.png) / [PDF](current_XTX/foil_p3cm/ytar_residual_ladder.pdf) | [PNG](current_XTX/foil_p3cm/yptar_residual_ladder.png) / [PDF](current_XTX/foil_p3cm/yptar_residual_ladder.pdf) |
| +8 cm | [PNG](current_XTX/foil_p8cm/xptar_residual_ladder.png) / [PDF](current_XTX/foil_p8cm/xptar_residual_ladder.pdf) | [PNG](current_XTX/foil_p8cm/ytar_residual_ladder.png) / [PDF](current_XTX/foil_p8cm/ytar_residual_ladder.pdf) | [PNG](current_XTX/foil_p8cm/yptar_residual_ladder.png) / [PDF](current_XTX/foil_p8cm/yptar_residual_ladder.pdf) |

The unscaled-`XTX` coefficients are from the local ROOT ladder and therefore reproduce the algorithm and truncation rule. They do not reproduce the exact ifarm summation order bit-for-bit; that accumulation-order ambiguity remains documented separately.

## Immediate pattern

For the rescaled direct-`X` ladder, the `+8 cm` foil is the only nominal foil showing a five-consecutive-`N`, 5% held-out degradation in all three fitted quantities:

| quantity | minimum held-out RMS for N>=35 | first five-N run >5% above minimum |
|---|---:|---:|
| xptar | N=61 | N=64 |
| yptar | N=115 | N=120 |
| ytar | N=163 | N=170 |

The most severe excursions occur later than these descriptive onset markers, particularly above roughly `N=180-200`.

The foil aggregation combines all settings that share a nominal foil position. Because earlier diagnostics localized the most extreme behavior to the final `theta=12.495 deg, +/-8 cm` setting, a setting-by-foil breakdown remains necessary before a production decision.

## Tables

- `angular_per_foil_residual_ladder.tsv`: every `N`, method, foil, quantity, count, and RMS.
- `angular_per_foil_overfitting_onsets.tsv`: compact minimum/onset summary.
- `angular_per_foil_residual_ladder_summary.txt`: human-readable summary.
