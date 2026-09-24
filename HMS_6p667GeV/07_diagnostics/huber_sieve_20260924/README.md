# Before/after Huber sieve plots

Before: balanced full-basis SVD (68,535 training events). After: smooth Huber on all 396,320 supported development events. Every panel uses identical reserved cores plus shoulders; no fitting is repeated.

[All foils, all delta slices combined](plots/all_foils_all_delta.png) · [Browser gallery](gallery.html)

| Foil (cm) | Slice 0 | Slice 1 | Slice 2 | Slice 3 | Slice 4 | All slices |
|---|---|---|---|---|---|---|
| -8 | [0](plots/foil_m8cm_delta_0.png) | [1](plots/foil_m8cm_delta_1.png) | [2](plots/foil_m8cm_delta_2.png) | [3](plots/foil_m8cm_delta_3.png) | [4](plots/foil_m8cm_delta_4.png) | [combined](plots/foil_m8cm_all_delta.png) |
| -3 | [0](plots/foil_m3cm_delta_0.png) | [1](plots/foil_m3cm_delta_1.png) | [2](plots/foil_m3cm_delta_2.png) | [3](plots/foil_m3cm_delta_3.png) | [4](plots/foil_m3cm_delta_4.png) | [combined](plots/foil_m3cm_all_delta.png) |
| +0 | [0](plots/foil_0cm_delta_0.png) | [1](plots/foil_0cm_delta_1.png) | [2](plots/foil_0cm_delta_2.png) | [3](plots/foil_0cm_delta_3.png) | [4](plots/foil_0cm_delta_4.png) | [combined](plots/foil_0cm_all_delta.png) |
| +3 | [0](plots/foil_p3cm_delta_0.png) | [1](plots/foil_p3cm_delta_1.png) | [2](plots/foil_p3cm_delta_2.png) | [3](plots/foil_p3cm_delta_3.png) | [4](plots/foil_p3cm_delta_4.png) | [combined](plots/foil_p3cm_all_delta.png) |
| +8 | [0](plots/foil_p8cm_delta_0.png) | [1](plots/foil_p8cm_delta_1.png) | [2](plots/foil_p8cm_delta_2.png) | [3](plots/foil_p8cm_delta_3.png) | [4](plots/foil_p8cm_delta_4.png) | [combined](plots/foil_p8cm_all_delta.png) |

All-foil overview for each slice: [0](plots/all_foils_delta_0.png), [1](plots/all_foils_delta_1.png), [2](plots/all_foils_delta_2.png), [3](plots/all_foils_delta_3.png), [4](plots/all_foils_delta_4.png).

Delta boundaries (%): -10.0, -8.0, -5.0, 0.0, 5.0, 10.0.

Identical reserved cores + shoulders in both panels; run groups pooled by physical foil.
Saved xtar; existing sieve projection. Same axes, bins and log color scale throughout.

Bin sizes are 0.05 cm in Y and 0.08 cm in X. Outside-frame counts are printed on each panel and saved in `counts.tsv`. Coordinates follow the existing `elastic_diagnostics.py` sieve-cloud projection; this is not a fresh iterative replay. This before/after comparison includes both sample expansion and the loss change.
