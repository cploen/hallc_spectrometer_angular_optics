# Reading the comparison plots

- `*_center.png`: signed residuals near zero. Each bin shows a percentage of
  **all events in that pool**, with common bins for every matrix. The window is
  ±beam P90; the legend in each panel reports the percentage visible. No curve
  is recentered. A narrower, taller peak is not by itself evidence of less bias.
- `*_tails.png`: percentage of events with an absolute residual at least as
  large as the horizontal threshold. For example, 1% at 3 mrad means 1% of the
  pool has |residual| ≥3 mrad. Lower is better. The horizontal axis is linear,
  from zero to the largest P99.9 among the five matrices; the vertical axis is
  logarithmic, labelled in percent with 10%, 1%, and 0.1% guides.
- `*_extremes.png`: a separate full-range tail view, with logarithmic axes.
  This retains the rare extremes beyond the main tail view. Both tail figures
  use the entire pool as denominator; zooming never renormalizes the data.
- `*_foil_delta.png`: absolute RMS, all five matrices, all populated physical
  foils and delta slices. Each target uses a common color scale across matrices.
- `*_change.png`: RMS change, 100 × (new/reference − 1), for each foil/delta
  cell. Blue/negative means smaller RMS; red/positive means larger RMS. Colors
  saturate at ±20%; printed values remain exact to the displayed precision.
  Zero reference RMS or absent cells have no defined percentage and show “—”.

Angular residuals are at the target in **mrad**; ytar and derived ztar are in
**cm**. N is the number of evaluation events, identical for both matrices in
each comparison. An asterisk marks N<10 for caution; it never excludes events.
RMS includes bias and is not a fitted Gaussian resolution. Consult bias and
mean-subtracted spread in the unchanged TSVs as well. The historical old-matrix
offset convention remains provisional; changing the plots does not resolve it.

Protected cores and protected noncore events are separate. Surplus cores are
development coverage. These figures change presentation only, not matrices,
offsets, event membership, or residuals.
