# Standard matrix diagnostics

Run from the project root after the three beam-20 results are complete:

```bash
bash run_standard_diagnostics.sh --check
bash run_standard_diagnostics.sh
```

Campaigns run sequentially. Append a campaign name to process only that campaign.
New outputs use `07_diagnostics/beam20_<tag>_standard_2024flags/`, so earlier
zero-offset comparisons remain separately identifiable.

## Matrices and offsets

Every protected comparison applies all four complete matrices through the same
fixed-vertex iterative-xtar afterburner on identical reserved event IDs:

1. The 2024 replay matrix plus the historical momentum-specific HMS flags.
2. Full-basis SVD, with fitted constants already embedded.
3. The EN-selected reduced-basis SVD, with fitted constants embedded.
4. The final beam matrix, with fitted constants embedded.

The modern matrices receive zero additional offsets. The 2024 reference receives
these additions, in `[xptar mrad, ytar cm, yptar mrad]` order:

| Campaign | xptar | ytar | yptar |
|---|---:|---:|---:|
| HMS_5p878GeV | 0.660916405 | 0 | 0.453977455 |
| HMS_6p117GeV | 0.581654738 | 0 | 0.647613850 |
| HMS_6p667GeV_joint_min5 | 0.711802209 | 0 | 0.584959117 |

`historical_offsets_2024.json` pins the values, units, source URLs and snapshot
checksums. The flags snapshots are in `references/nps_replay_2024_flags/` and
come from JeffersonLab/nps_replay commit
`82f6fa3bc8aae79298b649a8e83523c6ee39797a`. The loader verifies that xptar is
1000 times `hphi_offset` and yptar is 1000 times `htheta_offset`. No external
YTAR addition is specified by these flags.

The values are applied inside reconstruction, not by recentering residuals.
`tsv/matrices.tsv`, `manifest.json` and the code snapshots record the actual
values and sources. Every comparison table and heatmap uses those predictions.

## Heatmap pages and event populations

Each population has three pages:

- `<pool>_foil_delta.png`: absolute RMS for all four matrices.
- `<pool>_change.png`: percentage RMS change against the offset-corrected 2024
  reference, plus beam versus full SVD and beam versus the compact starting basis.
- `final_beam_<pool>_foil_delta.png`: absolute RMS of the final beam matrix alone.

The populations are:

| Pool | Events |
|---|---|
| `core` | Protected quality-2 cores |
| `noncore` | All protected supported quality-1 shoulders |
| `core_shoulder` | Union of those cores and all supported shoulders |
| `unsupported` | Protected quality-0 stress sample |

Blocked positions have counts only and are excluded from all four scored
populations. The combined page includes all supported shoulders, including
those outside the inner-shoulder training band. It retains their actual event
counts; it is not resampled to the training 90/10 or 80/20 allocation.

Combined statistics are calculated from the concatenated event residuals.
Combined RMS is therefore the square root of the count-weighted mean squared
residual, not the average of the two pool RMS values. Each event appears once
in a combined cell. Original pool identities remain in the event archive and
`tsv/ids.tsv`; the latter also marks combined-page membership.

Angles are in mrad; ytar and derived ztar are in cm. N<10 is flagged and retained.
Empty cells remain blank. RMS includes bias. Mean-subtracted spread, bias,
median and P90 remain separately available in `tsv/residuals.tsv` down to
setting/hole. `tsv/cell_comparisons.tsv` stores exact RMS ratios and percentage
changes. Negative percentage change means smaller RMS; zero means unchanged.

## Fold pages and closure

`plots/reduced_basis_folds.png` and `plots/beam_folds.png` use saved training
validation results. No fit or beam search is rerun. Their modern matrices have
embedded fitted constants, so no historical flags are added to those scores.

`tsv/replay_closure.tsv` separately reconstructs the **zero-offset 2024 reference**
and compares it with the saved zero-offset replay coordinates. It deliberately
does not use the offset-corrected historical predictions from the heatmaps.
Adding the historical flags to the comparison should not change this check.
Afterburner iteration-limit counts for the compared matrices are recorded in
`tsv/afterburner_iterations.tsv`.

These remain fixed-vertex afterburner diagnostics. Momentum and reaction vertex
are not recomputed. The 2024 matrix's original training membership is unknown;
the historical comparison describes behavior on the same protected sample.
