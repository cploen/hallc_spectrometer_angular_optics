# Sieve afterburner

Make GMM-versus-beam sieve plots directly from full replay trees. No new tracking,
labeling, GMM selection, core selection or fit is run. The default 6.667 setup
selects exactly 12.490 degrees: rg01/run 1544 (0 cm) and rg03/run 1540 (−8/+8 cm).
It produces three PNGs per matrix, six total.

```bash
./run_sieve.sh HMS_6p667GeV --check
./run_sieve.sh HMS_6p667GeV
```

Outputs: `HMS_6p667GeV/07_diagnostics/sieve/`. Existing outputs are never
replaced. Choose another short name to repeat:

```bash
./run_sieve.sh HMS_6p667GeV sieve2 --source beam10 --threads 8
./run_sieve.sh HMS_6p667GeV sieve15 --angle 15.195
```

The second example selects all configured groups at 15.195 degrees, with
separate plots per rungroup/foil. The angle comparison is exact to 1e-6 degrees;
12.495 groups are not mixed into the requested 12.490 sample.

## Inputs and cuts

The campaign rungroup TSV supplies the existing ifarm replay paths. Foil
positions come from the existing optics-ID metadata. The GMM matrix comes from
`config/comparison.json`; beam defaults to
`06e_beam_search/beam10/matrices/beam.dat`. Optional `--gmm path` and `--beam path`
let this tool compare other saved matrices without changing campaign defaults.
Custom relative paths are relative to the current directory.

`config/sieve.json` holds the defaults; `--angle`, `--delta LOW HIGH`, `--cer`,
`--cal`, `--foil-width` and `--source` override them. Default electron cuts are
H.cer.npeSum > 2 and H.cal.etottracknorm > 0.65, following the existing ztar
resolution-by-delta diagnostic. Delta is the original replay H.gtr.dp, strictly
between −10 and +10 percent. Foils use the original H.react.z within ±2 cm of
nominal. Foil windows must not overlap. These identical cuts and assignments
are made once, before either matrix is applied. No sieve-hole masks are used.
Invalid/nonfinite replay rows are counted and excluded before selection.
These plots describe the fixed original electron/foil selection, including
its tails and possible foil contamination; they do not promise pure truth foils.

## Reconstruction and offsets

This is a **fixed-vertex afterburner**, not a complete replay. It keeps the
saved golden-track focal-plane coordinates, reaction y/z and delta selection.
For each matrix it repeats the extended-target correction from
[THcExtTarCor](https://github.com/JeffersonLab/hcana/blob/master/src/THcExtTarCor.cxx):

1. Evaluate the angular/ytar polynomial with x_tg = −react.y − xmis.
2. Set xtar = x_tg − xptar × react.z × cos(theta), then reevaluate.
3. Repeat until the xptar change is at most 2 mrad, with at least one update
   and at most five, following the upstream implementation.
4. Project xsieve = xtar + 168 × xptar and ysieve = 100 × ytar + 168 × yptar,
   with target slopes in rad, matrix ytar in m and sieve positions in cm.

The xmis default is the HMS angle-dependent formula used by HCANA and the
repository; `xmis` in sieve.json can specify an explicit cm value. The default
focal-plane convention has zero extra focus/detector/angular transformations,
consistent with the polynomial solves. A differently configured replay must
be checked rather than assumed equivalent. Reaction z is not recomputed from
the new matrix, tracking is not rerun, and momentum is not refitted. A full
replay may therefore change vertex and event selection beyond this comparison.
SHMS is explicitly rejected until its corresponding adapter is implemented.

All valid angular matrix terms, including xtar-dependent terms and constants,
are used. The known corrupt historical GMM row uses the same explicit,
checksum-pinned exclusion as the matrix-comparison step. The original files
are unchanged. Supplying a custom GMM file does not inherit this exclusion.

The configured source directory identifies these as zero-offset replays, but
the exact matrix, hmsflags and HCANA version used to produce them have not been
verified. The closure reference uses config/oldfit.dat with zero external offsets,
not the legacy corrections used in the residual comparison. GMM and beam already
include fitted constants and also receive no external addition by default. Explicit `replay_offsets`,
`gmm_offsets`, `beam_offsets` use [xptar mrad, ytar cm, yptar mrad]. Settings are
saved in the output manifest. No offset is inferred or fitted to the data.

## Replay check and products

The first actual farm run failed xsieve closure for rg03 (both outer foils
combined). See the [campaign closure recap](../HMS_6p667GeV/07_diagnostics/SIEVE_GEOMETRY_RECAP.md)
for the reported numbers and outstanding provenance checks; synthetic closure
tests did not establish agreement with the historical replay.

`--check` inspects branches and reconstructs the first 50,000 input entries per
rungroup, applying the cuts above. It prints xsieve/ysieve RMS differences from
the saved replay coordinates and writes no output. This first-entry sample is
a compatibility check, not a statistical survey. A full run computes closure
on every selected event and saves `tsv/closure.tsv` (bias, RMS, max differences,
invalid-row count, and number of events still above the iteration threshold at
five iterations). The default closure tolerance is 0.05 cm RMS per coordinate.
It is an engineering check, not a resolution requirement. A mismatch prints a
warning and marks every plot as needing replay-check review; it does not stop
plot production. Do not interpret such plots as validated reconstruction until
the mismatch is understood. A successful check does not validate effects omitted
from this fixed-vertex approach.

`plots/` contains six PNGs with Y sieve horizontal and X sieve vertical, common
0.05 cm bins, equal physical aspect and a single log color scale across all
plots. Each title reports foil, angle, delta, N and events outside the frame.
`tsv/counts.tsv` records those counts. `histograms.npz` retains the plotted bin
counts for inexpensive replotting. No large event-level ROOT export is needed.
The manifest records settings, matrix/config hashes, and replay path, size,
mtime and entry count; large input ROOT files are not fully hashed.

The tool processes 50,000 entries at a time, prints progress for each chunk and
precomputes the focal-plane portion of the polynomial before iterating xtar.
Threads default to 8; this does not promise eight-core utilization throughout.
There is no expensive per-hole statistics pass.

For review, expose `plots/*.png`, `tsv/*.tsv`, and `manifest.json`. Keep the
histogram archive on ifarm for future changes to plot presentation. Other HMS
campaigns use their own rungroup table, metadata, comparison matrix paths and
verified offset settings; do not assume 6.667 offsets or angles apply elsewhere.
