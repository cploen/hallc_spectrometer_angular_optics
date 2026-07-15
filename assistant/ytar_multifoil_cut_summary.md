# Hall C / NPS-HMS Optics: Final `ytar` Foil-Ridge Selection Summary

## To Run

$ hcana -l
$ .x ridge_width_central_guard_multifoil_cut.C(1544,"multifoil_noexpert","-1","none")

## Purpose

This note summarizes the final direction of the `ytar` foil-selection work so it can be carried into a fresh chat before moving on to sieve-hole ownership / selection.

The current goal was to replace hand-drawn `ytar` vs. `delta` foil cuts with an automated, data-driven method that still respects the known Hall C target geometry.

The final working script is conceptually:

```text
ridge_width_central_guard_multifoil_cut.C
```

It detects one or more foil ridges, assigns foil indices, builds one `TCutG` per foil, and overlays expert hand cuts only for benchmarking.

---

## Core principle

The final method is **not** a global density cut.

Instead, it treats each foil as a coherent ridge in `ytar` that can move smoothly with `delta`.

The key idea is:

```text
Find the foil ridge center from the data.
Measure the ridge width from the data.
Use the best central/high-density region to prevent unphysical widening elsewhere.
```

This keeps the method mostly unsupervised while adding a simple geometry-informed guardrail:

```text
The target is a foil, not a wedge.
```

---

## Final method label

Recommended method language:

```text
Data-adaptive, central-width-guarded ridge-envelope selection
```

Longer paper-friendly phrasing:

```text
Foil ridges are identified from the central-delta ytar density distribution and indexed by increasing ytar. For each ridge, a smooth ytar-vs-delta centerline is estimated, local widths are measured from density falloff, and a run-specific central-width guardrail prevents low-density tails or boot-like structures from producing unphysical envelope widening.
```

---

## Production event gate

The event population starts with the same simple production-style cuts used throughout this workflow:

```cpp
H.cer.npeSum > 2.0
-10 < H.gtr.dp < 10
finite ytar
finite delta
```

In code terms:

```cpp
if (!(sumnpe > 2.0)) continue;
if (!(delta > deltaMin && delta < deltaMax)) continue;
if (!std::isfinite(ytar)) continue;
if (!std::isfinite(delta)) continue;
```

This defines the base event population for both ridge detection and final acceptance testing.

---

## Ridge detection for 1, 2, or 3 foils

The multifoil version first builds a 1D `ytar` density histogram using events in the central delta region:

```cpp
abs(delta) < centralDeltaAbsMax
```

Current central region:

```cpp
const double centralDeltaAbsMax = 1.5;
```

Then it smooths the central `ytar` distribution and finds dominant local maxima.

Current dominant-ridge criteria:

```cpp
const double minPeakFracOfMax = 0.20;
const double minPeakSeparation = 0.20;
```

Interpretation:

```text
A detected foil ridge must be at least 20% of the tallest central-delta peak
and separated from neighboring detected ridges by at least 0.20 in ytar.
```

This is intentionally simple because Hall C foil ridges are expected to be obvious. If many weird peaks appear, the run should be flagged rather than over-interpreted.

---

## Foil indexing convention

Detected ridge peaks are sorted by increasing `ytar`.

Indexing convention:

```text
foil 0 = most negative ytar ridge
foil 1 = next ridge
foil 2 = highest ytar ridge
```

The output cut names follow this convention:

```text
delta_vs_ytar_cut_foil0
delta_vs_ytar_cut_foil1
delta_vs_ytar_cut_foil2
```

This is important for later sieve-hole selection because each event can inherit foil ownership from the `ytar` cut it passes.

---

## Per-ridge independence

Each detected foil ridge gets its own data-driven parameters.

The following are computed separately for each ridge:

```text
ridge center vs delta
raw left/right widths
smoothed left/right widths
central reference left/right widths
max allowed left/right widths
final TCutG
event counts
expert-overlap counts
```

No central width, peak height, or final width cap is shared between ridges.

---

## How the ridge center is tracked

For each detected foil, the script builds local `ytar` histograms in overlapping delta windows.

Current settings:

```cpp
const double deltaStep = 0.25;
const double deltaHalfWindow = 0.35;
```

Meaning:

```text
Every 0.25 in delta, build a local ytar histogram using events within ±0.35 delta.
```

Within each local delta window, the script searches near the previous ridge center. This prevents ridge identity from swapping between foils.

For multifoil detection, each ridge also receives a `searchHalfWidth` based on nearest-neighbor foil separation:

```text
searchHalfWidth ≈ 0.45 × nearest ridge separation
```

with bounds:

```text
minimum 0.18
maximum 0.75
```

For a single detected ridge, the default search half-width is:

```text
0.75
```

---

## Width definition: `widthFracLevel`

The ridge width is measured using a local density-fraction crossing.

Current setting:

```cpp
const double widthFracLevel = 0.25;
```

For each local delta row:

```cpp
target = widthFracLevel * peakVal;
```

Then the code walks left and right from the local ridge peak and finds where the smoothed local `ytar` density falls to 25% of the local peak.

Conceptually:

```text
left width  = distance from ridge peak to 25%-of-peak crossing on the left
right width = distance from ridge peak to 25%-of-peak crossing on the right
```

This is partly data-driven and partly a method choice:

```text
Data-driven:
  peakVal, crossing location, and width are measured from the data.

Chosen criterion:
  0.25 defines what we call the ridge edge.
```

The value `0.25` was chosen empirically because the default auto cut looked excellent and matched the expert cut very well near central delta.

---

## Central-width guardrail

The important final improvement was replacing a fixed hard-coded `maxWidth` with a run-derived, ridge-specific central-width guardrail.

For each ridge, the script finds a central/high-density subset of rows:

```cpp
abs(deltaCenter) < centralDeltaAbsMax
peakVal >= minPeakFracOfGlobal * globalPeak
```

Current settings:

```cpp
const double centralDeltaAbsMax = 1.5;
const double minPeakFracOfGlobal = 0.75;
const double guardScale = 1.25;
```

For each ridge, it computes:

```text
refLeft  = median raw left width in central/high-density region
refRight = median raw right width in central/high-density region
```

Then:

```text
maxLeft  = guardScale × refLeft
maxRight = guardScale × refRight
```

With the current `guardScale = 1.25`, the maximum allowed width is:

```text
125% of the central reference width for that ridge
```

This is ridge-specific and run-derived.

If too few rows pass the high-density requirement, the code falls back to using central rows based only on:

```text
abs(deltaCenter) < centralDeltaAbsMax
```

---

## Why the guardrail is defensible

The guardrail is not a hand cut.

It is a geometry-informed regularization:

```text
The central/high-density region tells us the measured foil width.
The foil should not become arbitrarily wider at large |delta|.
```

This prevents low-density tails, boot-like features, or poor separation from causing unphysical acceptance growth.

Important observation from testing:

```text
Even when guardScale was loosened dramatically, e.g. 5.0, the auto cut barely exceeded the expert cut.
```

Interpretation:

```text
The local density-width criterion is already doing most of the work.
The guardrail is mostly a safety net, not an active hand-tuned constraint.
```

---

## Final boundary construction

For each foil/ridge:

1. Find raw ridge center and raw left/right widths in each delta row.
2. Smooth the ridge center and widths across delta.
3. Cap the smoothed widths using that ridge’s central-width guardrail.
4. Build boundaries:

```cpp
left(delta)  = ridgeCenter(delta) - widthLeftFinal(delta)
right(delta) = ridgeCenter(delta) + widthRightFinal(delta)
```

5. Build a closed `TCutG` polygon:

```text
left boundary from low delta to high delta
right boundary from high delta back to low delta
```

---

## Final event acceptance

An event is accepted for a given foil if:

```cpp
autoCut_foilN->IsInside(ytar, delta)
```

after passing the base event gate.

For multifoil data, each detected foil gets its own `TCutG`.

---

## Expert hand cuts

The expert-human cut file is imported only for overlay and benchmarking.

Expected file pattern:

```text
cuts/ytar_delta_<run>_<handFileID>_cut.root
```

The script attempts to load:

```text
delta_vs_ytar_cut_foil0
delta_vs_ytar_cut_foil1
delta_vs_ytar_cut_foil2
...
```

If standard names are not found, it collects all `TCutG` objects and sorts them by mean `ytar`.

Expert cuts are never used to build the auto cuts.

They are used for:

```text
visual overlay
nHand
nBoth
comparison diagnostics
```

---

## Output files

The multifoil version writes:

```text
cuts/ytar_delta_<run>_<tag>_multifoil_cut.root
plots/ridge_width_central_guard_multifoil_run<run>_<tag>.pdf
plots/ridge_width_central_guard_multifoil_run<run>_<tag>.csv
```

The ROOT cut file contains:

```text
delta_vs_ytar_cut_foil0
delta_vs_ytar_cut_foil1
delta_vs_ytar_cut_foil2
```

depending on how many ridges are detected.

It also writes diagnostic graphs:

```text
left_boundary_vs_delta_foilN
right_boundary_vs_delta_foilN
ridge_center_vs_delta_foilN
raw_ridge_center_vs_delta_foilN
left_width_vs_delta_foilN
right_width_vs_delta_foilN
```

---

## Diagnostic PDF contents

The multifoil PDF includes:

1. Central-delta `ytar` density plot with detected ridges labeled.
2. Base-cut `ytar` vs. `delta` histogram with all auto cuts and expert cuts overlaid.
3. One detail page per foil:
   - base `ytar` vs. `delta` histogram
   - auto cut
   - matching expert cut if available
   - event counts
   - central-width calibration values
4. One boundaries/center page per foil:
   - left boundary
   - right boundary
   - smoothed ridge center
   - raw ridge center

A bug was fixed by using an explicit dummy frame or explicit axis range so the boundary plots span the correct `ytar` range instead of autoscaling too tightly.

---

## Main configurable values still present

Current main values:

```cpp
const double deltaStep = 0.25;
const double deltaHalfWindow = 0.35;

const double centralDeltaAbsMax = 1.5;

const double minPeakFracOfMax = 0.20;
const double minPeakSeparation = 0.20;

const double widthFracLevel = 0.25;

const int smoothPasses1D = 2;
const int smoothPassesEnvelope = 8;

const double minPeakFracOfGlobal = 0.75;
const double guardScale = 1.25;

const int minRowEntries = 40;
```

However, the goal is to keep reducing unnecessary tunable parameters.

The most important remaining method choices are:

```text
widthFracLevel = where the ridge edge is defined
guardScale = how much wider than central reference width is allowed
minPeakFracOfMax / minPeakSeparation = what counts as a dominant ridge
```

---

## Future improvement: remove hard-coded fit range

Earlier single-ridge versions had a hard-coded `fitYtarMin/fitYtarMax`.

The multifoil version moves away from this by detecting central ridges across a broad `ytar` range and then assigning each ridge a local search window based on nearest-ridge separation.

This is better because it avoids hard-coding the foil location.

Still, future cleanup should continue removing or documenting any hidden geometry assumptions.

---

## Future step: sieve-hole selection

The next project step is to use the selected foil events as input for sieve-hole ownership.

Likely next-stage direction:

```text
For each foil-indexed event set:
  project into xsieve vs ysieve
  work delta-slice-by-delta-slice
  detect sieve-hole clusters / ellipses
  assign hole ownership
  compute quality scores
```

Important continuity from this ytar work:

```text
Each event can now be associated with a foil index based on which ytar TCutG it passes.
Foil index ordering is deterministic: increasing ytar.
```

The same philosophy should carry forward:

```text
data-driven first
geometry-informed guardrails second
expert cuts only for validation
avoid unnecessary tunable parameters
track all assumptions explicitly
```

---

## Key lesson from the ytar work

The best-performing version was not the most complicated KDE/contour version.

The robust solution was:

```text
dominant ridge detection
smooth ridge center
density-defined width
central-width guardrail
simple foil indexing
```

This is probably the right pattern to keep in mind for sieve-hole ownership too:

```text
Let the data find the structure.
Use known geometry only to prevent obviously unphysical behavior.
Do not overfit visual hand-cut aesthetics.
Judge final quality by downstream optics stability.
```
