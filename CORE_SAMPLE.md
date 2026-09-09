# Core samples for angular optics

This prototype selects the central population of each labeled sieve hole, reserves
events for independent evaluation, and limits how many cores enter the fit sample.
It reads the existing X and Y candidate trees, before GMM cleanup. It writes to new
campaign directories and does not change graphical cuts, GMM outputs, or matrices.

Start by looking at the sieve plots. The thresholds below are starting values,
not validated physics cuts.

## 1. Run on HMS 6.667

Use the `core-sample-campaign` branch, based on `hms-shms-campaign-switch`.
The earlier `core-sample` branch was based on `main` and lacks the campaign branch's
tracked setup files. Do not use that earlier branch for the ifarm test.

Run these commands from the repository directory in your usual analysis environment.
The selector needs Python, NumPy, SciPy, Matplotlib, uproot, and scikit-image. It does
not require PyROOT. Install missing dependencies once in your Python environment:

```bash
python3 -m pip install -r requirements_core.txt
```

Both candidate stages must already be complete. If they are not:

```bash
./run_make_yscol_candidates_all.sh HMS_6p667GeV
./run_make_xscol_candidates_all.sh HMS_6p667GeV
```

Those stages apply the campaign's existing pre-GMM cut vetoes. Rebuild both after
changing cuts, vetoes, delta slices, or the source reconstruction. The new selector
does not reinterpret those cuts.

Check the inputs, then produce the sample:

```bash
./run_core_sample.sh HMS_6p667GeV --check
./run_core_sample.sh HMS_6p667GeV
```

`--check` opens both candidate trees, verifies their required columns and agreement,
and reports exclusions. It writes nothing. Selection reads the campaign table
`config/rungroups_*_inputs.tsv` and these existing files:

```text
04a_candidate_trees_y/root/YscolCandidates_<rungroup>.root  (TYCand)
04b_candidate_trees_x/root/XscolCandidates_<rungroup>.root  (TXCand)
```

The default output is `HMS_6p667GeV/05c_core_sample/core/`. The final directory is
published only after the entire selection and plotting step succeeds.

## 2. Review the plots

Open the PNGs in `05c_core_sample/core/plots/`. Each foil/delta slice has six panels:

1. Unambiguous joint labels, with accepted core contours and X,Y hole IDs.
2. Cores allocated to the fit.
3. Cores in the protected random holdout.
4. Surplus cores.
5. Development events outside accepted cores, including unsupported holes.
6. Protected holdout events outside accepted cores, including unsupported holes.

Horizontal is Y sieve; vertical is X sieve. Panels share axes, binning, and the
logarithmic count scale. One-count bins remain visible. Counts outside the plotted
frame are printed in the panel title. Plot bounds come from development data;
extreme holdout outliers cannot enlarge or otherwise alter the learned regions.

Check that sparse holes survive, cores sit inside the right populations, and the
fit and held-out cores occupy similar regions. A clean-looking picture alone does
not establish that the matrix fit has adequate focal-plane or delta coverage.

The tables are small enough to inspect directly:

| File under `tsv/` | What it tells you |
|---|---|
| `counts.tsv` | Fit, holdout, holdout-core, surplus, and unsupported counts per hole. |
| `regions.tsv` | Accepted or flagged holes, peak positions, core means, and stability. |
| `allocation.tsv` | Shared pool, available cores, and fit/surplus allocation. |
| `core_<rungroup>_excluded.tsv` | Original entries excluded for missing axes, multiple labels, or invalid values. |

`no_development` means the entire slice was reserved; it has no fitted region or
plot. `low_stats`, `weak_peak`, `seed_collision`, and `unstable` are explicit failures
to establish a core. Their events remain recorded, not silently relabeled.

## 3. What the selector does

**Join the labels.** An event is identified by its rungroup and original entry in
the merged source tree. Entry numbers remain 64-bit. An X label and Y label must
agree on foil, delta slice, coordinates, and the other shared selection fields.
Multiple labels on either axis are excluded rather than resolved by guessing.
The holdout therefore represents the unambiguous, finite joint-label population,
not all triggers or all events that passed a single-axis cut.
Campaign rows must name disjoint original runs and source files; duplicate original
runs across rungroups are rejected to avoid counting the same event in both fit
and holdout. This check relies on the campaign table's provenance being accurate.

**Reserve the holdout first.** In each rungroup/foil/delta/X-hole/Y-hole stratum,
reserve 20%, rounded up, using a stable hash of the seed, rungroup, and entry.
A singleton is reserved. Density estimation, spacing estimation, peak location,
and contour tightening use only the remaining development events. The same frozen
regions then classify holdout events as core, noncore, or unsupported.

**Find locally seeded density regions.** Coordinates are scaled by hole spacing,
estimated from labeled development populations. This does not impose an ideal
centroid or a fixed hole count. A smoothed histogram uses about 30 bins per spacing
and a smoothing width of 0.12 spacing. A hole's own labeled density seeds its peak,
within 0.4 spacing of its robust initial position (the coordinate-wise median).
A watershed of the combined development density separates the seeded populations.

**Use local contrast.** For each seed, measure the peak density `p` and a local
background `b`: the median density in an annulus 0.45–0.70 spacings away. Accept the
connected part of that hole's watershed region satisfying:

```text
q = (density - b) / (p - b)
q >= core
```

The default `core` is 0.50. This gives each hole a different absolute density
threshold. There is no common event-density cutoff across the sieve plane.
The annular background is a practical first estimator; it is not a fitted
background model or an exact neighboring saddle level. The score is not a
probability of a correct label. The contour is also restricted to 0.65 spacing
from its seed and is applied only to events bearing that hole's label.

**Check support and stability.** Require at least 20 development events for a seed,
a peak above twice the local background, and at least five development core events.
Compare the event means at thresholds `core - 0.15` and `core + 0.15`, bounded to
0.05–0.95. Both comparison contours need five events. Their separation must be at
most 0.15 in spacing-scaled coordinates. The reported core center is the arithmetic
mean inside the accepted contour; the density peak is reported separately.
These observed centers select events. They do not replace the physical hole
coordinates used as fit targets.

This version shares a smoothing scale across holes, but does not yet borrow
centroids from adjacent delta slices, resolve overlapping labels, or impose a
focal-plane quality cut. It cannot establish a reliable core from arbitrarily few
events. These limits are visible in the status table.

## 4. Core quality and fit allocation are separate

Being a core makes an event eligible. It does not automatically put it in the fit.
The default allocation has three limits:

- At most 80% of a rungroup's development cores in an allocation cell can enter
  the fit. Rounding down leaves at least one surplus core whenever cores exist.
- At most 400 fit events per pool/physical-foil/delta-interval/X-hole/Y-hole cell.
- At most 200,000 fit events across the complete campaign.

First distribute the campaign budget fairly among populated cells, capped at 400.
Then distribute each cell's budget fairly among participating rungroups, giving
unused shares from sparse groups to groups with more cores. Within each group,
choose events using a second stable random order, independent of the holdout order.
This spreads fit and surplus events throughout the accepted core. It does not
always put the most central events into the fit. Small budgets can leave some
individual rungroups underrepresented; inspect `allocation.tsv`.

Surplus cores are useful validation data, but they helped establish the regions.
They are kept separate from the protected holdout. Changing the fit cap, total
budget, or pools leaves the protected holdout membership unchanged. Rebuilding
candidates, changing the seed, or changing the holdout fraction defines a new split.
Adding a new rungroup leaves existing rungroup reserves unchanged, provided their
candidate membership and identities are unchanged.

By default each rungroup has its own pool. Different angles/settings are not
automatically assumed to be interchangeable. To share a cap between known
comparable groups, use a one-time campaign configuration, for example:

```json
{
  "pools": {
    "rg03_theta12p490_foilpm8": "theta12p49",
    "rg07_theta12p495_foilpm8": "theta12p49"
  }
}
```

This example treats the two angles as comparable; make that choice deliberately.
Only identical physical `zfoil`, delta edges, and X/Y hole IDs share the cap.
Even groups in the same pool keep separate density maps and holdout reserves.

The current producers write `optics_id` to `run`, even when several original runs
were merged. We therefore cannot balance the original runs within a merged group
from these trees. That requires preserving original run/event provenance upstream.

## 5. Adjust and repeat with short commands

Repository defaults live in `core_sample.json`. Put only campaign-specific changes
in `HMS_6p667GeV/config/core_sample.json`. For a tighter trial:

```json
{
  "core": 0.65
}
```

Then run with a short tag:

```bash
./run_core_sample.sh HMS_6p667GeV tight
```

Outputs go to `05c_core_sample/tight/`. Existing tags are never overwritten. The
tag changes only the output location; it does not change random assignments.

The main controls are `core` (higher is tighter), `smooth` (in hole spacings),
`fit_cap`, `fit_max`, and `fit_fraction`. `min_events` changes the minimum evidence
for attempting a hole. Lowering it is not a guarantee of a credible sparse core.
Leave `seed` and `holdout` fixed while comparing selectors.

`spacing` normally stays `null`. If there are too few rows/columns to infer spacing,
or the inference is visibly wrong, set `[X spacing in cm, Y spacing in cm]`.
For the current HMS geometry the nominal starting values are `[2.54, 1.524]`.
SHMS campaigns can run the selector with their labels and spacing; no nine-hole
assumption is made by selection.

Each completed tag saves the resolved `config.json`, the campaign table, frozen
maps in `models/`, selector source and resolved defaults in `code/`, exact versions
of the five direct Python dependencies in `requirements.txt`, and a `manifest.json`
with input/code/output SHA-256 hashes and counts. ROOT files include every usable
joint event, so sample membership is preserved without rerunning the algorithm.

To reproduce an earlier tag after edits, restore its saved `config.json` as the
campaign's `config/core_sample.json`, use its saved dependencies in an isolated
Python environment, and run its archived selector with an absolute campaign path
and a new output tag. Keep the candidate files unchanged; compare their hashes to
the old manifest. Restoring numerical dependency versions and code matters for
exact boundary reproducibility. ROOT file container metadata need not be bytewise
identical even when all event assignments are identical.

The code archive can be invoked as:

```bash
python3 HMS_6p667GeV/05c_core_sample/core/code/core_sample.py \
  "$PWD/HMS_6p667GeV" repeat
```

During visual tuning, regard plots involving protected holdout events as validation.
If those plots or residuals guide repeated decisions, this reserve is no longer an
untouched final test set. Use fresh withheld data for the final unbiased assessment.

## 6. Build compatible HMS fit trees

After reviewing a sample, build its `TFit` trees:

```bash
./run_build_core_fit.sh HMS_6p667GeV
./run_build_core_fit.sh HMS_6p667GeV core holdout
./run_build_core_fit.sh HMS_6p667GeV core surplus
```

The default tag is `core`, and the default sample is `fit`. For the tighter trial:

```bash
./run_build_core_fit.sh HMS_6p667GeV tight
```

This requires the original replay ROOT files and HCANA. Set `PYTHON` to choose a
Python interpreter, or `HCANA=root` to use ROOT where HCANA is unnecessary.

The builder uses the existing HMS target calculations in
`make_fit_ntuple_from_gmm.C`, with an optional core-mask input appended to its
arguments. Existing calls retain their behavior. The adapter verifies that all
requested event IDs, and only those IDs, appear in the output. It checks reconstructed
sieve coordinates against the source tree and refuses a mismatched replay.
It also checks that foil positions and delta edges still match the optics metadata.

The separate outputs are under `06c_core_ntuple/<tag>/<sample>/`. Filenames retain
`Optics_<optics_id>_-1_fit_tree_gmm.root` for compatibility with the existing SVD
reader. The containing directory identifies their actual source. The `TFit` trees
also carry `sample` and `core_keep`, allowing the protected holdout's core events
to be evaluated separately. `build.json` records metadata, macro, and geometry-header hashes.

**This command builds trees; it does not launch SVD.** The existing SVD reader can
read the `fit/root` directory, but its internal count cuts and input-order limits
can further subsample it. Do not interpret the allocated fit count as the exact
number used by that solver without checking its report. Surplus here means events
not allocated to the input fit sample, not all events ultimately unused by SVD.
Never point the optimizer at `holdout/root` or `surplus/root`.

The core fit command currently accepts HMS campaigns only. It now calls the
campaign-selected geometry in the shared HMS/SHMS macro, preserving the existing
SHMS GMM workflow. Extending and validating the core fit command for SHMS is a
separate step; the selector itself already supports other hole counts.

## 7. ROOT output definitions

`root/CoreSample_<rungroup>.root` contains a real TTree named `CoreSample`, the
numeric Y-candidate fields, the joined `xscol`, and these added fields:

| Field | Meaning |
|---|---|
| `quality` | 0 unsupported; 1 noncore in an accepted hole; 2 core. |
| `core_keep` | 1 if quality is core, independent of sample assignment. |
| `core_score` | Local contrast score; NaN when unsupported or outside the map. |
| `sample` | 0 development noncore/unsupported; 1 fit; 2 protected holdout; 3 surplus cores. |

Use `(rungroup, entry)` for identity across files, not `entry` alone. No event
appears in two sample assignments. Exclusions are recorded separately in TSV.

## 8. Development checks

Run the bounded regression suite:

```bash
python3 -m unittest discover -s tests -p 'test_core_sample.py' -v
```

It exercises sparse-hole retention, a missing hole, holdout independence, repeatable
assignment, ambiguous labels, 64-bit entries, shared budgets, and ROOT output.
To create a small synthetic campaign and its plots without real replay data:

```bash
python3 tests/test_core_sample.py --demo /tmp/HMS_core_demo
HCANA=root ./run_build_core_fit.sh /tmp/HMS_core_demo
```

The fixture has eight populated holes, one missing hole, and a 24-to-1 occupancy
ratio. It checks mechanics, not performance on real HMS or SHMS data. Real
HMS 6.667 candidate trees were absent from the development workspace; the first
campaign run and visual acceptance check remain to be done on the analysis host.
