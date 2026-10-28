# Continuous-delta surface selection on ifarm

This branch constructs **pre-GMM X and Y populations from original replay**.
It fits one ordered surface family per foil/axis across continuous delta, then
writes the existing candidate-tree format. It does not run GMM, the separate
density cleanup, fit-ntuple construction, or SVD automatically.

Branch: `codex/joint-surface-ifarm`, based on `505e656`.
The `codex/provisional-geometric-clustering` branch remains the fallback.
Use a new campaign directory; never mix its numbered outputs with an earlier
campaign. The original selectors and downstream programs remain available.

The illustrated [method guide](docs/JOINT_SURFACE_METHOD_GUIDE.pdf) explains the
prior information, mathematics, boundaries, tunable parameters and pilot results.
The [viewer guide](tools/membership_viewer/README.md) covers native display,
offline Firefox, sampling and the saved-data format.

### Farm-test correction on `joint-surface-test`

Population testing now uses **all eligible geometry-training events** in each
reporting interval by default (`occupancy_cap: null`). The previous 2,000-event
cap discarded enough evidence to miss the sparse outer Y ridge in rg02 over
delta [-8,0)%. An explicit positive `occupancy_cap` remains available for an
intentional capped comparison; omit that setting or use `null` for this rerun.

Use top-level `coverage_min_events: 5` for the planned count-sensitivity rerun.
The counting and delta-gap rules themselves are unchanged. Run in a fresh
campaign so the old frozen outputs remain available for comparison. The PDF
version 1 describes the original 2,000-event default; this correction supersedes
that one population-test setting. A full refit can also choose a different
family because population densities enter model-choice scoring.

## 1. Get the branch and check the environment

From an existing clean checkout of this repository:

```bash
git fetch origin
git switch -c codex/joint-surface-ifarm --track origin/codex/joint-surface-ifarm
```

If the branch already exists locally, use `git switch codex/joint-surface-ifarm`
and `git pull --ff-only`. A separate checkout is appropriate if another branch
is in active use. Run the remaining commands from this checkout's root.

Use the usual HCANA/ROOT environment. The numerical process needs Python 3.9+
with NumPy and SciPy. ROOT I/O needs PyROOT and NumPy. The native viewer needs
Matplotlib and Tk; offline HTML export optionally needs the Python `plotly`
package. GMM later needs its normal combined PyROOT/scikit-learn environment.

```bash
command -v hcana
python3 -c 'import numpy, scipy; print("Numerical dependencies OK")'
python3 -c 'import ROOT, numpy; print("ROOT", ROOT.gROOT.GetVersion())'
python3 -c 'import matplotlib, tkinter; print("Native viewer imports OK")'
```

`requirements-joint-surface.txt` lists the numerical/viewer Python packages.
It deliberately does not install ROOT, HCANA or a different system Python.
If the numerical Python and PyROOT Python differ, set:

```bash
export HALLC_ROOT_PYTHON=/path/to/the/PyROOT/python3
```

The driver uses that executable only for export and candidate ROOT I/O.
The executable invoking the driver runs the numerical fits. Alternatively pass
`--root-python /path/to/the/PyROOT/python3` on each stage command.
`--root-command hcana` selects the stage-01 executable; `root` can be used in
an environment where the existing foil macro runs under ROOT directly.
For an interactive single-process farm test, setting `OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1` avoids excess BLAS threads during
the many small coefficient solves.

## 2. Create a fresh same-arm campaign

The source campaign must already contain its rungroup input TSV and
`config/oldfit.dat`. Replay files named by the input table must be original
replay files with arm branches, not candidate or TFit files. The metadata
`optics_id` records must be present in `DATfiles/list_of_optics_run.dat`.

```bash
python3 setup_joint_surface_campaign.py \
  /path/to/stable/HMS_6p667GeV HMS_6p667GeV_joint_trial

C=HMS_6p667GeV_joint_trial
python3 run_joint_surface_campaign.py "$C" --plan
```

For SHMS, substitute a stable `SHMS_...` source and a new `SHMS_...` destination.
The folder prefix selects the arm automatically:

| Profile | Replay prefix / Cherenkov | Acceptance | Sieve |
|---|---|---|---|
| HMS | `H.*` / `H.cer.npeSum` | `-10 < delta < 10` | 9 rows, 9 columns |
| SHMS | `P.*` / `P.ngcer.npeSum` | `-10 < delta < 22` | 11 rows, 11 columns; centered sieve only |

Reporting delta edges and foil positions come from run metadata, not hardcoded
HMS intervals. A different arm does not reuse an HMS fitted model.

Setup copies the input table and seed matrix, plus `core_sample.json`,
`core_balance.json` and `sieve_mask.json` from the source config when present.
For `HMS_6p667GeV`, these carry the explicit `min_events: 10` core threshold,
equal15 policy and blocked holes `(2,3)` and `(5,5)`. Relative replay paths are
resolved against the source checkout. Old component overrides and cut vetoes
are archived as reference material, not activated. The new campaign marker is
`config/joint_surface_campaign.json`; settings are in
`config/joint_surface.json`. Existing destination directories are refused.

Run metadata is validated for every campaign optics ID and copied to
`config/optics.dat`. Setup uses the source campaign's `config/optics.dat` when
present, otherwise repository `DATfiles/list_of_optics_run.dat`; pass
`--metadata /path/to/list_of_optics_run.dat` to setup to choose another source.
The campaign marker records the copied files' hashes. The surface driver uses
this snapshot by default and refuses a different metadata override.

The plan reads configuration and metadata but does not access replay files or
write stage outputs. Restrict an initial test to one familiar rungroup:

```bash
RG=rg01_theta12p490_foil0
python3 run_joint_surface_campaign.py "$C" --plan --rungroup "$RG"
```

Omit `--rungroup` to process all rows; repeat it to choose several groups.

## 3. Execute the five explicit stages

```bash
python3 run_joint_surface_campaign.py "$C" ridge --rungroup "$RG"
python3 run_joint_surface_campaign.py "$C" export --rungroup "$RG"
python3 run_joint_surface_campaign.py "$C" fit --rungroup "$RG"
python3 run_joint_surface_campaign.py "$C" label --rungroup "$RG"
python3 run_joint_surface_campaign.py "$C" candidates --rungroup "$RG"
```

Add `--dry-run` to print one stage without executing it. For a new campaign,
choose custom metadata during setup; every stage then uses its saved snapshot.
The stage-level `--metadata` option remains available for older campaigns that
have no snapshot, or for a byte-identical copy of a new campaign's metadata.

1. **ridge:** runs the existing foil-selection macro and writes the conventional
   `01_ytar_cuts/cuts/ytar_ridge_cut_TAG.root`. Inspect the foil polygons before
   admitting events for the surface study.
2. **export:** reads replay `T`, or `Tout` if `T` is absent. It applies the existing
   strict electron PID, arm acceptance, metadata delta intervals and foil
   polygons. If foil polygons overlap, the highest matching foil index wins,
   as in the candidate makers. No old X/Y band cut, GMM mask, density-core
   selection or historical row/column label is applied. Original global
   replay `entry` is retained.
3. **fit:** fits X and Y from their measured focal-plane position/slope and
   event-level delta. Default event roles are 80% for frame/coefficients/
   occupancy/envelopes and independent 20% for initializer/family choice.
   All exported events then receive frozen predictions. One primary fit is
   used; the exploratory two-fit audit is not required to make candidates.
4. **label:** verifies geometry hashes, attaches one physical orientation and
   integer origin per foil/axis, and writes the viewer manifest. By default,
   reconstructed sieve medians align the already frozen bands with the known
   hardware grid. Manual global maps avoid reconstructed sieve input here.
5. **candidates:** verifies frozen geometry and memberships, then writes TXCand
   and TYCand with the conventional branches, types, original entries and
   summary TSVs. X and Y populations remain independent until the downstream
   join. No 2D TCutG approximation is substituted for the continuous selection.

All delta intervals include their lower boundary and exclude their upper
boundary. Strict overall arm acceptance still excludes the acceptance limits.
Nonfinite values and out-of-range events have explicit rejection reasons.

## 4. Inspect the saved outputs

| Location | Contents |
|---|---|
| `00_joint_surface_work/TAG/preband/` | Original replay event export, admission provenance and counts |
| `00_joint_surface_work/TAG/foilN/{x,y}_model.json` | Fixed frame, selected family/coefficients, periodic seed, occupancy, training-envelope tables, full settings |
| `.../foilN/{x,y}_fit_metadata.json` | Training/model-choice indices, candidate scores, initialization decisions |
| `.../foilN/{x,y}_prediction.npz` | All admitted events' latent assignments, proposed index, coordinates, flags and reasons |
| `.../geometry_frozen.json` | Hashes fixed before physical naming |
| `.../physical_maps.json` | Global names, reconstructed median alignment or manual-map provenance |
| `.../foilN/membership.npz` | Measured coordinates, exact entries, physical/latent assignments, acceptance/reasons and event-role masks |
| `.../labels_frozen.json` | Physical-map and membership hashes checked before candidate publication |
| `00_joint_surface_work/viewer_manifest.json` | Reusable index of all labeled foil cases |
| `04a_candidate_trees_y/{root,tsv}/` | `YscolCandidates_TAG.root`, TYCand, conventional summary TSV |
| `04b_candidate_trees_x/{root,tsv}/` | `XscolCandidates_TAG.root`, TXCand, conventional summary TSV |

The candidate schema retains 64-bit original entries. The unchanged downstream
GMM mask format uses signed 32-bit entries, so this adapter explicitly refuses
original entries above 2,147,483,647 rather than allowing truncation.

Source replay and prior campaigns are opened read-only. Written stage outputs
are protected from replacement. An interrupted or failed attempt remains for
inspection; use a fresh campaign rather than deleting pieces and mixing results.

## 5. Rotate the measured clouds without a browser

Connect with X forwarding and a running local X server (for example XQuartz on
macOS). From a graphical session on ifarm:

```bash
python3 -m joint_surface.viewer "$C" --desktop --max-points 4000
```

Select the foil/run, saved configuration, delta interval, coloring and event
filter. Drag to rotate. Click a point for its exact entry and both memberships.
The X-row or Y-column shared-color modes follow one population through both
planes; legend hiding is linked in those modes. Rejection reasons have distinct
cross colors. Event detail includes saved training-role information.

If the window cannot open, check `echo "$DISPLAY"` and
`python3 -c 'import tkinter'`. Tk is part of the Python/system installation,
not supplied by pip installing Matplotlib. Run the viewer from a graphical
login session rather than an ordinary headless batch job. `--max-points 0`
plots all finite-coordinate events but can be slow over X forwarding.

For offline Firefox instead:

```bash
# Optional one-time addition to the analysis environment:
python3 -m pip install plotly
python3 -m joint_surface.viewer "$C" --html "$C/membership.html"
firefox "$C/membership.html"
```

The exported file embeds the JavaScript and plotted data; it needs neither a
web server nor internet access. A local bundle can replace the optional Python
package via `--plotly-js /path/to/plotly.min.js`. Firefox needs a usable WebGL
context. If forwarded graphics cannot provide one, use the native viewer or
copy the HTML to a local machine with working graphics. Native/Firefox runtime
interaction on ifarm is not yet verified; that is part of the first farm test.

The default 12,000-event display is sampled, with deliberate representation of
rejection categories. It is not an unbiased density estimate and does not
necessarily show every rejection. Full-population counts are shown separately.
Changing a viewer filter never refits a surface or changes a candidate mask.

## 6. Review and tune deliberately

The model uses the full delta range to predict each eligible ridge's q center
at the event's own r and delta. Population evidence remains measured separately
in reporting intervals. Overlap in a 2D projection can be disentangled using
continuous delta and the other FP coordinate; coincident populations at the
same measured coordinates cannot be uniquely identified.

Edit `config/joint_surface.json` before running fit. Example:

```json
{
  "schema_version": 1,
  "selection": {"cal_cut": 0.65, "cer_cut": 6.0, "tree": "auto", "hash_replay": false},
  "coverage_min_events": 10,
  "fit": {
    "common": {"families": ["bilinear", "curved", "static", "shift", "spacing", "tilt", "divergence"]},
    "x": {"theta_deg": 135.0},
    "y": {"theta_deg": 135.0}
  },
  "global_maps": {}
}
```

The engine's complete defaults are in `joint_surface/engine.py` and are saved
with every model. Axis keys override `fit.common`; unknown keys are errors.

| Observation | First relevant settings |
|---|---|
| Wrong periodic spacing/count, doubled or missed ridge | `min_q_frequency`, `max_q_frequency`, `phase_starts`, `phase_iterations`; inspect raw density and initialization alternatives |
| Unnecessary curvature or delta variation | `families`, `high_order_penalty` (0.00002), `extra_penalty` (0.002), `model_choice_tolerance` (0.01 nat) |
| Geometry-fit center sample too narrow/broad | `seed_radius` (0.24), `refine_radius` (0.30); these are training radii, not final event cuts |
| A predicted position lacks population evidence | `occupancy_likelihood_gain` (10), `occupancy_min_effective_events` (8), `occupancy_min_coherent_intervals` (2), `occupancy_min_events` (100) |
| Ridge events rejected by training-envelope counts | Top-level `coverage_min_events` (10; 5 is the sensitivity setting), plus envelope diagnostics |
| Tail or delta-gap rejections | `coverage_quantiles` ([0.005,0.995]), `coverage_padding` (0.04 standardized r), `coverage_bin_width` (1%), `maximum_bridge` (2%) |
| Whole physical family shifted/reversed | `global_maps`, after inspecting the frozen geometry |

The label stage allows only `global_maps` changes after fit. To specify a manual
map before label:

```json
"global_maps": {
  "rg01_theta12p490_foil0/foil0": {
    "x": {"orientation": 1, "origin": 0},
    "y": {"orientation": 1, "origin": 1}
  }
}
```

This means `physical_index = orientation * latent_index + origin` for every
event of that foil/axis across delta. It is not an event-specific relabel.
An incorrect periodicity cannot be fixed with a naming offset. Old
`pre_gmm_veto.tsv` entries name TCutGs that do not exist in the direct surface
path; active legacy cut-veto entries are refused rather than silently ignored.

For the first count-sensitivity campaign, use another fresh destination,
change only top-level `coverage_min_events` from 10 to 5, then run the same stages
with the same input, numerical settings and seed. The driver recomputes the
fit for that fresh campaign, but the count is not part of the coefficient
objective. The engine's `predict_axis(..., coverage_min_events=5)` API can
reevaluate a saved model without fitting; the protected stage driver does not
replace published primary predictions in place. Compare coefficients as well
as counts when verifying that only the envelope rule changed.

The primary classifier has no extra posterior/core cut. Large-q tails can be
assigned if they pass its geometric and training-envelope rules. Do not confuse
assignment with electron purity or foil-vertex provenance. The separate density
cleanup branch and sparse-tail origin study are not folded into this branch.

## 7. Validation and the next handoff

Completed locally for this branch:

- Seven targeted numerical tests: ordering, independent prediction, missing-band
  control, SHMS delta through +22%, nonfinite/edge handling, and count-only changes.
- Frozen-model port check: all 46,516 central X pilot evaluation predictions
  reproduced exactly, including q/r/u, assignments and acceptance.
- Campaign unit checks for fresh setup, planning, naming and overwrite protection.
- Reproducible 14,000-event numerical adapter smoke through fit, naming and
  candidate payload preparation: both axes accepted 13,916 events, all matching
  the generated band identities. Altering a pre-band coordinate with entries
  unchanged correctly blocked publication through the frozen input hash.
- Synthetic replay ROOT I/O for HMS/SHMS and T/Tout: strict PID/acceptance,
  half-open delta, foil overlap precedence, original entries, candidate branch
  types/counts, overflow refusal, overwrite refusal and unchanged input hashes.
- Eight actual unchanged GMM consumer fits (X/Y, both arms, both tree names)
  accepted the new candidate format and retained exact fixture entries with
  keep fraction 1. This checks software compatibility, not the default 90%
  cleanup's physics behavior.
- Seven viewer contract/render tests and a real offline HTML export with syntax
  checks. Desktop/X forwarding and Firefox interaction still need ifarm testing.
- The PDF was rendered and visually inspected.

Reproduce the focused checks when useful:

```bash
OPENBLAS_NUM_THREADS=1 python3 tests/test_joint_surface_engine.py
python3 tests/test_joint_surface_campaign.py
OPENBLAS_NUM_THREADS=1 python3 tests/joint_surface_numeric_smoke.py
MPLBACKEND=Agg python3 -m unittest discover -s tests -p test_joint_surface_viewer.py

# PyROOT + NumPy interpreter:
python3 tests/joint_surface_root_smoke.py
# Add --gmm-python /path/to/combined/python to test the unchanged GMM consumer.
# If separate, invoke the preceding command with "$HALLC_ROOT_PYTHON" instead.
```

These checks do not establish full-replay band quality, electron purity,
another nominal-momentum setting, SHMS physics performance or SVD improvement.
The local pilot used previously selected candidate samples. The ifarm replay
run is the first test before old band selection; no prior-band fallback is used
if discovery fails. A failure leaves the stage diagnostics and stops the driver.

For the HMS 6.667 training sample, feed the completed sheet-labeled candidate
trees directly into the existing density-core selector and frozen allocator:

```bash
./run_core_sample.sh "$C" min10 --check
./run_core_sample.sh "$C" min10
./run_core_balance.sh "$C" min10 equal15
```

`min10` is an output tag; the threshold comes from `config/core_sample.json`.
The surface `coverage_min_events` setting is independent of this density-core
threshold. Setup inherits the existing seed and holdout settings. Core selection
reserves events before estimating density; equal15 reallocates only development
cores and preserves the reserved IDs. The core tag saves and hashes
`metadata/optics.dat`, so allocation needs no earlier SVD diagnostics. Existing
tags are never replaced. Additional GMM cleanup is not required for this handoff.

Stop after allocation until the chosen `xptarT` correction is consistently
implemented and verified before ntuple export. These commands do not export
ntuples or run SVD, elastic net or beam search.

PDF regeneration uses distributed figures and requires ReportLab and Pillow:
`python3 docs/source/build_joint_surface_guide.py`. The optional scientific
figure-generation source requires the original frozen pilot data at an explicit
`--project` path; it is not needed to run the farm software.
