# Provisional geometric X/Y selection

Branch: `codex/provisional-geometric-clustering`, based on `ff2247b`.
The original angle-scan selectors and split-X runner remain available. Use a
fresh campaign to compare methods; do not mix their numbered stage outputs.

This alternative replaces step 2 and supplies metadata-driven step-3 runners.
The foil cuts, physical sieve relabeling algorithms, candidate trees, GMM
cleanup, fit ntuples and SVD remain the existing workflow. The six replay-reading
macros now accept `T` or `Tout`, preferring `T` when both exist.

## Create and run a campaign

Run from this repository's root in the normal HCANA/ROOT environment. Python
needs NumPy, SciPy and Matplotlib for the selectors; the unchanged GMM scripts
also need scikit-learn and PyROOT in the same Python environment. The existing
software requirements in README still apply. No research-directory imports or
local development paths are required.

Copy configuration from an existing campaign into a **new** same-arm directory:

```bash
python3 setup_geometric_campaign.py \
  /path/to/stable/repository/HMS_6p667GeV HMS_6p667GeV_geometry_trial

CAMPAIGN=HMS_6p667GeV_geometry_trial
python3 run_geometric_campaign.py "$CAMPAIGN" --plan --tag geometry_trial
```

Setup copies the rungroup input table and `config/oldfit.dat`; it copies no
numbered outputs. Historical override/veto files are archived under
`config/previous_selection_reference/`, since their component IDs may differ.
The new campaign's marker prevents accidentally targeting a stable campaign.
As before, each `optics_id` must exist in `DATfiles/list_of_optics_run.dat`, and
the rungroup replay paths must be accessible. The source campaign is read-only.

Execute stages in order, inspecting the foil cuts, component plots and physical
labels at their usual places:

```bash
python3 run_geometric_campaign.py "$CAMPAIGN" ridge
python3 run_geometric_campaign.py "$CAMPAIGN" selection
python3 run_geometric_campaign.py "$CAMPAIGN" relabel-y
python3 run_geometric_campaign.py "$CAMPAIGN" relabel-x
python3 run_geometric_campaign.py "$CAMPAIGN" candidates-y
python3 run_geometric_campaign.py "$CAMPAIGN" candidates-x
python3 run_geometric_campaign.py "$CAMPAIGN" gmm-y
python3 run_geometric_campaign.py "$CAMPAIGN" gmm-x
python3 run_geometric_campaign.py "$CAMPAIGN" fit --tag geometry_trial
```

`--dry-run` prints one stage; `--plan` prints the complete sequence without farm
access or output writes. `--python /path/to/python3` selects the Python runtime
for all stages. `--root-python /path/to/PyROOT/python3` can provide separate I/O
for the new selectors; GMM still needs its combined PyROOT/scikit-learn runtime.
The driver refuses an existing output stage and forces `OVERWRITE=0`.

For SHMS, use `SHMS_<campaign>` source/destination names. The shared profile
selects `P.*`, `P.ngcer`, 11 sieve columns/rows and acceptance through 22%; HMS
uses `H.*`, `H.cer`, 9 columns/rows and acceptance through 10%. Foil positions
and delta edges still come from run metadata. Centered SHMS sieve is required.
Real SHMS geometric selection remains to be validated on replay data.

## Selection and output contracts

**X:** standardize the full distribution, apply one rotation (default 135°),
discover a smooth ordered polynomial band family from periodic density, and
test occupancy separately. There is no low/high X split. For the first delta
slice, learn geometry from the first two adjacent delta intervals combined;
test occupancy and export cuts only on the original target slice. With the HMS
edges used in the study this is geometry on `[-10,-5)`, selection on `[-10,-8)`.
The remaining slices use their own events. Final delta binning is unchanged.
Dividers are halfway between accepted bands in the unfolded coordinate; they
are not independently measured valley minima. An absent slot need not be kept.

**Y:** strong peaks establish spacing, locally supported weak peaks become
candidates, and regularized two-dimensional Gaussian components refine the
partition. The automatic rotation scan remains; an explicit angle can replace
it. The final boundaries are equal-density decision curves. These are not a
fit of Y divergence. Disconnected regions and holes are retained in the cut
conversion. Nine/eleven are possible physical slots, not a required count of
detected bands.

Both fits use **pre-band replay events**, after the existing PID, foil and delta
selections. No stored X/Y labels enter discovery. Original replay `entry` is
preserved for downstream masks. Baseline PID and acceptance are strict,
delta intervals are half-open, and overlapping foil cuts assign to the highest
matching foil index, consistent with candidate construction.

| Stage | Conventional output consumed by the next stage |
|---|---|
| 01 | `01_ytar_cuts/cuts/ytar_ridge_cut_TAG.root` |
| 02a | `02a_angle_scan_y/cuts/YpFpYFp_TAG_auto_band_cut.root` |
| 02b | `02b_angle_scan_x/cuts/XpFpXFp_TAG_auto_band_cut.root` |
| 03a | `03a_relabel_y/cuts/YpFpYFp_TAG_auto_band_cut_yscolRelabeled.root` |
| 03b | `03b_relabel_x/cuts/XpFpXFp_TAG_auto_band_cut_xscolRelabeled.root` |
| 04a/b | Existing `TYCand` / `TXCand`, branch names and summary TSVs |
| 05a/b | Existing per-foil/delta `GMMClean` event masks |
| 06a/b | Existing `TFit` and `nps_hms_newfit_TAG.dat` / `nps_shms_newfit_TAG.dat` |

X provisional cuts use
`hXpFpXFp_cut_candidate_XX_zone_global_nfoil_F_ndel_D`; Y uses
`hYpFpYFp_cut_yscol_ID_nfoil_F_ndel_D`. These are **component IDs**, not final
physical labels. Y disconnected pieces receive separate provisional IDs.
Step 3 can assign several pieces to the same physical `xscol`/`yscol` and the
candidate makers OR those pieces. The new X wrapper enumerates metadata rather
than extracting a split from diagnostic filenames.

Additional artifacts use unique attempt directories:

* `00_geometric_work/TAG/ATTEMPT/`: pre-band NPZ, provenance, per-slice model
  JSON/assignments, ROOT payload and membership report, logs or unresolved notes.
* `02a_angle_scan_y/{plots,root,tsv}/ATTEMPT/` and the corresponding `02b` paths:
  diagnostic PDF, overview PNG, ROOT histograms/canvas/assignment tree, and a
  per-rungroup geometric summary TSV. Basenames retain axis, tag, foil and
  `delta_LOW_to_HIGH` conventions; `geometric` distinguishes the method.

Diagnostics include raw focal-plane density and cuts, standardized rotated
density, actual 1D projections, component/sieve views and per-component pages.
X also shows nonlinear unfolding and the pooled geometry-window histogram.
Its display rotation makes the tracks nearly vertical without changing the
selection. The fitted and display angles are both printed. Unfolding depends
on the fitted model and is not an independent validation of it.

ROOT stores the arm, model, component mapping and membership QA. Every target
event is checked against the actual `TCutG`s before publishing canonical cuts:
no multiple polygon matches, and at most 0.2% numerical-classifier disagreement.
That tolerance measures contour approximation, **not physical accuracy**.
Unresolved discovery or a failed cut check leaves diagnostic work and no
canonical cut file for that axis/rungroup. The runners fail explicitly rather
than silently advancing an unresolved slice.

## Manual correction and review

Optional `config/geometric_overrides.json` addresses individual slices:

```json
{
  "slices": {
    "rg01_theta12p490_foil0/foil0/delta0": {
      "x": {"angle": 135.0, "learning_delta": [-10, -5]}
    }
  }
}
```

The angle is the **standardized-coordinate rotation**, not the spectrometer
angle. Supported per-axis keys:

| Key | Meaning |
|---|---|
| `angle` | Replace automatic/default rotation before fitting |
| `frozen_model` | Use an existing per-slice model JSON without refitting |
| `polygons_file` | Use reviewed raw-coordinate polygons from JSON; takes precedence over discovery |
| X `learning_delta` | Two limits for the geometry-training window; target delta slice is unchanged |
| X `accepted_indices` | Explicit reviewed occupancy hypotheses |
| X `boundaries_u` | Increasing internal dividers, one fewer than accepted bands |

Paths are relative to the override file's directory unless absolute. A frozen
model already fixes its rotation/window and cannot be combined with those two
overrides. Saved model coefficients can be edited in a separate copy and frozen
for geometric adjustments. A standalone `polygons_file` bypasses discovery
entirely, including initialization failures and samples below the normal
100-event fitting minimum. Its polygons must have distinct X component IDs and
non-overlapping membership. Physical labels are still assigned in step 3.

Preview any slice without changing canonical cuts:

```bash
python3 run_geometric_selection.py "$CAMPAIGN" --axis x \
  --rungroup rg01_theta12p490_foil0 --foil 0 --ndel 0 --preview
```

The printed work directory contains `preview_cuts.root` and the exact model.
To reuse cuts edited with ROOT, save a separate reviewed ROOT file and extract
the selected slice (preserve the provisional cut names):

```bash
python3 geometric_selection/root_io.py --extract-cuts /path/to/reviewed.root \
  --axis x --foil 0 --ndel 0 --output /path/to/manual_x_f0_d0.json
```

Reference that JSON as `polygons_file`, then preview it. The same operation works
for Y. To change an already published selection, make a fresh campaign; the
selector never replaces its canonical cut files. Failed attempts can be
retried directly with the selector, which uses a new attempt directory. The
stage driver is deliberately stricter and refuses any existing stage files.

Physical-label corrections use the existing text formats:

```text
# X: optics_id foil ndel zone oldIndex newXscol action comment
666701 0 0 global 2 4 force reviewed_component
# Y: optics_id foil ndel oldYscol newYscol action comment
666701 0 0 2 4 force reviewed_component
```

Put X and Y rows in separate files. Assign multiple components the same final
column to merge them; use `-1 reject` to discard one. Pass the relevant file:

```bash
python3 run_geometric_campaign.py "$CAMPAIGN" relabel-x \
  --label-overrides /path/to/new_x_labels.tsv
```

Use `relabel-y` for Y. The exact component mapping is in the ROOT membership
report; do not transfer historical component IDs without checking. The normal
`config/pre_gmm_veto.tsv` remains available after relabeling.

## Validation and remaining campaign test

Run the targeted numerical tests and the synthetic ROOT integration:

```bash
python3 tests/test_geometric_selection.py
python3 tests/geometric_root_smoke.py --gmm-python "$(command -v python3)"
```

The smoke also supports `--numerical-python` when numerical fitting uses a
separate interpreter. Supply `--gmm-python` to test actual unchanged GMM cleanup;
without it, the test explicitly substitutes software mask fixtures.

Software validation covers both arms, `T`/`Tout`, PID/foil/delta edge behavior,
unsplit relabeling, candidate trees, unchanged GMM, fit ntuples and a **toy** SVD.
Stage-1 ridge construction is checked on a separate controlled distribution;
the main downstream fixture uses fixed foil cuts for deterministic counts.
Additional controls cover missing X slots, 11-slot geometry, Gaussian holes,
manual model/boundary/polygon precedence and overwrite refusal. The packaged
engines reproduce the original HMS prototype numerics on a saved central
slice; its 67,428 Y assignments also match exported ROOT cuts exactly.

Verified locally on 2026-10-26: nine numerical/override tests passed under both
Matplotlib 3.4.3 and 3.11.2. The full synthetic flow retained 4,950 HMS events
through `TFit` and SVD, and 6,050 SHMS events through `TFit`; 5,000 entered SHMS
SVD under the existing per-column allocation cap. Actual unchanged GMM cleanup
was used with `--keep-frac 1 --n-init 1` for deterministic event accounting;
this does not test the default 90% cleanup's physics performance. All four
exported X/Y cut files had zero assignment discrepancies,
multiple matches or unassigned fixture events. A separate two-foil/two-delta
test preserved all eight components per axis, including an early manual merge,
through cumulative relabeling on both arms. The existing spectrometer regression
also passed, including HMS `TFit` equality with the baseline implementation.

The cumulative relabel test can be rerun with:

```bash
python3 tests/geometric_relabel_cumulative_smoke.py
```

This establishes software compatibility. It does not establish calibration
quality, universal first-try success or SHMS physics performance. The next
scientific validation is a new full-replay campaign, inspection of its physical
labels and GMM selections, followed by SVD/replay comparison with the stable
matrix. The farm replay paths are not mounted in the local development session.
