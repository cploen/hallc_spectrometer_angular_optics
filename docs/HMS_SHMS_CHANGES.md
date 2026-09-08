# Campaign-selected HMS/SHMS implementation

Local branch: `hms-shms-campaign-switch`, based on `057b908956286099c49a8d453856b0a6feb846c7` of `prelim-6p667-data`.

This is a software-tested implementation on the `hms-shms-campaign-switch` review branch. Real SHMS replay validation is pending; no calibration matrix is being recommended for replay. Publication of this branch was requested before ifarm testing. No merge or direct ifarm change is part of this work.

## Selection and independent run properties

The existing campaign argument is sufficient: a whole path component named `HMS_…` selects HMS; `SHMS_…` selects SHMS. Paths to files within that campaign also work. Unknown or conflicting names fail explicitly. The code does not infer the spectrometer from the target, replay filename, or available branches. There is no added spectrometer option or prompt.

**The foil target is independent of the spectrometer.** The −10, 0, +10 cm three-foil target is available for both arms; NPS simply used other targets. No arm profile contains foil positions or a foil count. These come from each run's existing `DATfiles/list_of_optics_run.dat` entry. Momentum slice boundaries are also run metadata; the nominal acceptance limits belong to the arm profile. The campaign's existing rungroup TSV supplies the replay path and run grouping.

Centered SHMS sieve is the user-approved assumption. Geometry stages reject an explicit SHMS `SieveFlag` other than 1. A shifted-sieve position option is deferred. The existing metadata field can support that later without adding an arm argument.

| Profile setting | HMS | SHMS centered |
|---|---|---|
| Replay branch prefix | `H.` | `P.` |
| Cherenkov branch | `H.cer.npeSum` | `P.ngcer.npeSum` |
| Sieve rows × columns | 9 × 9 | 11 × 11 |
| Sieve distance | 168 cm | 253 cm |
| x positions | `(i−4) × 2.54 cm` | `(i−5) × 2.5 cm` |
| y positions | `(j−4) × 0.6 × 2.54 cm` | `(j−5) × 1.64 cm` |
| Nominal momentum window | −10% < δ < 10% | −10% < δ < 22% |
| Mispointing | Existing angle-dependent equations | xMP = −0.126 cm, yMP = −0.05 cm |
| HB correction to truth projection | None | `−0.0398 δ + 0.000398 δ² cm` |

δ is expressed in percentage points. The SHMS survey constants are fixed defaults, not something requested for each campaign. Hardware acceptance and narrower analysis slices are distinct. Existing PID thresholds are retained (ridge Cherenkov >2; most subsequent stages Cherenkov >6 and calorimeter >0.65); the switch does not claim those cuts are validated for Gema's data.

The SHMS truth geometry changes the foil-z sign, beam/vertex treatment, offsets and HB projection as well as the sieve distance. See [the reference review and equations](HMS_SHMS_REVIEW.md#shms-truth-equations-to-review-with-holly-and-mark). `P.extcor.ysieve` is read without adding another HB correction. Its actual replay convention still needs checking with Gema's configuration.

## File-by-file changes

Each changed executable also has a local comment or import identifying its use of the shared configuration.

| File(s) | Change |
|---|---|
| `.gitignore` | Exposes only the new Gema campaign input metadata and README; generated outputs remain ignored. |
| `DATfiles/list_of_optics_run.dat` | Appends runs 3283–3286 with reported +8.915° angle and approved −10/0/+10 cm, centered-sieve assumptions. Existing rows remain unchanged. |
| `SHMS_Gema_3283_3286/README.md` | Records reported settings, approved assumptions, file-listing evidence, initial delta bins and pending seed/real-data checks. |
| `SHMS_Gema_3283_3286/config/rungroups_gema_inputs.tsv` | Four initial single-run groups using the confirmed replay paths. |
| `SHMS_Gema_3283_3286/config/ztar_runlist.txt` | Three-foil run list with exact replay paths for standalone ztar checks. |
| `spectrometer_profiles.def` | Single constants table used by C++ and Python; no foil target settings. |
| `spectrometer_config.h` | Campaign resolver, centered hole coordinates, acceptance, survey/mispointing and truth equations; run metadata reader independent of arm. |
| `spectrometer_config.py` | Same profile table, campaign routing, run metadata, interval/tag matching and GMM defaults. CLI is used internally by shell runners. |
| `spectrometer_config.sh` | Shared shell adapter; validates the existing campaign name and exposes its profile. |
| `spectrometer_root.h` | ROOT branch checks, profile logging/provenance, metadata access and explicit centered-sieve guard. |
| `ytar_ridge_cut.C` | Selected arm branches and momentum window. SHMS peak ordering maps back to metadata foil IDs using the negative foil-z sine term. Peak count must match the target metadata. HMS expert-cut lookup remains HMS-specific. |
| `run_ytar_ridge_all_rungroups.sh` | Campaign-name validation; explicit campaign replay paths retained. |
| `assign_yfp_ypfp_angleScanBands.C` | Selected branches; default physical band limit from profile; delta slices from metadata; centered guard. |
| `assign_xfp_xpfp_angleScanBands_split.C` | Same routing and interval work for X, including the physical X band limit. |
| `run_xfp_xpfp_dynamicSplit_batch.C` | Selected branches, physical band limit and metadata interval loop; existing split method retained. |
| `run_y_angle_scan.sh` | Campaign validation, nominal arm window and physical band count passed to the batch macro. |
| `run_y_multifoils_fixedtheta.py` | Target/reference slices come from metadata. Rejects mismatched boundary arrays before reusing stored angles by interval index. Foil count remains run-specific. |
| `run_x_multifoils_fixedtheta.py` | Same metadata checks, plus selected branches in the replay-based split estimate. Existing split-estimate PID thresholds retained. |
| `make_angle_rerun_commands.py` | Campaign-based band counts and metadata interval lookup, including fractional boundaries; generated command context retained. |
| `relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch.C` | Selected branches, centered-sieve check, physical guides, nearest-column assignment, dimensions and plot extent. |
| `relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch.C` | Same X/Y profile geometry and branch changes; existing shell-generated campaign context retained. |
| `run_relabel_y_all.sh`, `run_relabel_x_all.sh` | Validate campaign identity before running existing chains. |
| `run_relabel_xfp_from_plot_splits_v2.sh` | Converts plot interval tags using the run's actual boundaries instead of five fixed cases. |
| `make_yscol_candidate_tree.C` | Selected branches, profile-sized arrays/loops, centered guard, nominal arm momentum window and output profile marker. |
| `make_xscol_candidate_tree.C` | Same changes; default delta limits now resolve from profile, while explicit existing overrides remain available. |
| `run_make_yscol_candidates_all.sh`, `run_make_xscol_candidates_all.sh` | Validate campaign and retain explicit replay inputs. |
| `gmm_cleanup_yscol_candidates_v2.py` | Selected raw-branch fallbacks; x-component limit, physical guides/extent and output profile marker. Generic candidate schema and algorithm retained. |
| `gmm_cleanup_xscol_candidates_v2.py` | Selected raw-branch fallbacks; y-component limit, physical guides/extent and output profile marker. |
| `run_gmm_cleanup_y_all.sh`, `run_gmm_cleanup_x_all.sh` | Default component counts from profile; existing explicit tuning overrides retained. |
| `make_fit_ntuple_from_gmm.C` | Selected replay branches; mask filenames use actual interval boundaries, including slices beyond index 4. Profile truth equations, physical column bounds, momentum limits and QA extents. Rejects mask/output arm disagreement and writes profile markers. |
| `fit_opt_matrix_gmm.C` | Selected physical Y grid, per-run foil/slice array sizes, validated run metadata and centered guard. Matrix filename uses `nps_hms_…` or `nps_shms_…`; QA file records profile. Numerical solver, basis and coefficient parser unchanged. |
| `run_build_gmm_fit_and_matrix.sh` | Campaign validation and matching arm-specific matrix output name. |
| `diagnostics/run_diagnostics.sh` | Campaign validation before running the diagnostic chain. |
| `diagnostics/fit/plot_xptar_residuals.C` | Profile-sized sieve arrays/guides and arm-neutral angle labels; separate metadata interval correction below. |
| `diagnostics/fit/plot_yptar_residuals.C`, `plot_ytar_residuals.C` | Profile-sized column arrays/guides and angle labels; separate metadata interval correction below. |
| `diagnostics/validation/ztar/fit_ztar_peaks_from_replay.C` | Arm from run-list campaign path; selected branch/formula labels and branch checks. Optional fourth run-list field is the exact replay path; SHMS requires it. |
| `diagnostics/validation/ztar/fit_ztar_resolution_by_delta_slice.C` | Same run-list/branch work. SHMS diagnostic slices also cover 10–15, 15–20, 20–22%; these plotting bins are distinct from workflow metadata slices. |
| `diagnostics/validation/ztar/plot_ztar_stability.py` | Arm from TSV campaign path; campaign run angles take precedence over historical HMS fallback values. All target arrangements are included; three-foil targets show the outer pair plus the center foil. Generic `angle_deg` summary field and angle plot name. |
| `diagnostics/conditioning/preliminary_angular_conditioning.py` | Selected physical Y grid in event eligibility; target positions still supplied independently. |
| `diagnostics/conditioning/angular_term_ladder.py` | Passes the campaign profile into eligibility selection. |
| `diagnostics/conditioning/angular_per_foil_residual_ladder.py` | Same profile propagation; foil residual comparisons stay run-specific. |
| `diagnostics/conditioning/angular_solver_comparison.py` | Selected physical Y grid plus the SHMS beam recovery/inverse-vertex identity. Numerical comparison methods retained. |
| `diagnostics/validation/replay/probe_replay_inputs.C` | Read-only schema/finite-sample/cut-count probe for the campaign-selected arm. |
| `diagnostics/run_gema_replay_probe.sh` | Bounded probe of the four provided replay paths; default campaign selects SHMS, with no additional input needed. |
| `tests/check_workflow_syntax.py` | Reproducible fresh-interpreter ROOT loads plus Python/shell syntax checks. |
| `tests/test_replay_probe.py` | Probe routing, cut counts, invalid/missing/empty/nonfinite inputs, wrapper and read-only regression. |
| `tests/test_spectrometer_switch.py` | Profile/geometry, both-arm triple-target, mask→TFit, baseline HMS comparison, toy SVD, shifted rejection and residual interval tests. |
| `tests/test_switch_python_workflow.py` | Both cleanup entry points on both arms, target eligibility and plots, and SHMS inverse geometry. |
| `README.md`, `assistant/` navigation documents | Link the active implementation and distinguish historical copied snippets. |
| `docs/HMS_SHMS_REVIEW.md`, this file | Sources, decisions, per-file changes, validation and deferred work. |

`analyze_fit_sample_balance.C` was reviewed and remains unchanged because its generic TFit grouping already works with arbitrary foil/slice/column IDs. Existing campaign data, cuts, matrices and results are unchanged.

## Separately approved existing residual-plot bug

The user approved this correction after it was identified during the conversion.

In all three `diagnostics/fit/plot_*_residuals.C` files, `ndelcut` was read as the number of boundaries and then used as an interval count. The reader also requested one extra boundary, potentially consuming text from the next metadata entry. For a six-boundary run this produced a spurious sixth interval.

These plots now use the shared metadata reader and set their local interval count to `edges.size()−1`. No residual formula, numerical matrix fit or coefficient parser was changed by this correction. Regression tests execute all three plots with a following metadata entry present and check that histogram interval indices end at 4 for six boundaries and at 7 for nine boundaries.

The inherited matrix-reader guard is a different issue: it checks nine successfully parsed values and skips incomplete rows. More explicit malformed-coefficient diagnostics remain deferred; see the reference review.

## Running and checking

Continue using the existing commands with an `HMS_…` or `SHMS_…` campaign name. Existing explicit replay paths are required for SHMS; the legacy NPS filename fallback is not used to guess a SHMS replay. Supply an appropriate seed matrix in the campaign's existing configuration, as before.

A metadata entry uses the existing format (this is a synthetic format example, not a Gema run assignment):

```text
99001,example,12.5,3,1,9,0.0
-10,0,10
-10,-8,-5,0,5,10,15,20,22
```

That target line works identically in either arm's campaign. The final header field is legacy metadata; it does not replace the shared survey/mispointing profile.

For standalone ztar validation, place the run list and result TSV inside the named campaign. The run list retains its existing first three fields (`run opticsID comma-separated-foils`) and may append the exact replay pathname as a fourth field. SHMS needs that pathname. Do not populate actual Gema run angles, seed settings or grouping from the synthetic examples.

Local checks (ROOT/PyROOT, NumPy and the usual plotting/GMM dependencies required):

```sh
python3 tests/check_workflow_syntax.py
python3 tests/test_replay_probe.py
python3 tests/test_spectrometer_switch.py
python3 tests/test_switch_python_workflow.py
```

All listed checks passed locally. All 17 ROOT macros loaded in fresh interpreters; Python entry-point imports/help and shell syntax checks passed. Shell and fixed-theta dry runs also verified campaign routing and eight metadata intervals for the same three-foil target on both arms. The HMS regression compares every TFit column against the baseline macro on synthetic replay/mask inputs. The toy fit verifies event accounting and output routing, not calibration quality.

Once the files are installed on ifarm, the bounded real-data check is:

```sh
bash diagnostics/run_gema_replay_probe.sh > gema_replay_probe.txt 2>&1
```

The user has approved assuming −10, 0, +10 cm for all four Gema runs, and the
initial campaign is now included. Their existence and common reported
kinematics are established from user terminal output. The Mac has no direct ifarm access. Real-data checks of branches, conventions, PID, ridge/band selection and physics residuals remain outstanding. Check Gema's replay matrix/offset configuration and review the documented SHMS geometry with Holly and Mark before relying on fitted optics. Shifted-sieve support remains a later change.
