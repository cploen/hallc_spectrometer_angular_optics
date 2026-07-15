# Current Optics Calibration Workflow

The repository root is referred to as:

<repo>

The campaign directory is referred to as:

<campaign>

The current implementation is HMS-focused, but the workflow should be described in a spectrometer-aware and campaign-general way whenever the code permits.

## Stage 0 — Campaign inputs and metadata

Required inputs:

- spectrometer;
- campaign name;
- central momentum;
- replayed ROOT-file location;
- run or run-group definitions;
- foil configuration;
- central spectrometer angle;
- delta boundaries;
- run metadata in:
  <repo>/DATfiles/list_of_optics_run.dat

Two valid modes exist:

### Direct ROOT-file mode

Use an existing replayed ROOT file directly.

### Hadd/run-group mode

Combine multiple replayed ROOT files into one logical run-group file.

The hadd step is optional. It is required only when the user chooses to combine runs.

Quality checks:

- all input files exist;
- ROOT trees are readable;
- expected branches are present;
- run metadata agree with the replay;
- foil positions and delta boundaries are correct;
- bad runs are excluded explicitly;
- user-specific paths are recorded in campaign configuration, not assumed globally.

Expected outputs:

- campaign configuration;
- run or run-group manifest;
- optional hadd ROOT files;
- README documenting provenance.

Next stage:

- ytar ridge selection.

## Stage 1 — Ytar ridge selection

Primary script:

- ytar_ridge_cut.C

Purpose:

Select events along the reconstructed ytar-versus-delta foil ridge.

Inputs:

- replay ROOT file;
- run or run-group identifier;
- campaign;
- foil metadata;
- delta boundaries;
- ridge-finder parameters.

Typical outputs:

- ROOT cut file;
- diagnostic PDF or image;
- TSV summary.

Quality checks:

- each expected foil ridge is found;
- ridge width follows the physical population;
- sparse delta regions are not clipped;
- neighboring foil ridges are not merged;
- low-statistics tails are treated cautiously;
- all delta values are interpreted as percent.

Adjustment guidance:

- change ridge width and guard parameters only after inspecting the diagnostic plot;
- increase acceptance when legitimate ridge tails are cut;
- decrease acceptance when background or neighboring structures are included;
- adjust smoothing cautiously because excessive smoothing can merge physical structure.

Next stage:

- YFP/YPFP and XFP/XPFP angle scans.

## Stage 2a — YFP/YPFP angle scan and provisional Y-band assignment

Primary script:

- assign_yfp_ypfp_angleScanBands.C

Purpose:

Rotate the YFP/YPFP distribution so sieve-column-related structures become separated peaks or bands.

Inputs:

- replay ROOT file;
- ytar ridge cut;
- foil index;
- delta slice;
- scan range or fixed angle;
- peak-finding and smoothing parameters.

Outputs:

- provisional Y-band cut file;
- diagnostic PDF and ROOT plots;
- selected-angle summary.

Quality checks:

- expected structures are separated;
- shoulders are not silently merged;
- noise peaks are not promoted;
- sparse physical bands are not discarded;
- adjacent delta slices behave reasonably consistently;
- fixed-angle overrides are documented.

Adjustment guidance:

- use fixed theta only when the automatic scan selects a visibly inferior nearby angle;
- reduce smoothing when shoulders merge;
- increase smoothing only when noise dominates;
- modify separation, height, or prominence thresholds one or two at a time;
- preserve physical distortion rather than forcing an idealized grid.

Next stage:

- Y relabeling.

## Stage 2b — XFP/XPFP angle scan and provisional X candidate assignment

Primary scripts:

- assign_xfp_xpfp_angleScanBands_split.C
- run_xfp_xpfp_dynamicSplit_batch.C

Purpose:

Rotate the XFP/XPFP distribution and identify provisional candidate pieces for later geometric xscol assignment.

The candidate cuts are not yet trusted final xscol labels.

The current preferred approach may use a dynamic low/high X' split per delta slice.

Inputs:

- replay ROOT file;
- ytar ridge cut;
- foil index;
- delta slice;
- scan or fixed angle;
- optional X' gates;
- peak-finding parameters;
- dynamic-split quantiles and minimum event requirements.

Outputs:

- candidate-cut ROOT file;
- diagnostic PDFs and ROOT files;
- selected-angle or split information.

Quality checks:

- low- and high-zone candidate pieces are sensible;
- no obvious shoulder is merged into the wrong candidate;
- sparse candidates are not lost solely because of low density;
- the dynamic split follows the occupied X' range;
- candidate labels are not mistaken for final xscol values.

Adjustment guidance:

- use a fixed angle only when diagnostics justify it;
- change smoothing, minimum separation, height, or prominence cautiously;
- inspect both low- and high-X' zones;
- preserve candidate pieces even when multiple pieces later map to the same physical xscol.

Next stage:

- X geometric relabeling.

## Stage 3a — Y geometric relabeling

Primary scripts:

- relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch.C
- run_relabel_y_all.sh

Purpose:

Map provisional Y bands to physical yscol labels using geometry and metadata.

Inputs:

- replay ROOT file;
- Y-band cuts;
- ytar cuts;
- run metadata;
- foil index;
- delta slice.

Outputs:

- relabeled yscol cut file;
- diagnostic plots.

Quality checks:

- labels follow physical sieve ordering;
- neighboring bands are not swapped;
- each foil and delta slice is handled correctly;
- metadata come from the intended source.

Next stage:

- yscol candidate-tree generation.

## Stage 3b — X geometric relabeling

Primary scripts:

- relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch.C
- run_relabel_x_all.sh
- run_relabel_xfp_from_plot_splits_v2.sh

Purpose:

Map provisional X candidate pieces to physical xscol labels.

Multiple candidates may map to the same xscol.

Inputs:

- replay ROOT file;
- X candidate cuts;
- ytar cuts;
- run metadata;
- foil index;
- delta slice;
- low/high-zone information when used.

Outputs:

- relabeled xscol cut file;
- diagnostic plots;
- generated chain macros in some workflows.

Quality checks:

- candidate pieces are merged only when physically appropriate;
- xscol ordering follows sieve geometry;
- generated chain files do not contain stale absolute paths;
- zone-aware labels remain consistent;
- no candidate is accepted solely because it is dense.

Next stage:

- xscol candidate-tree generation.

## Stage 4a — Y candidate-tree generation

Primary scripts:

- make_yscol_candidate_tree.C
- run_make_yscol_candidates_all.sh

Purpose:

Write event-level Y candidates for later density cleanup.

Inputs:

- replay ROOT file;
- ytar cuts;
- relabeled yscol cuts;
- run metadata;
- campaign and run-group information.

Outputs:

- yscol candidate ROOT tree.

Quality checks:

- expected candidate rows are written;
- foil and delta indices are correct;
- no required branch is missing;
- output path matches the campaign structure.

Next stage:

- Y GMM cleanup.

## Stage 4b — X candidate-tree generation

Primary scripts:

- make_xscol_candidate_tree.C
- run_make_xscol_candidates_all.sh

Purpose:

Write event-level X candidates for later density cleanup.

Inputs:

- replay ROOT file;
- ytar cuts;
- relabeled xscol cuts;
- run metadata;
- campaign and run-group information.

Outputs:

- xscol candidate ROOT tree.

Quality checks:

- expected candidate rows are written;
- overlapping X columns are investigated;
- overlapping X components inside one xscol are investigated;
- foil and delta indices are correct.

Next stage:

- X GMM cleanup.

## Stage 5a — Y GMM cleanup

Primary scripts:

- gmm_cleanup_yscol_candidates_v2.py
- run_gmm_cleanup_y_all.sh

Purpose:

Remove low-density background and leakage from yscol candidate events.

Inputs:

- candidate ROOT tree;
- campaign;
- run group;
- foil index;
- delta-slice index;
- minimum events;
- maximum components;
- keep fraction.

Outputs:

- cleaned ROOT tree;
- CSV;
- summary TSV;
- diagnostic PDF.

Quality checks:

- legitimate sparse holes remain;
- dense background bridges are removed;
- rejection fraction is physically reasonable;
- component count is plausible;
- adjacent delta slices are compared;
- aggressive cleanup is not accepted without visual review.

Parameter guidance:

- higher keep-frac is less aggressive;
- lower keep-frac is more aggressive;
- raise keep-frac when sparse legitimate holes disappear;
- lower keep-frac only when obvious contamination remains;
- minimum-event and component limits should reflect statistics and expected geometry.

Next stage:

- fit ntuple, after both Y and X cleanup are complete.

## Stage 5b — X GMM cleanup

Primary scripts:

- gmm_cleanup_xscol_candidates_v2.py
- run_gmm_cleanup_x_all.sh

Purpose:

Remove low-density background and leakage from xscol candidate events.

Inputs and outputs follow the same general pattern as Y cleanup.

Quality checks:

- sparse physical sieve holes are preserved;
- candidate shoulders are not removed merely because they are low density;
- rejection fraction is reviewed;
- expected xscol populations remain;
- plots are checked before proceeding.

Next stage:

- fit ntuple, after both Y and X cleanup are complete.

## Stage 6a — Fit-ntuple creation

Primary script:

- make_fit_ntuple_from_gmm.C

Purpose:

Build the event sample used for the optics matrix fit from the cleaned Y and X selections.

Inputs:

- replay ROOT files;
- cleaned Y outputs;
- cleaned X outputs;
- run metadata;
- campaign;
- selected runs or run groups.

Outputs:

- fit ntuple;
- diagnostic plots and summaries.

Quality checks:

- accepted sieve pattern is physically sensible;
- expected foils and delta slices are represented;
- vetoed or missing regions are documented;
- input selections correspond to the intended cleanup version.

Next stage:

- SVD matrix fit.

## Stage 6b — SVD optics fit

Primary script:

- fit_opt_matrix_gmm.C

Purpose:

Fit updated optics matrix elements using the selected calibration events.

Inputs:

- fit ntuple;
- source optics matrix;
- fit configuration;
- run metadata.

Outputs:

- fitted matrix or matrix corrections;
- old-versus-new diagnostic plots;
- residual summaries.

Quality checks:

- matrix output format is understood;
- correction terms are not confused with full replacement coefficients;
- fit residuals improve without unphysical behavior;
- results are compared across foil and delta slices;
- the matrix is not accepted solely because the numerical fit converged.

Next stage:

- replay and independent validation.

## Stage 7 — Replay and validation

Purpose:

Replay calibration runs with the candidate matrix and test physical reconstruction quality.

Recommended checks:

- ztar reconstruction by foil;
- resolution by delta slice;
- sieve pattern preservation;
- comparison with previous matrices;
- central and outer foil behavior;
- angle and momentum consistency;
- independent runs not used to tune parameters when available.

A successful script run is not sufficient evidence of a successful calibration.

Final acceptance remains a human physics judgment.

## Core principle

The workflow is not merely classifying sieve holes.

It is preserving the physical distortion pattern needed for the optics fit while rejecting background and ambiguous events.
