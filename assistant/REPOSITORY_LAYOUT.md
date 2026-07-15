# Repository Layout

The repository root is referred to as:

<repo>

The repository name is not assumed to be permanent.

## Active top-level areas

### <repo>/HMS_6p117GeV/

Current campaign directory for the HMS 6.117 GeV/c optics calibration.

Its stage layout is:

- 01_ytar_cuts
- 02a_angle_scan_y
- 02b_angle_scan_x
- 03a_relabel_y
- 03b_relabel_x
- 04a_candidate_trees_y
- 04b_candidate_trees_x
- 05a_gmm_cleanup_y
- 05b_gmm_cleanup_x
- 06a_fit_ntuple
- 06b_svd_fit
- issues

This directory is an example of the current campaign organization.

Do not assume that all future campaigns are HMS or use 6.117 GeV/c.

### <repo>/DATfiles/

Shared repository metadata and run information.

Important file:

<repo>/DATfiles/list_of_optics_run.dat

The relationship between this file and campaign-specific or historical run-list objects is not yet fully normalized.

Do not assume duplicated run metadata are synchronized unless explicitly verified.

### <repo>/assistant/

Version-controlled documentation for the LibreChat optics assistant.

This directory contains:

- assistant behavior instructions;
- workflow documentation;
- naming conventions;
- campaign setup guidance;
- script catalog;
- parameter guidance;
- troubleshooting notes;
- upload manifest;
- LibreChat reproduction instructions.

### <repo>/demo/

Presentation-only output collected for run group:

rg01_theta12p385_foil0

This directory mirrors the workflow stage structure for demonstration purposes.

It contains selected plots and figures used to present the workflow.

It is not a campaign template and should not be used as the authoritative source for creating a new campaign.

### <repo>/ML_angular_6p117_20260630/

Historical and operational support directory used during the 6.117 GeV/c calibration.

It contains:

- config
- cuts
- logs
- lookup
- notes
- plots
- root_out
- scripts

This directory was referenced during the calibration.

Some run-list-like or configuration content may overlap with:

<repo>/DATfiles/list_of_optics_run.dat

Treat this overlap as unresolved until the metadata sources are compared and consolidated.

Do not assume this directory is required for future campaigns without checking the current scripts.

### <repo>/logs/

Top-level logs not yet fully assigned to campaign-local stage directories.

Prefer campaign-local logs for future generalized workflows when practical.

## Top-level source files

The repository currently keeps active ROOT macros, Python scripts, and shell runners at the repository root.

Major functional groups include:

### Ytar ridge selection

- ytar_ridge_cut.C

### Y-side angle scan and assignment

- assign_yfp_ypfp_angleScanBands.C
- run_y_angle_multifoils_fixedtheta_6p117.sh
- run_relabel_y_all.sh
- relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch.C

### X-side angle scan and assignment

- assign_xfp_xpfp_angleScanBands_split.C
- run_xfp_xpfp_dynamicSplit_batch.C
- run_x_multifoils_fixedtheta_6p117.py
- run_relabel_x_all.sh
- run_relabel_xfp_from_plot_splits_v2.sh
- relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch.C

### Candidate-tree generation

- make_yscol_candidate_tree.C
- make_xscol_candidate_tree.C
- run_make_yscol_candidates_all.sh
- run_make_xscol_candidates_all.sh

### GMM cleanup

- gmm_cleanup_yscol_candidates_v2.py
- gmm_cleanup_xscol_candidates_v2.py
- run_gmm_cleanup_y_all.sh
- run_gmm_cleanup_x_all.sh

### Fit ntuple and SVD fit

- make_fit_ntuple_from_gmm.C
- fit_opt_matrix_gmm.C
- run_build_gmm_fit_and_matrix_6p117.sh

### Support scripts

- make_angle_rerun_commands.py
- rerun_open_pdf_commands.sh
- gmm_cleanup_early_y_multislice.py

Some runners are explicitly 6.117-specific and should not be treated as generalized campaign runners.

## Generated or non-authoritative areas

The following directories should normally be excluded from assistant knowledge and campaign setup guidance unless a specific recovery question requires them:

- <repo>/backups/
- <repo>/rescue_20260702/
- <repo>/script_attempts/
- <repo>/__pycache__/

These may contain obsolete, temporary, recovery, compiled, or experimental material.

## Current design direction

Future repository development should move toward:

- generalized campaign runners;
- campaign-local inputs and outputs;
- consistent run metadata;
- reduced redundancy between run lists;
- spectrometer-aware naming;
- future SHMS support;
- less dependence on momentum-specific runner names.

Until that refactor is complete, the assistant must distinguish between:

- current general scripts;
- HMS-specific implementations;
- 6.117-specific runners;
- historical or presentation-only material.
