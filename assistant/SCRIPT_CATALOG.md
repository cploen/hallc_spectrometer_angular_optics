# Script Catalog

This catalog maps workflow tasks to the current active scripts.

The repository root is referred to as:

<repo>

The campaign directory is referred to as:

<campaign>

For exact commands, inspect the current source signature or argument parser. Comments and examples may be historical.

# Stage 1 — Ytar ridge selection

## ytar_ridge_cut.C

Status:

Current active ROOT macro.

Purpose:

- apply baseline production cuts;
- locate one or more ytar ridges;
- track each ridge versus delta;
- build one TCutG per foil ridge;
- write diagnostic plots and TSV summaries.

Public entry point:

    ytar_ridge_cut(...)

Primary inputs:

- metadata run;
- output tag;
- optional legacy file identifiers;
- optional explicit replay ROOT path;
- campaign directory.

Primary outputs:

- <campaign>/01_ytar_cuts/cuts/
- <campaign>/01_ytar_cuts/plots/
- <campaign>/01_ytar_cuts/tsv/

Important behavior:

- foil 0 is assigned to the most negative ytar ridge;
- each ridge receives independent width calibration;
- hand cuts are diagnostic overlays only;
- explicit inputRootOverride should be used for run-group ROOT files or user-specific paths.

Path status:

Partly generalized.

The campaign directory and input ROOT path can be passed explicitly.

# Stage 2a — YFP/YPFP angle scan

## assign_yfp_ypfp_angleScanBands.C

Status:

Current active ROOT macro.

Purpose:

- rotate normalized YFP/YPFP coordinates;
- scan projection angle;
- identify strong and weak projected peaks;
- prefer clean separation valleys;
- permit fixed-theta overrides;
- write provisional Y-band cuts and diagnostic pages.

Public entry points:

    run_yfp_ypfp_angleScan_oneDelta(...)
    assign_yfp_ypfp_angleScanBands(...)

Primary inputs:

- metadata run;
- delta range in percent;
- ytar cut tag;
- foil index;
- angle-scan parameters;
- optional fixed theta;
- campaign directory;
- optional explicit input ROOT path.

Important behavior:

- angle should be tuned before peak-finder thresholds;
- fixedThetaDeg values above the sentinel enable fixed-angle mode;
- old cuts for the same foil and delta slice are deleted before rewrite;
- provisional band labels are not yet final physical yscol labels.

Known issue:

The current wrapper default is:

    campaignDir="HMS_6p117"

The active campaign directory is:

    HMS_6p117GeV

Pass the campaign explicitly.

Path status:

Partly generalized through campaignDir and inputRootOverride.

# Stage 2b — XFP/XPFP angle scan

## assign_xfp_xpfp_angleScanBands_split.C

Status:

Current active ROOT macro.

Purpose:

- rotate normalized XFP/XPFP coordinates;
- scan projection angle;
- identify provisional X candidate pieces;
- support low/high X' gates;
- permit fixed-theta overrides;
- write candidate cuts and diagnostic pages.

Public entry points:

    run_xfp_xpfp_angleScan_oneZone(...)
    assign_xfp_xpfp_angleScanBands_split(...)

Primary inputs:

- metadata run;
- delta range in percent;
- ytar cut tag;
- foil index;
- peak-finding parameters;
- optional fixed theta;
- optional X' minimum and maximum gates;
- campaign directory;
- optional explicit input ROOT path.

Important behavior:

- candidate labels are provisional;
- low/high zone candidates may later map to the same physical xscol;
- angle tuning should precede smoothing and threshold tuning;
- multiple candidate pieces are acceptable when they do not merge distinct physical regions.

Path status:

Partly generalized through campaignDir and inputRootOverride.

## run_xfp_xpfp_dynamicSplit_batch.C

Status:

Current active ROOT batch driver.

Purpose:

- process standard delta slices;
- estimate occupied X' range;
- choose a low/high X' split;
- run the X angle scan separately in both zones.

Public entry point:

    run_xfp_xpfp_dynamicSplit_batch(...)

Primary inputs:

- metadata run;
- run-group tag;
- campaign directory;
- explicit input ROOT file;
- foil index;
- angle and peak-finding controls;
- split quantiles;
- minimum events per zone.

Important behavior:

- dynamic split is based on occupied X' quantiles;
- each zone may select its own angle;
- slices or zones may be skipped for insufficient statistics;
- exact signature should be inspected before generating a command.

Path status:

Generalized when the explicit input ROOT file is supplied.

# Stage 3a — Y geometric relabeling

## relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch.C

Status:

Current active ROOT macro.

Purpose:

- project provisional Y-band events into HMS sieve coordinates;
- assign physical yscol labels;
- write relabeled cuts;
- generate colored-density diagnostics.

Public entry points:

    relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch_oneDelta(...)
    relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch(...)

Primary inputs:

- metadata run;
- run-group tag;
- campaign directory;
- input replay ROOT file;
- foil index;
- delta index;
- optional override file;
- geometry and minimum-event controls.

Optional override format supports:

- forced relabeling;
- rejection of bad candidate pieces;
- comments documenting manual decisions.

Important behavior:

- provisional q-band ordering is converted to physical yscol;
- geometry and metadata determine the assignment;
- exact argument order should be read from the active source before command generation.

Path status:

Generalized in the current run-group interface.

# Stage 3b — X geometric relabeling

## relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch.C

Status:

Current active ROOT macro.

Purpose:

- project X candidate pieces into HMS sieve coordinates;
- assign physical xscol labels;
- preserve multiple disconnected components for one physical xscol;
- generate low-zone, high-zone, overlay, and combined diagnostics.

Public entry points:

    relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch_oneDelta(...)
    relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch(...)

Primary inputs:

- metadata run;
- foil index;
- delta index;
- ytar tag or run-group context;
- optional override file;
- geometry and minimum-event controls;
- zone selection;
- optional X' split.

Important behavior:

- multiple candidate cuts may map to the same physical xscol;
- disconnected components are stored separately;
- downstream code must treat same-xscol components as an OR/additive collection;
- low/high X' gates must be preserved to avoid double-counting.

Optional override format supports:

- forced xscol assignment;
- rejected candidates;
- zone-specific corrections;
- comments documenting manual decisions.

Path status:

Mixed.

The macro supports overrides, but generated chain wrappers may embed absolute paths.

# Stage 4a — Y candidate tree

## make_yscol_candidate_tree.C

Status:

Current active ROOT macro.

Purpose:

- apply ytar and relabeled Y cuts to the replay tree;
- write one row per accepted event and yscol assignment;
- prepare the input tree for Y GMM cleanup.

Public entry point:

    make_yscol_candidate_tree(...)

Primary inputs:

- metadata run;
- run-group tag;
- campaign directory;
- input replay ROOT file;
- optional cut and output path overrides.

Primary output:

- TTree named TYCand;
- overview histograms by foil and delta slice.

Important behavior:

- does not require X-side row cuts;
- overlapping Y cuts may write the same event more than once with different yscol values;
- nYMatches is stored for overlap diagnostics;
- final user-facing yscol values are physical labels.

Path status:

Generalized through campaign and explicit input ROOT arguments.

# Stage 4b — X candidate tree

## make_xscol_candidate_tree.C

Status:

Current active ROOT macro.

Purpose:

- apply ytar and final relabeled X component cuts;
- preserve low/high X' zone gates;
- write one row per accepted event and xscol assignment;
- prepare the input tree for X GMM cleanup.

Public entry point:

    make_xscol_candidate_tree(...)

Primary inputs:

- metadata run;
- run-group tag;
- campaign directory;
- input replay ROOT file;
- optional cut and output path overrides.

Primary output:

- TTree named TXCand;
- overview histograms by foil and delta slice.

Important behavior:

- multiple component TCutGs for the same xscol are treated as an OR/additive collection;
- zone gates are retained to prevent double-counting;
- singleton compatibility aliases are used only when component cuts are absent;
- overlap diagnostics should be reviewed.

Path status:

Generalized through campaign and explicit input ROOT arguments.

# Stage 5a — Y GMM cleanup

## gmm_cleanup_yscol_candidates_v2.py

Status:

Current active Python script.

Purpose:

- read TYCand;
- process candidate events in projected sieve space;
- choose GMM component count by BIC;
- remove low-density events;
- write cleaned ROOT, CSV, TSV, and PDF outputs.

CLI entry:

    argparse command-line interface

Required arguments:

- --rungroup
- --campaign
- --foil
- --ndel

Important optional arguments:

- --yscol
- --input-root
- --tree
- --outdir
- --min-events
- --max-components
- --keep-frac
- --min-points-per-component
- --reg-covar
- --n-init
- --random-state

Default input:

    <campaign>/04a_candidate_trees_y/root/YscolCandidates_<rungroup>.root

Default output stage:

    <campaign>/05a_gmm_cleanup_y

Default tree:

    TYCand

Path status:

Generalized relative to the current project-root working directory.

# Stage 5b — X GMM cleanup

## gmm_cleanup_xscol_candidates_v2.py

Status:

Current active Python script.

Purpose:

- read TXCand;
- process candidate events in projected sieve space;
- choose GMM component count by BIC;
- remove low-density events;
- write cleaned ROOT, CSV, TSV, and PDF outputs.

CLI entry:

    argparse command-line interface

Required arguments:

- --rungroup
- --campaign
- --foil
- --ndel

Important optional arguments:

- --xscol
- --input-root
- --tree
- --outdir
- --min-events
- --max-components
- --keep-frac
- --min-points-per-component
- --reg-covar
- --n-init
- --random-state

Default input:

    <campaign>/04b_candidate_trees_x/root/XscolCandidates_<rungroup>.root

Default output stage:

    <campaign>/05b_gmm_cleanup_x

Default tree:

    TXCand

Known issue:

The top-of-file example still shows an older interface using:

    --run

The current parser requires:

    --rungroup
    --campaign

Treat the argparse parser as authoritative.

Path status:

Generalized relative to the current project-root working directory.

# Stage 6a — Fit ntuple

## make_fit_ntuple_from_gmm.C

Status:

Current active ROOT macro with both single-setting and batch entry points.

Purpose:

- combine replay variables with cleaned Y and X selections;
- apply optional column masks and veto rules;
- write fit ntuples;
- generate fit-sample QA plots.

Public entry points include:

    plot_fit_ntuple_qa_from_existing(...)
    make_fit_ntuple_from_gmm(...)
    plot_fit_ntuple_qa_all_existing(...)
    make_fit_ntuples_from_gmm_all(...)

Primary inputs:

- metadata run;
- run-group or setting information;
- replay ROOT location;
- cleaned Y outputs;
- cleaned X outputs;
- optics metadata;
- optional veto and mask definitions.

Important behavior:

- reads `DATfiles/list_of_optics_run.dat`;
- exact current signatures should be inspected before command generation;
- current defaults still include user-specific and 6.117-specific paths in some entry points.

Path status:

Mixed.

Supports arguments, but some defaults remain campaign- and user-specific.

# Stage 6b — SVD optics fit

## fit_opt_matrix_gmm.C

Status:

Current active ROOT macro, but strongly 6.117-specific.

Purpose:

- read fit ntuples;
- construct and solve the optics matrix system with SVD;
- compare old and fitted reconstruction;
- write coefficient and diagnostic outputs.

Public entry point:

    fit_opt_matrix_gmm(...)

Current major arguments:

- tag;
- FileID;
- maximum number of fit events;
- number of settings;
- input fit-ntuple directory;
- output directory;
- old coefficient file.

Current hard-coded content includes:

- maximum foil and delta counts;
- synthetic metadata run list 611701 through 611708;
- 6.117 campaign defaults;
- default old matrix file.

Important behavior:

- this macro is not yet campaign-general;
- inspect coefficient-writing behavior before treating output as a full replacement matrix;
- do not assume convergence proves physical validity;
- replay and independent validation are required.

Path status:

Campaign-specific and partially hard-coded.

# Batch runners

The repository also contains shell and Python runners for processing all run groups.

Examples include:

- run_y_angle_multifoils_fixedtheta_6p117.sh
- run_x_multifoils_fixedtheta_6p117.py
- run_relabel_y_all.sh
- run_relabel_x_all.sh
- run_make_yscol_candidates_all.sh
- run_make_xscol_candidates_all.sh
- run_gmm_cleanup_y_all.sh
- run_gmm_cleanup_x_all.sh
- run_build_gmm_fit_and_matrix_6p117.sh

Status:

Current operational helpers, but many are 6.117-specific.

Hazards:

- hard-coded ROOT paths;
- hard-coded campaign names;
- campaign-specific run-group TSV defaults;
- stale assumptions after replay regeneration.

Required behavior:

Inspect each runner before reuse in another campaign.

# Generated chain macros

Files named:

    xrelabel_chain_<rungroup>.C

Status:

Generated operational files.

Purpose:

Run a configured X relabeling chain for one run group.

Hazard:

They may embed absolute replay ROOT paths and stale configuration.

Required behavior:

Inspect or regenerate before reuse.

# Historical examples and comments

Some source comments still contain:

- older macro names;
- older output paths;
- pre-campaign directory conventions;
- deprecated command-line options.

The current public function signature or argparse parser is authoritative.
