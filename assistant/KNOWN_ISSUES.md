# Known Issues and Workflow Hazards

This file records confirmed limitations and recurring hazards in the current Hall C optics code.

The repository root is referred to as:

<repo>

The campaign directory is referred to as:

<campaign>

The assistant should distinguish confirmed code issues from suspected problems that still require inspection.

## User-specific ROOT paths

Several current runners and generated macros contain absolute paths under one user's `/volatile` directory.

Examples include:

- run_x_multifoils_fixedtheta_6p117.py
- run_y_angle_multifoils_fixedtheta_6p117.sh
- run_make_xscol_candidates_all.sh
- run_relabel_x_all.sh
- make_angle_rerun_commands.py
- rerun_open_pdf_commands.sh
- generated xrelabel_chain_<rungroup>.C files
- default arguments in some ROOT macros

Hazard:

A command may run against the wrong replay, fail because the path does not exist, or silently use stale campaign input.

Required behavior:

- ask for the user's actual ROOT-file location;
- inspect the active runner before generating a command;
- prefer explicit ROOT path or ROOTDIR overrides;
- never assume `/volatile`, `/work`, or another filesystem;
- warn when a script still contains a hard-coded path.

## Y angle-scan campaign default mismatch

Confirmed current wrapper:

assign_yfp_ypfp_angleScanBands.C

uses:

    campaignDir="HMS_6p117"

The active campaign directory is:

    HMS_6p117GeV

Hazard:

Relying on the wrapper default may direct input or output lookup to a nonexistent or incorrect campaign directory.

Required behavior:

Pass the campaign directory explicitly in generated commands until the source default is corrected.

## Campaign metadata are split across files

The current workflow uses both:

    <repo>/DATfiles/list_of_optics_run.dat

and campaign-specific files under:

    ML_angular_6p117_20260630/config/

The campaign-specific files include:

- rootdir_6p117.txt
- rungroups_6p117.tsv
- rungroups_6p117_inputs.tsv
- run_manifest_6p117.tsv
- run_foil_manifest_6p117.tsv
- runlist_6p117_use.txt

Hazard:

These files overlap in purpose and are not known to synchronize automatically.

Required behavior:

- determine which file the active script actually reads;
- do not assume one manifest updates another;
- verify run, foil, angle, delta, use-status, and ROOT path information;
- treat `DATfiles/list_of_optics_run.dat` as the primary physics metadata source only where the active code confirms that usage.

## Run groups and hadd are optional

The current 6.117 GeV/c campaign uses hadd-based run groups, but this is not universally required.

Hazard:

The assistant may incorrectly force every user or campaign to construct run groups.

Required behavior:

Ask whether the user is using:

- direct ROOT-file mode; or
- hadd/run-group mode.

Do not insert a hadd stage when each logical calibration input already has a suitable ROOT file.

## Generated X relabel chains contain absolute paths

Files named like:

    xrelabel_chain_<rungroup>.C

may embed the replay ROOT path directly.

Hazard:

A regenerated replay or another user's filesystem can leave these chains pointing at stale inputs.

Required behavior:

Before rerunning a generated chain:

- inspect its input override;
- confirm the path exists;
- confirm it corresponds to the intended replay;
- regenerate or patch the chain when necessary.

## Angle choice should be tuned before peak-finder thresholds

In YFP/YPFP and XFP/XPFP angle scans, apparent merging may be caused by a poor projection angle rather than bad smoothing or peak thresholds.

Hazard:

Loosening peak-finder parameters first can admit noise while failing to solve the physical overlap.

Required behavior:

1. Inspect nearby angle-scan pages.
2. Try nearby integer or fractional angles.
3. Prefer an angle that physically separates neighboring regions.
4. Use a fixed angle when diagnostics justify it.
5. Tune smoothing and thresholds only when the regions are already visibly separated but not identified.

## Over-splitting and merging are not equally harmful

Multiple candidate divisions of one physical region are acceptable at the provisional angle-scan stage when enough legitimate events survive cut writing.

A single candidate color spanning two or more physical regions is a merge problem.

Hazard:

The assistant may incorrectly optimize for the fewest candidates and thereby merge distinct physical structures.

Required behavior:

- prefer recoverable over-splitting to irreversible merging;
- verify that the minimum event requirement does not discard small legitimate pieces;
- allow later geometric relabeling to map multiple candidate pieces to one physical yscol or xscol;
- correct merged regions during angle or peak-finder tuning;
- record unresolved merged columns for a later veto or cleanup step.

## Minimum event thresholds can remove entire physical regions

Candidate-tree and GMM stages use minimum-event thresholds.

The GMM scripts currently default to:

    --min-events 30

Hazard:

A sparse but legitimate sieve hole or candidate piece can disappear before density cleanup is meaningfully performed.

Required behavior:

When a region disappears completely:

1. check whether it was present in the upstream candidate tree;
2. check the cut-writing minimum event requirement;
3. check GMM `--min-events`;
4. distinguish threshold removal from density rejection.

Do not assume every missing region was rejected by the GMM likelihood cut.

## GMM keep fraction can be too aggressive

The current GMM default is:

    --keep-frac 0.95

A value of:

    --keep-frac 0.90

rejects 10 percent of the selected population.

Hazard:

Aggressive cleanup can remove sparse physical holes even when the retained points are very clean.

Required behavior:

- inspect every affected foil and delta-slice diagnostic;
- raise `keep-frac` when legitimate sparse regions disappear;
- distinguish low-density physical structure from background;
- do not judge success solely by visual cleanliness.

## Candidate cuts are not final physical labels

The X angle-scan stage writes provisional candidate pieces.

Hazard:

Candidate numbering may be mistaken for final xscol identity.

Required behavior:

Do not treat provisional candidate labels as trusted physical xscol labels.

Final physical assignment occurs during geometric relabeling.

## Delta values and indices are different

The workflow uses both:

- delta boundaries in percent;
- integer delta-slice indices such as `ndel`.

Hazard:

The assistant may infer an `ndel` value from a delta range without confirming the active mapping.

Required behavior:

- always state delta units as percent;
- ask for `ndel` when it is required and not supplied;
- inspect metadata or source code before mapping a percent range to an index;
- never insert a guessed `ndel` into a runnable command.

## Script names and documentation may be historical

The repository README contains older commands, paths, macro names, and directory conventions.

Hazard:

A historically valid command may no longer match the active source.

Required behavior:

For exact commands:

- inspect the current script or macro;
- use its present argument order;
- verify current path construction;
- treat README examples as historical context unless confirmed against source.

## Successful execution does not prove physical validity

A script can complete, produce files, and still merge physical regions, discard sparse holes, use stale metadata, or fit an inappropriate event sample.

Required behavior:

Require diagnostic review at each stage.

Final matrix acceptance remains a human physics judgment based on replay and validation.
