# Troubleshooting Guide

This guide maps common workflow symptoms to the first checks that should be performed.

The repository root is referred to as:

<repo>

The campaign directory is referred to as:

<campaign>

Do not skip upstream checks. A failure visible in a later stage may have been introduced earlier.

# General command failures

## Script cannot find the campaign directory

Check:

- current working directory;
- value passed to `--campaign` or `campaignDir`;
- whether the campaign path is relative or absolute;
- known Y angle-scan default mismatch between `HMS_6p117` and `HMS_6p117GeV`.

For the GMM Python scripts, `--campaign` is resolved relative to the current working directory.

Run them from the repository root unless intentionally using another layout.

## Input ROOT file cannot be opened

Check:

- exact path;
- filename;
- filesystem location such as `/volatile` or `/work`;
- permissions;
- whether the file is empty or incomplete;
- whether a generated chain still contains a stale absolute path;
- whether the intended replay was regenerated in a new directory.

Do not substitute another ROOT file merely because its name is similar.

## Output directory is missing

Many current scripts create their own output subdirectories, but not all historical scripts do.

Check:

- campaign directory exists;
- stage directory exists;
- script path construction;
- write permissions;
- whether the script is using an old directory convention.

# Metadata problems

## Wrong foil count or foil positions

Check:

- `DATfiles/list_of_optics_run.dat`;
- campaign-specific run manifest;
- run-group composition;
- whether the representative run has the same foil configuration as all grouped runs;
- whether foil positions were confirmed interactively.

Do not infer foil positions from a run-group name alone when metadata are uncertain.

## Wrong delta slice or `ndel`

Check:

- delta boundaries in `DATfiles/list_of_optics_run.dat`;
- script convention for slice indexing;
- any explicit `ndelOverride`;
- whether the command uses percent boundaries, an integer index, or both.

Never guess `ndel`.

## Wrong angle or campaign tag appears in output

Check:

- run-group manifest;
- selected wrapper default;
- explicit `campaignDir`;
- generated filenames;
- stale scripts or chain macros;
- whether the wrong run-group tag was copied from an earlier command.

# Ytar ridge selection

## Expected foil ridge is missing

Check:

- correct ROOT file;
- correct metadata;
- foil count;
- central peak-finding range;
- `minPeakFracOfMax`;
- `minPeakSeparation`;
- statistics in the central delta region.

Before lowering thresholds, confirm the ridge is visible in the raw diagnostic.

## Ridge envelope clips real events

Check:

- `widthFracLevel`;
- `guardScale`;
- weak-ridge width protection;
- envelope smoothing;
- local statistics;
- whether the clipping is confined to one delta region.

Remember that Stage 1 should be generous because later GMM cleanup removes background.

## Ridge envelope is excessively wide

Check:

- whether the central reference width is inflated;
- whether neighboring ridges were confused;
- `widthFracLevel`;
- `guardScale`;
- excessive local widening;
- whether the apparent background is intended to be removed later.

Do not tighten the ridge solely to make the plot look cleaner.

## Two foil envelopes overlap

Check:

- central peak identification;
- foil ordering;
- width padding;
- non-overlap protection;
- metadata foil positions.

# YFP/YPFP angle scan

## Two physical regions appear as one color

This is a merge problem.

First:

1. inspect nearby angle pages;
2. test nearby integer or fractional angles;
3. choose the angle that separates the regions physically;
4. use `fixedThetaDeg` when justified.

Only if the regions are already visibly separated should you tune:

- `smoothPasses`;
- `minPeakSep`;
- `minPromFrac`;
- weak-peak controls.

If unresolved, record the affected column for a later veto or cleanup step.

## One physical region is split into several colors

This can be acceptable.

Check:

- whether all legitimate points are retained;
- whether the cut-writing minimum event threshold removes small pieces;
- whether later geometric relabeling can map the pieces to one physical yscol;
- whether any piece actually belongs to a neighboring physical region.

Prefer recoverable over-splitting to merging.

## A visible shoulder is not identified

Check in this order:

1. angle choice;
2. smoothing;
3. minimum peak separation;
4. prominence threshold;
5. weak-peak thresholds.

Do not lower all thresholds together.

## Too many noise peaks are identified

Check:

- angle choice;
- smoothing;
- `minPeakFrac`;
- `minPromFrac`;
- weak-peak thresholds;
- event statistics;
- whether the ytar cut is active.

# XFP/XPFP angle scan

## Global projection cannot separate structures

Use low/high X' zones or the dynamic-split driver.

Check:

- occupied X' range;
- dynamic split value;
- event count in each zone;
- selected angle in each zone;
- whether one zone falls below its minimum event requirement.

## One color spans several physical X regions

This is a merge problem.

Tune:

1. zone split;
2. angle;
3. smoothing;
4. peak separation and prominence.

If it cannot be resolved, record the affected xscol for a future veto.

## Several candidate pieces map to one physical xscol

This is not automatically a failure.

The geometric relabeling stage may assign multiple provisional pieces to the same physical xscol.

Check that no candidate piece crosses into another physical region.

# Geometric relabeling

## Column labels appear reversed or shifted

Check:

- foil index;
- delta index;
- metadata source;
- sieve geometry convention;
- candidate ordering;
- generated chain input path;
- whether low/high X zones are combined correctly.

Do not repair labels by hand until the geometry convention is confirmed.

## A generated X relabel chain uses the wrong ROOT file

Inspect:

    xrelabel_chain_<rungroup>.C

Check the input override near the top of the file.

Regenerate or patch the chain before rerunning.

# Candidate-tree generation

## Candidate tree has zero rows

Check:

- input ROOT path;
- tree name;
- ytar cut file;
- relabeled cut file;
- foil index;
- delta index;
- metadata;
- branch names;
- whether upstream cuts contain objects with the expected names.

## Fewer candidate rows than expected

Check:

- cut-writing minimum event requirement;
- missing provisional candidates;
- merged physical regions;
- wrong foil or delta index;
- candidate overlap logic;
- upstream relabeling.

## Overlapping X columns are reported

Investigate:

- geometry assignment;
- duplicated candidate pieces;
- low/high zone overlap;
- incorrect xscol mapping;
- whether multiple candidates intentionally map to one xscol.

Do not ignore overlap counts without examining the affected events.

# GMM cleanup

## Entire physical column disappears

Check in this order:

1. Was it present in the candidate tree?
2. Did the candidate-tree or cut-writing minimum event requirement remove it?
3. Is it below `--min-events`?
4. Did the GMM density cut reject it?
5. Is the run-group, foil, or `ndel` wrong?

## Sparse holes are removed

Consider:

- raising `--keep-frac`;
- lowering `--min-events`;
- checking `--max-components`;
- checking whether the sparse hole was already damaged upstream.

Do not assume the cleanest-looking plot is the best physical selection.

## Background bridges remain

Check:

- whether the bridge is truly background;
- `--keep-frac`;
- number of GMM components;
- upstream candidate merging;
- whether separate physical populations were forced into one candidate.

## GMM runs slowly

Check:

- number of candidates;
- event count;
- `--max-components`;
- `--n-init`;
- whether all columns are being processed;
- whether one problematic column can be tested separately with `--yscol` or `--xscol`.

Do not reduce physics quality blindly to shorten runtime.

## GMM fails with no candidate events

Check:

- `--rungroup`;
- `--campaign`;
- `--foil`;
- `--ndel`;
- optional `--yscol` or `--xscol`;
- default candidate-tree path;
- tree name;
- whether the campaign directory is resolved from the intended working directory.

# Fit ntuple

## Expected holes or delta slices are missing

Trace upstream:

1. ytar ridge selection;
2. angle scan;
3. geometric relabeling;
4. candidate tree;
5. GMM cleanup;
6. fit-ntuple inclusion or veto logic.

Do not start by changing the SVD fit.

## Fit ntuple looks clean but statistics are unexpectedly low

Check:

- GMM keep fraction;
- minimum event thresholds;
- missing candidate pieces;
- column vetoes;
- excluded runs;
- foil and delta coverage.

# SVD fit

## Fit converges but reconstructed pattern is worse

Check:

- fit sample composition;
- missing physical regions;
- merged columns;
- matrix term selection;
- source matrix;
- whether the output contains corrections or full replacement coefficients;
- replay configuration.

Numerical convergence is not sufficient.

## New matrix file gives unexpectedly large changes

Check:

- whether it contains corrections only;
- coefficient ordering;
- matrix formatting;
- source matrix combination;
- units;
- whether the replay expects a full matrix or an update.

# Replay and validation

## Ztar improves in one delta slice but worsens elsewhere

Check:

- delta coverage of the fit sample;
- missing outer slices;
- sparse foil coverage;
- over-weighting of dense regions;
- matrix term flexibility;
- whether one slice was tuned manually.

## Sieve plot looks cleaner but physical distortion is lost

This is a serious warning.

The workflow must preserve the physical distortion pattern needed by the optics fit.

Do not prefer an idealized grid merely because it looks cleaner.

# Escalation rule

When the first check does not resolve the problem:

- inspect the active source;
- inspect the exact input and output files;
- compare with the previous successful stage;
- rerun only the affected run group, foil, and delta slice;
- document the observed symptom and retained parameter change.

Do not jump directly to a broad campaign rerun.
