# Parameter Guide

This guide describes the current parameters in the active Hall C optics workflow.

The repository root is referred to as:

<repo>

The campaign directory is referred to as:

<campaign>

Parameter advice must be grounded in the active source code and diagnostic plots.

Do not recommend changing several sensitive parameters at once unless the user explicitly requests a broad scan.

All delta values are expressed in percent.

## General tuning rules

Before changing a parameter:

1. Identify the specific failure visible in the diagnostic plot.
2. Identify which parameter directly controls that failure.
3. Change one or two parameters at a time.
4. Rerun only the affected foil and delta slice when practical.
5. Compare before and after.
6. Document any retained nondefault value.

Do not optimize for a visually idealized sieve grid.

The goal is to preserve the physical distortion pattern required for the optics fit.

# Ytar ridge selection

Primary script:

ytar_ridge_cut.C

## Function arguments

Current wrapper signature:

    void ytar_ridge_cut(
        Int_t nrun=1544,
        const char *tag="",
        const char *inputFileID="-1",
        const char *handFileID="-1",
        const char *inputRootOverride="",
        const char *campaignDir="HMS_6p117GeV"
    )

### nrun

Run or representative run identifier used for metadata and labeling.

Do not invent this value.

### tag

Output tag used in generated filenames.

The current campaign commonly uses the full run-group tag.

### inputFileID

Legacy or campaign-specific file identifier used when an explicit input ROOT override is not supplied.

### handFileID

Optional identifier for comparison with a hand-selected cut.

### inputRootOverride

Explicit replay ROOT-file path.

Use this when the input is a hadd run-group file or another user-specific file outside the default path logic.

### campaignDir

Campaign directory relative to the repository root.

Example:

HMS_6p117GeV

## Internal ridge settings

These are currently constants inside the macro rather than function arguments.

### deltaMin and deltaMax

Current values:

- deltaMin = -10.0
- deltaMax = 10.0

Units:

percent.

They define the full delta range considered by the ridge finder.

### deltaStep

Current value:

0.25 percent.

Controls the spacing between ridge-envelope evaluation points.

Smaller values provide finer sampling but increase computation and may follow statistical fluctuations more closely.

Larger values reduce resolution along delta.

### deltaHalfWindow

Current value:

0.35 percent.

Controls the half-width of the local delta window used to estimate each ridge point.

Increasing it gives more statistics but averages over more delta variation.

Decreasing it gives more local behavior but may become unstable in sparse regions.

### centralDeltaAbsMax

Current value:

1.5 percent.

Defines the central delta region used for reference behavior.

### centralSmoothPasses

Current value:

3.

Controls smoothing of the central ytar density used for dominant-ridge identification.

More smoothing suppresses noise but can merge nearby structures.

### minPeakFracOfMax

Current value:

0.20.

A candidate central peak must reach at least this fraction of the maximum peak height.

Lowering it admits weaker peaks and more noise.

Raising it can remove legitimate weak foil populations.

### minPeakSeparation

Current value:

0.20.

Minimum accepted ytar separation between central peaks.

Lowering it allows closer peaks to remain distinct.

Raising it can merge or suppress neighboring structures.

### widthFracLevel

Current value:

0.25.

The ridge width is estimated at this fraction of the local peak height.

Lower values generally produce wider envelopes.

Higher values generally produce narrower envelopes.

Because later GMM cleanup removes background, the ytar ridge cut should remain generous enough not to clip physical ridge events.

### smoothPasses1D

Current value:

2.

Smoothing applied to local one-dimensional ridge projections.

Increasing it may stabilize noisy width measurements but can erase shoulders or merge nearby structures.

### smoothPassesEnvelope

Current value:

8.

Smoothing applied to the ridge center and width envelopes as functions of delta.

Increasing it produces a smoother envelope but may fail to follow real local widening.

Decreasing it follows local structure more closely but may produce jagged boundaries.

### minPeakFracOfGlobal

Current value:

0.75.

Used in the ridge-strength logic relative to the global reference behavior.

Do not tune this without inspecting how the macro classifies strong and weak ridge regions.

### guardScale

Current value:

1.75.

Sets a generous upper width guard relative to reference ridge widths.

Increasing it permits wider envelopes and more background.

Decreasing it can clip legitimate local widening.

### weakRidgeFracThreshold

Current value:

0.60.

Controls when a ridge region is treated as weak.

Weak regions receive more conservative width protection.

### strongRidgeMinWidthScale

Current value:

1.10.

Minimum allowed width scale for a strong ridge region.

### weakRidgeMinWidthScale

Current value:

1.35.

Minimum allowed width scale for a weak ridge region.

Weak regions are deliberately kept wider to avoid losing legitimate low-statistics events.

### weakRidgeGuardBoost

Current value:

1.25.

Additional guard boost for weak ridge regions.

### rawWideningRetain

Current value:

0.85.

Preserves a fraction of real local widening after envelope smoothing.

### minRowEntries

Current value:

40.

Minimum local statistics required for a row-level ridge estimate.

Raising it suppresses unstable low-statistics estimates but can remove sparse delta regions.

Lowering it accepts noisier estimates.

## Ytar tuning symptoms

### Legitimate ridge tails are clipped

Consider:

- lowering widthFracLevel;
- increasing guardScale;
- increasing weakRidgeMinWidthScale;
- increasing weakRidgeGuardBoost;
- reducing excessive envelope smoothing if the widening is local.

Change cautiously.

### Envelope is too wide and includes background

Consider:

- increasing widthFracLevel;
- reducing guardScale;
- checking whether the central reference width is itself inflated;
- checking whether neighboring ridges are being confused.

Do not tighten the ridge merely to perform the job intended for later GMM cleanup.

### Ridge center is jagged

Consider:

- increasing smoothPassesEnvelope slightly;
- increasing deltaHalfWindow;
- checking local statistics.

### Neighboring foil ridges merge

Consider:

- reducing central smoothing;
- reducing excessive width padding;
- checking whether distinct peaks are being merged;
- checking the non-overlap protection already applied later in the macro.

# YFP/YPFP angle scan

Primary script:

assign_yfp_ypfp_angleScanBands.C

## Current wrapper signature

    void assign_yfp_ypfp_angleScanBands(
        Int_t nrun=1544,
        Double_t deltaMin=-8.0,
        Double_t deltaMax=-5.0,
        const char *ytarTag="ML_dev",
        Int_t maxBands=9,
        Double_t thetaStepDeg=1.0,
        Double_t minPeakSep=0.18,
        Double_t minPeakFrac=0.08,
        Double_t minPromFrac=0.25,
        Int_t smoothPasses=2,
        Double_t weakMinPeakFrac=0.005,
        Double_t weakMinPromFrac=0.04,
        Double_t weakMinSepFactor=0.45,
        Bool_t useYtarCut=true,
        Int_t foilIndex=0,
        Long64_t maxEvents=-1,
        Int_t FileID=-1,
        Int_t ndelOverride=-999,
        Bool_t runAllDeltaSlices=false,
        Double_t fixedThetaDeg=-999.0,
        const char *campaignDir="HMS_6p117",
        const char *inputRootOverride=""
    )

Important:

The current default campaignDir is HMS_6p117, while the active campaign directory is HMS_6p117GeV.

Pass the campaign explicitly rather than relying on this default.

## deltaMin and deltaMax

Delta-slice boundaries in percent.

Current standard slices include:

- -10 to -8;
- -8 to -5;
- -5 to 0;
- 0 to 5;
- 5 to 10.

## ytarTag

Tag identifying the Stage 1 ytar cut.

For the current generalized campaign workflow, this is commonly the full run-group tag.

## maxBands

Current default:

9.

Maximum number of accepted Y bands.

Lower values can truncate real structures.

Higher values permit more candidates but may admit noise if the physical geometry does not support them.

## thetaStepDeg

Current default:

1.0 degree.

Angular step used during the scan.

Smaller values provide finer angular resolution but increase runtime.

Larger values scan faster but may miss the best orientation.

When fixedThetaDeg is active, the implementation uses the fixed angle rather than the ordinary scan.

## minPeakSep

Current default:

0.18.

Minimum projected separation between strong peaks.

Lowering it permits nearby shoulders to be resolved as separate peaks.

Raising it suppresses closely spaced peaks and may merge physical structures.

## minPeakFrac

Current default:

0.08.

Minimum strong-peak height relative to the relevant maximum.

Lowering it admits weaker peaks and more noise.

Raising it removes weaker structures.

## minPromFrac

Current wrapper default:

0.25.

Minimum strong-peak prominence fraction.

Lowering it admits shoulder-like peaks.

Raising it requires peaks to stand more clearly above their local surroundings.

## smoothPasses

Current default:

2.

Number of histogram smoothing passes before peak finding.

Increasing it suppresses noise but can merge shoulders.

Decreasing it preserves local structure but may create nuisance peaks.

When a visible shoulder is being merged, reducing smoothing is usually more directly justified than globally lowering every threshold.

## weakMinPeakFrac

Current default:

0.005.

Minimum weak-peak height fraction.

This is much looser than the strong-peak threshold.

Lowering it further can admit noise.

Raising it can remove sparse physical bands.

## weakMinPromFrac

Current default:

0.04.

Minimum weak-peak prominence fraction.

Lowering it accepts weaker shoulders.

Raising it rejects ambiguous shoulders and background fluctuations.

## weakMinSepFactor

Current default:

0.45.

Scales the separation requirement for weak peaks relative to the main separation threshold.

Lower values allow weak peaks closer to stronger neighbors.

Higher values demand more isolation.

## useYtarCut

Current default:

true.

Controls whether the Stage 1 ytar cut is applied.

The normal workflow should use the ytar cut.

Disabling it is a diagnostic or development action and should be stated explicitly.

## foilIndex

Foil index within the metadata-defined foil list.

Do not substitute foil position in centimeters for this index.

## maxEvents

Current default:

-1.

A negative value means no event limit.

Use a positive value only for testing or rapid diagnostics.

## FileID

Legacy input file identifier.

Current run-group workflows often use inputRootOverride instead.

## ndelOverride

Current default:

-999.

Sentinel meaning no explicit delta-index override.

Use only when the source code and workflow require an index separate from deltaMin and deltaMax.

## runAllDeltaSlices

Current default:

false.

When false, process one requested delta slice.

When true, process the macro's standard delta slices.

## fixedThetaDeg

Current default:

-999.0.

Sentinel meaning scan normally.

A real angle forces fixed-theta mode.

Use fixed theta only when diagnostic plots show that the automatic choice is visibly worse.

## inputRootOverride

Explicit user-specific input ROOT path.

Prefer this over relying on hard-coded or legacy path construction in generalized workflows.

## Y angle-scan tuning symptoms

### Shoulder merges into a neighboring peak

First determine whether this is an angle-selection problem or a peak-identification problem.

1. Inspect the angle-scan diagnostic pages.
2. Try nearby integer or fractional angles.
3. Prefer an angle that physically separates neighboring regions before loosening peak-finding thresholds.
4. Use fixedThetaDeg when a nearby tested angle clearly preserves the separation better than the automatically selected angle.
5. Only if the regions are already visibly well separated but the peak finder fails to identify them should you tune smoothPasses, minPeakSep, minPromFrac, or weak-peak controls.

When tuning the peak finder after the angle is satisfactory:

1. Reduce smoothPasses.
2. Reduce minPeakSep slightly.
3. Reduce minPromFrac if the shoulder is visibly real but insufficiently prominent.
4. Reduce weakMinSepFactor if the structure is being considered only as a weak peak.

Do not change several peak-finding controls at once.

Avoid simultaneously changing:

- smoothPasses;
- minPeakSep;
- minPeakFrac;
- minPromFrac;
- weakMinPeakFrac;
- weakMinPromFrac;
- weakMinSepFactor.

Otherwise it becomes difficult to tell whether the physical shoulder was recovered or the finder was merely made more permissive.

### Acceptable over-splitting versus unacceptable merging

Multiple candidate divisions of one physical peak are acceptable at this stage, provided the minimum event requirement for writing cuts does not discard too many legitimate events.

Two or more physical regions merged into one candidate are not acceptable.

A single color spanning multiple physical regions indicates a merge problem.

That merge should be addressed by:

- tuning the angle at this stage;
- then tuning peak-finding parameters if the regions are already well separated;
- or, when necessary, recording the affected column for a later veto or cleanup step.

Prefer temporary over-splitting to irreversible merging.

### Too many noise peaks

Consider:

- increasing smoothPasses;
- increasing minPeakFrac;
- increasing minPromFrac;
- increasing weakMinPeakFrac;
- increasing weakMinPromFrac.

### Automatic theta is close but visibly inferior

Use a fixed theta for that foil and delta slice, and document it.

Do not force one angle across all slices without checking each diagnostic.

# XFP/XPFP angle scan and split

Primary script:

assign_xfp_xpfp_angleScanBands_split.C

## Current wrapper signature

    void assign_xfp_xpfp_angleScanBands_split(
        Int_t nrun=1544,
        Double_t deltaMin=-8.0,
        Double_t deltaMax=-5.0,
        const char *ytarTag="ML_dev",
        Int_t maxBands=9,
        Double_t thetaStepDeg=1.0,
        Double_t minPeakSep=0.18,
        Double_t minPeakFrac=0.08,
        Double_t minPromFrac=0.25,
        Int_t smoothPasses=2,
        Double_t weakMinPeakFrac=0.005,
        Double_t weakMinPromFrac=0.04,
        Double_t weakMinSepFactor=0.45,
        Bool_t useYtarCut=true,
        Int_t foilIndex=0,
        Long64_t maxEvents=-1,
        Int_t FileID=-1,
        Int_t ndelOverride=-999,
        Bool_t runAllDeltaSlices=false,
        Double_t fixedThetaDeg=-999.0,
        Double_t xpfpMinGate=-999.0,
        Double_t xpfpMaxGate=999.0,
        Bool_t writeAutoCuts=true,
        const char *campaignDir="HMS_6p117GeV",
        const char *inputRootOverride=""
    )

The shared peak-finding parameters have the same basic interpretation as the Y scan.

## xpfpMinGate and xpfpMaxGate

Current defaults:

- xpfpMinGate = -999.0
- xpfpMaxGate = 999.0

These effectively disable a restrictive X' gate.

For split processing, the low and high zones receive explicit bounds.

Use the split when the unsplit XFP/XPFP distribution contains structures that cannot be separated reliably with one global projection.

## writeAutoCuts

Current default:

true.

Controls whether provisional candidate cuts are written.

These are candidate pieces, not final trusted xscol labels.

## X-specific tuning symptoms

### One global scan merges different X structures

Use low/high X' zones or the dynamic-split driver.

Inspect both zones independently.

### A shoulder is lost in one zone

Use the same diagnostic logic as the Y scan.

First inspect nearby integer or fractional angles and prefer the angle that physically separates neighboring regions.

Only if the regions are already well separated but not identified should you:

1. Reduce smoothing.
2. Reduce minimum separation.
3. Lower prominence cautiously.
4. Adjust weak-peak controls one at a time.

Multiple candidate divisions of one physical region are acceptable if the cut-writing minimum event requirement does not remove too many legitimate points.

One candidate color spanning multiple physical regions is a merge problem and must be corrected here or recorded for a later column veto or cleanup step.

### Candidate pieces appear duplicated

Duplication is not automatically an error.

Multiple candidate pieces may later map to the same physical xscol during geometric relabeling.

### Dynamic split falls in a bad location

Inspect:

- the occupied X' range;
- low- and high-zone event counts;
- quantiles used by the dynamic driver;
- whether one zone is below its minimum event requirement.

Do not hand-select a split without recording the value and affected slice.

# GMM cleanup

Primary scripts:

- gmm_cleanup_yscol_candidates_v2.py
- gmm_cleanup_xscol_candidates_v2.py

The Y and X scripts use parallel parameter structures.

## Required arguments

### --rungroup

Full run-group tag used for campaign input/output paths and labels.

Example:

rg01_theta12p385_foil0

### --campaign

Campaign path relative to the current repository root.

Example:

HMS_6p117GeV

The scripts resolve this relative to the current working directory and require the resulting campaign directory to exist.

Run them from the project root unless explicitly using another intended layout.

### --foil

Foil index.

Required.

### --ndel

Delta-slice index.

Required.

Do not infer this value unless the slice-to-index mapping is confirmed.

## Optional selection arguments

### --yscol

Y cleanup only.

Process one yscol instead of all yscols present.

### --xscol

X cleanup only.

Process one xscol instead of all xscols present.

### --input-root

Override the default candidate ROOT input.

Default Y input:

<campaign>/04a_candidate_trees_y/root/YscolCandidates_<rungroup>.root

Default X input:

<campaign>/04b_candidate_trees_x/root/XscolCandidates_<rungroup>.root

### --tree

Defaults:

- Y: TYCand
- X: TXCand

### --outdir

Override the default output stage.

Defaults:

- Y: <campaign>/05a_gmm_cleanup_y
- X: <campaign>/05b_gmm_cleanup_x

## Cleanup controls

### --min-events

Current default:

30.

Candidate columns with fewer than this many events are removed.

Increasing it removes more sparse candidates.

Decreasing it preserves lower-statistics candidates but may retain unstable populations.

This threshold can remove an entire legitimate sparse sieve hole, so review diagnostics before raising it.

### --max-components

Current default:

9.

Upper bound on GMM components per yscol or xscol.

The scripts reject values above 9.

Too low a value can force distinct physical populations into one model component.

Too high a value can encourage nuisance components, depending on statistics and model selection.

### --keep-frac

Current default:

0.95.

Fraction retained after the GMM log-density cut.

Valid range:

greater than 0 and less than or equal to 1.

Higher values are less aggressive.

Lower values are more aggressive.

Examples:

- 0.95 retains 95 percent;
- 0.90 retains 90 percent.

A 0.90 setting rejects 10 percent and can remove sparse legitimate holes.

Use a lower keep fraction only when the resulting physical acceptance remains intact.

### --min-points-per-component

Current default:

5.

Guards against tiny nuisance components.

The scripts require a value of at least 2.

Increasing it discourages very small components.

Decreasing it allows smaller components but may model noise.

### --reg-covar

Current default:

1e-6.

Covariance regularization used for numerical stability.

Change only when there is evidence of singular or unstable covariance fitting.

### --n-init

Current default:

5.

Number of GMM initializations.

Increasing it may improve robustness at additional computational cost.

### --random-state

Current default:

13.

Controls reproducibility of stochastic initialization.

Keep fixed when comparing parameter changes.

## Plotting-only controls

The following parameters affect the diagnostic view, not the GMM selection:

- --ysieve-min
- --ysieve-max
- --xsieve-min
- --xsieve-max
- --guide-spacing
- --guide-origin
- --guide-center-index
- --guide-first
- --guide-last
- --no-guide-labels

Current plotting ranges:

- ysieve: -6.8 to 6.8
- xsieve: -12.5 to 12.5

Do not interpret a plotting-range change as a change in accepted events.

## GMM tuning symptoms

### Sparse legitimate holes disappear

Consider:

1. Raising keep-frac.
2. Lowering min-events if the entire candidate is removed before fitting.
3. Checking whether max-components is too restrictive.
4. Checking whether the candidate tree already lost the events.

### Background bridges remain

Consider:

1. Lowering keep-frac slightly.
2. Confirming the bridge is actually background.
3. Checking whether separate physical populations require more components.
4. Reviewing the candidate selection upstream.

### Runtime is high

Possible causes include:

- many columns or slices;
- high max-components;
- repeated GMM initializations;
- large event samples;
- processing all columns rather than one selected column.

Do not reduce quality parameters blindly merely to shorten runtime.

### Entire candidate column disappears

Check min-events first.

The script removes candidate columns below this threshold before meaningful cleanup can occur.

# Source-grounding requirement

When generating a command or recommending a parameter value, distinguish:

- confirmed current default;
- user-provided override;
- inferred value;
- recommended test value.

Never present a recommended test value as though it were the source-code default.
