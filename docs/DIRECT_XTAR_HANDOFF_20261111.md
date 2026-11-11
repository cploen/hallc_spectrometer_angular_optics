# Continue the HMS 2018 separate x-target predictor experiment

November 11, 2026. Handoff from the Mac Codex session to ChatGPT in the Linux
Mint browser, with Omnissa visible as another window.

## Read this first: the intended experiment

**Keep the complete existing 2018 transport functions fixed. Fit a separate
event-level x-target predictor, with geometric `xtarT` as its training target.**

An earlier agent interpretation was wrong: it refitted x-target-dependent
coefficients inside the three angular/position functions. That run was stopped
after the user explicitly rejected that interpretation. **Do not resume
`run_xtar_refit.sh` or treat that conditional fit as the requested experiment.**

No separate x-target predictor has been implemented, fitted or installed yet.
The user has now approved proceeding with the separate offline predictor and
specifically requests an implementation/SVD/coordinate-consistency audit.

## User requirements and approvals

1. Goal: improve the 2018-derived reconstruction before deciding how to refit
   NPS optics data with the improved starting reconstruction.
2. Preserve the original 2018 matrix as `oldfit.dat`. Do not overwrite it with
   the completed new 2018 fit or with the new x-target predictor. Confirm its
   provenance before giving a final comparison the label "original 2018".
3. Final comparison must show **original 2018**, **completed 2018 refit**, and
   **completed 2018 refit with the separate x-target predictor**, across all
   four requested variables: **xptar, yptar, ytar, xtar**.
4. All transport coefficients/functions stay fixed when adding the predictor.
   Feeding a different x-target input into those functions can change their
   outputs numerically. This is the effect the user wants to inspect.
5. Use the accepted event selections, exact IDs, masks, exclusions and existing
   80% core / 20% inner-shoulder training mixture. This is a population mixture,
   **not** the training/holdout split.
6. Select a limited predictor basis with training-only elastic net and beam
   search, **beam width 30**. Use scaled direct-X SVD for unpenalized refits.
   Preserve train-only scaling and the saved folds. No normal-equation solve
   is needed for comparison.
7. Explicitly double-check the x-target SVD/design/units/coordinate construction
   against the other established SVD paths and the actual HCANA/replay source.
   An earlier missing-mispointing error in XPTAR was a real bug; do not repeat
   such a mismatch. Prepare a factual review handoff for another agent and,
   if warranted, Mark and Charles. No messages to them have been sent.
8. The user chose to **skip afterburner closure as a prerequisite**. Do not
   reinstate reproduction of the unavailable historical 2018 analysis
   environment as a gate. Keep the known closure limitation visible. This does
   not waive the explicitly requested formula/implementation checks.
9. Ultimately validate through matched baseline/candidate HCANA replays. No
   matrix or new predictor should silently be installed into replay now.
10. Keep the workstation/session awake while actively working; Mac wakefulness
    does not prevent the remote RHEL session from locking independently.

## Required residual heatmaps

Compare identical reserved events, keeping original selection and delta
coordinates/bin assignments fixed. Do not reselect events with the candidate.
Show all accepted physical foils and these edges in percent:

```text
-15, -12, -10, -8, -5, 0, 5, 8, 10, 12, 14, 15
```

In particular, do not combine -15 to -10 into one bin. Retain separate
[-15,-12) and [-12,-10) cells. Keep the extended positive-delta intervals.
The usual saved edges were -15,-10,-8,-5,0,5,8,10,15; the extra diagnostic
edges can be added without changing selection or allocation.

Report bias, RMS including bias, change in RMS and counts for xptar, yptar,
ytar and xtar. Show both absolute baseline/candidate values and clearly signed
changes, so improvement cannot be confused with a larger error. Separate core
and supported noncore. Count unsupported/blocked populations explicitly, flag
sparse cells and leave uncovered cells blank. The user accepts very wide
images and specifically asks for adequately spaced/readable text.

Existing inspected holdouts are development evidence, not a fresh independent
test. A whole-run or setting stress check can help identify a model that learns
only the finite foil/angle pattern. Foil/hole labels must not be input features.

## What the proposed predictor is, and is not

The standard transport matrix takes four focal-plane coordinates **plus an
externally supplied x-target input**, and outputs xptar, ytar, yptar, delta.
The fourth output is **delta**, not xtar. There is no independent x-target
column to refit in this matrix format.

The approved new experiment is therefore a **separate model**, stored alongside
the matrix with its own input names, units, coefficients, source hashes and
validity range. It must not be inserted into oldfit.dat's delta column.

Proposed input contract: use focal-plane measurements and measured beam
information available during ordinary replay. Known central settings such as
angle may be used explicitly if needed; define and document the chosen basis
before fitting. No true foil position, foil ID, sieve-hole ID, xptarT, xtarT or
other truth-derived values may be predictor inputs. Geometry labels define
the training target and diagnostic truth only.

One proposed physical parameterization is a measured-beam/pointing anchor
plus a learned longitudinal-projection contribution. Its precise monomial
basis/order, angle dependence and intercept convention are **not yet chosen**.
Do not claim those details were already approved or implemented. The user
approved the model class and the controlled fitting/validation procedure.
Check that the data constrain the predictor; many events alone do not prove a
unique inverse from focal-plane coordinates to x-target.

## Production paths on ifarm2402

```bash
A=/u/group/nps/cploen/ML_HMS_OPTICS_JOINT_SURFACE
C="$A/HMS_LowMom_LowDelta_2018_Accepted"
```

Frozen sample tag: **inner20**.

- Selection/allocation: `$C/05c_core_sample/inner20/`.
- Training trees: `$C/06c_core_ntuple/inner20/fit/root/`.
- Reserved trees: `$C/06c_core_ntuple/inner20/holdout/root/`.
- Saved training folds: `$C/06e_beam_search/beam30_inner20/tsv/folds.tsv`.
- Original-matrix candidates to verify:
  `$C/config/oldfit.dat` and `$C/06c_core_ntuple/inner20/fit/oldfit.dat`.
- Completed 2018 raw Beam result:
  `$C/06e_beam_search/beam30_inner20/matrices/beam.dat`.
- Actual zero-offset export used in the completed new-coefficient NPS replay:
  `/u/home/cploen/nps_replay_coeff2018_zero_20261110/DATFILES/hms_recon_coeff2018_new_zero_offsets.dat`.

Use care with the completed-refit baseline: its source Beam fit had a free
YTAR intercept (an earlier mistake), then the **entire active 00000 row was
removed after the fit** for the approved zero-offset replay experiment. Other
coefficients were unchanged. All selected external HMS offsets were zero.
Do not silently restore those constants or describe this export as a fit that
had constrained every intercept during optimization. Future YTAR fits must
explicitly fix ytar_00000=0.

SHA256 of that exact zero-offset replay export:

```text
f715879845f777f2da2aa62b7e66756907929bdf4b7b40809ac55370156f8801
```

The preview file
`$A/HMS_LowMom_2018_SievePreview_20261110/config/beam_zero_offsets.dat`
has a different byte hash:

```text
a19a894ffb6310ee603a57266bdb4f83c881129fb900b411e7ad56371ba7064
```

The two row sets and every parsed coefficient compared **exactly equal** on
ifarm. Their byte representations differ; use explicit hashes/provenance.
The last actual run used the first file, matching the replay snapshot.

The accepted sample contains **118,850 training events**, with 23,770 per
physical foil across five foils. The previous preflight verified the complete
frozen membership and geometric xtarT equations on the real sample. The user
describes broad 2018 coverage across three angles, two momentum settings and
many runs. Tabulate actual accepted coverage rather than assuming every
nominal run/foil contributes. In particular, run 1347's -10 cm foil was excluded.

## Verified event-coordinate chain

In `make_fit_ntuple_from_gmm.C`:

- `xtar` is copied from replay `H.gtr.x` (around line 494).
- `reactycalc` stores replay `H.react.y` (496, 558).
- `ztar` stores saved replay reaction z.
- `xtarT` is calculated using the assigned foil/hole geometry and is exported
  separately (545, 626–630).
- `reactxcalc` and raster/BPM inputs also exist in the export/source chain;
  verify the actual available branches before choosing predictor inputs.

Current geometry, positions cm and slopes rad:

```text
L = 168 cm
x_vertex = -reactycalc - x_mispointing
z_vertex = zfoil * cos(angle)
xptarT = (xsT - x_vertex) / (L - z_vertex)
xtarT  = x_vertex - xptarT * z_vertex
```

The HMS pointing function currently used by the fit code is
`0.1*(2.37 - 0.086*d + 0.0012*d*d)` cm, where
`d=min(abs(angle_deg),50)`. Compare the actual replay parameters/source rather
than assuming this profile exactly describes all historical binaries.

Both `xtar` and `xtarT` are stored in cm; polynomial position inputs are in m.
Existing `fit_elastic.design` uses `xtar/100` as its fifth input. COSY supplies
the x-dependent coefficients; it does not supply an event-level geometric
training coordinate independently of the data/geometry.

The original Python `load_sample` dropped xtarT from its retained fields. The
published conditional extension adds an explicit `xtar_coordinate='xtarT'`
option that retains it and checks the equation above against saved geometry.
This loader option can be reused for the new direct predictor without using
the superseded conditional fit runner.

An actual run-1338 training file was checked:

```text
$C/06c_core_ntuple/inner20/fit/root/Optics_1338_-1_fit_tree_gmm.root
```

It contains 21,058 events; RMS(xtar-xtarT) is
0.08768237814780522 cm (about 0.877 mm). This includes bias and is **not** an
independent resolution measurement or proof of a bug.

## Installed HCANA source checks already made

Source directory inspected directly on ifarm:

```text
/u/group/halla/apps/hcana/1.1/src/src
```

`THcHallCSpectrometer.cxx` around lines 519–567:

- x-target is passed into `CalculateTargetQuantities`.
- `hut[4] = xtar/100.0` is the fifth polynomial input (around 537).
- `delta = sum[3] + fDeltaOffset` is the fourth output (around 567).

`THcExtTarCor.cxx` around 117–141:

- It obtains the beam/vertex transport anchor from `-vertex[1]-pointing_off[0]`.
- It reevaluates the transport functions while updating
  `xtar_new = x_tg - xptar*ztarg*cos(theta)` in cm.
- Up to five updates, with a 2 mrad xptar-change stopping condition.

The exact variable binding/stage represented by historical `H.gtr.x` has not
been certified. Nor has closure with historical 2018 replay been resolved.
The old afterburner evaluates saved vertices and has a known closure warning;
it must not be presented as a new HCANA replay.

Upstream references are useful leads, not substitutes for installed versions:

- https://github.com/JeffersonLab/hcana/blob/master/src/THcHallCSpectrometer.cxx
- https://github.com/JeffersonLab/hcana/blob/master/src/THcExtTarCor.cxx
- https://github.com/JeffersonLab/hcana/blob/master/src/THcReactionPoint.cxx

Check focal-plane translations/rotations/header parameters as well as the
polynomial loop. Check angle signs, pointing units and external offsets exactly
once. Do not infer unknown offsets by centering residual distributions.

## Source repository and implementation status

Repository: **https://github.com/cploen/hallc_spectrometer_angular_optics**.

Existing production source checkout on ifarm:
`/u/home/cploen/xtar_refit_2018_20261111`, detached at
`b184e5ad6c0d96c91135948b429aabae2a9c19a8`.

Published branch `codex/xtar-conditional-refit-20261111` contains the
**superseded conditional-coefficient runner**. Useful reusable pieces there:

- `fit_elastic.load_sample(..., xtar_coordinate='xtarT')`: coordinate retention
  and validation, preserving all existing membership checks.
- `elastic_net.scaling(..., fit_intercept=False)` and `tune(..., fit_intercept=False)`:
  optional uncentered no-intercept mode; old entry-point defaults are unchanged.
- `beam_search.Folds` and `fit_all` accept the same optional no-intercept mode.
- Direct-X/orthogonal-QR compression is used; no normal-equation solve.

**Important new-target implementation detail:** the existing numerical engines
still hardcode the three old TARGETS and SCALES=(1000,100,1000). A one-target
xtar predictor needs an explicit single-target interface (or a mathematically
equivalent dedicated wrapper), with correct native/output units. Do not relabel
an xptar fit as xtar, duplicate a response three times, or accidentally use
mrad scaling for an x-target position residual. Preserve defaults and test the
new path against the established scaled-X construction.

Tests passed for the superseded extension: geometric-coordinate rejection,
zero-intercept direct-X equivalence, frozen-block/ID preservation, and existing
centered solver/scaling behavior. They do not validate a new direct predictor
that has not yet been written.

A clean local Mac worktree was created at
`work/hallc-xtar-direct-20261111` on branch
`codex/xtar-direct-predictor-20261111`, initially based on b184e5a. Only handoff
documentation will be added/published there for this transfer. That branch
does **not** contain a completed direct predictor.

## Stopped run: preserve, do not resume

The wrong conditional fit ran on PID **2790948**. SIGTERM was sent and `ps`
confirmed it terminated. Partial outputs are preserved at:

```text
$C/06g_xtar_refit/xtarT_order4_beam30/
```

`SUPERSEDED.txt` states why the run was stopped. Its manifest may still say
running because SIGTERM interrupted it; the process check and note supersede
that status. No matrices from it were installed/adopted. Do not delete or
overwrite the evidence. Its log is:

```text
/u/home/cploen/xtar_refit_2018_20261111/xtarT_order4_beam30.log
```

## Recommended next actions in the Linux Mint session

The browser ChatGPT session may not have native computer-use tools. In that
case, have the user paste terminal commands and provide concise output/photos.
Do not assume it can control Omnissa merely because the window is open.

1. Establish access to ifarm and the accepted sample, then verify the original
   `oldfit.dat` provenance and the completed-refit export. Do not run a fit yet.
2. Read the source files listed above and audit target/feature units, sign,
   pointing correction, external offsets, matrix header transforms and the
   exact construction of the one-target SVD problem.
3. Define the compact predictor basis and available measured inputs explicitly.
   Preserve the exact data partitions. Implement the separate predictor with
   a clear human-facing CLI and a separate export. Do not use the conditional
   coefficient runner.
4. Run a read-only real-data preflight, then training-only EN/beam selection.
   Validate scaled-X SVD, exported prediction identity and unit/pointing
   consistency with targeted meaningful checks. Avoid exhaustive testing.
5. Generate the three-model, four-variable extended heatmaps described above
   on common event IDs. Label any offline vertex/iteration assumptions.
6. Pause for user review before replay installation or wider NPS refits. A
   representative matched 2018 replay may be requested from Josh later.

Read-only starting commands (ifarm terminal):

```bash
A=/u/group/nps/cploen/ML_HMS_OPTICS_JOINT_SURFACE
C="$A/HMS_LowMom_LowDelta_2018_Accepted"
ls "$C/06c_core_ntuple/inner20/fit/root"
sha256sum "$C/config/oldfit.dat" "$C/06c_core_ntuple/inner20/fit/oldfit.dat"
sha256sum /u/home/cploen/nps_replay_coeff2018_zero_20261110/DATFILES/hms_recon_coeff2018_new_zero_offsets.dat
```

Do not send messages to others without explicit authorization. If an
implementation inconsistency is found, produce a review note with matched
event IDs, source/matrix/flags hashes, expected and observed formulas/numbers,
and a focused question for Mark and Charles. A difference between xtar and
xtarT alone is not a proven analyzer defect.

## Other replay work, for context only

- Sieve workflow `hms_sieve_coeff2018_new_zero_offsets_20261110` (268650):
  all **45 succeeded**, as reported/observed on November 11.
- Separate high-momentum elastic workflow `hms_elastic_pcentral_20261110`
  (268691): 11 jobs were corrected/retried after a log-path launch failure.
  Later completion is **not verified in this handoff**. It uses a different
  high-momentum calibration branch, not the coeff2018 sieve matrix.
- Prior packaging bugs resolved `/mss` and `/farm_out` symlinks into namespaces
  rejected by SWIF2/Slurm. Preserve the literal `/mss` tape URI namespace and
  `/farm_out` log prefix; do not use Path.resolve() on those service paths.

The remote CUE RHEL9 session is currently locked. The user is away from the Mac
and is moving to Linux Mint. Do not attempt to bypass the lock.
