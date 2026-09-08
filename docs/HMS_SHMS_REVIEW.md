# HMS/SHMS campaign switch: reference review and implementation record

Review date: 2026-09-08. **Status: centered-sieve implementation and initial SHMS campaign inputs; software validation recorded below. Real SHMS replay validation remains pending.**

## Checkout and inherited matrix-reader fix

The supplied ifarm checkout is `/u/group/nps/cploen/ML_HMS_OPTICS_DEV`, on
`prelim-6p667-data` at `5cb66eee3fca4fbe485448cf2c07801bbce91943` with no local
changes shown by `git status`. It already contains GitHub `main` at
`ea832dc57d32445fd950b021b7013e9a35ee7211`; switching to main would go backward
relative to this work. At inspection, GitHub `prelim-6p667-data` is three commits
ahead, at `057b908956286099c49a8d453856b0a6feb846c7`.

The separate Mac checkout is `work/hallc-angular-hms-shms`, on local branch
`hms-shms-campaign-switch`, based on that latest preliminary branch. The older
Mac checkout and its uncommitted changes were left untouched. No push or merge
was performed during the initial review. The user subsequently requested branch
publication before testing; use `hms-shms-campaign-switch` for the current
implementation. The user was given the following preliminary-branch update; its
execution has not yet been confirmed:

```sh
cd /u/group/nps/cploen/ML_HMS_OPTICS_DEV
git pull --ff-only origin prelim-6p667-data
git status --short --branch
git rev-parse HEAD
```

Commit `1e7afe670c6e6192bce7a541caddc88cea1fcc7d` already changes the seed reader
in `fit_opt_matrix_gmm.C` to check that `sscanf` read all four coefficients and
five exponents (`parsed == 9`) before appending a matrix row. Previously a blank
line could append uninitialized values, or values retained from a preceding
iteration. The commit documents a blank line in the July 2026 seed matrix as the
trigger, with machine-dependent retained-xtar contributions to the fit.

This is a reader robustness fix, **not a validation of the matrix's physics**.
A blank line is not itself evidence of a bad matrix. The guard silently skips
any line with fewer than nine successfully parsed values, so a malformed term
can also be dropped without a diagnostic. It does not check finiteness,
duplicate exponents, extra trailing data, completeness, or correct spectrometer.

Proposed further improvement, discussed but **not implemented**: explicitly
skip allowed whitespace/comment lines, and fail with filename and line number
on malformed coefficient rows. This is separate from HMS/SHMS routing and
must preserve the valid seed format and recognized separators. The inherited matrix-coefficient reader is unchanged by the spectrometer work.
The separately approved residual-plot metadata correction is described below.

## Physics conventions and evidence

Campaign naming selects a fixed HMS or SHMS profile without a new
`--spectrometer` argument. Survey constants are not per-run tuning parameters.
The user clarified that new mispointing surveys are not normally performed.
Different literals in analysis code do not establish that the survey changed.

| Quantity | Existing HMS workflow | SHMS reference | Implementation consequence |
|---|---|---|---|
| Replay prefix | `H.` | `P.` | Route target, focal-plane, vertex, raster and sieve branches together. |
| Electron Cherenkov observable | `H.cer.npeSum` | `P.ngcer.npeSum` in Holly's delta and later angular codes | Do not mechanically create `P.cer.npeSum`, or silently substitute `P.hgcer.npeSum`. |
| Target-to-sieve distance | 168 cm | 253 cm | Centralize the geometry. |
| Sieve x rows | 9, spacing 2.54 cm | 11, spacing 2.5 cm | Update allocation, loops, labels, plots and fit bookkeeping. |
| Sieve y columns | 9, spacing 1.524 cm | 11 centered, 10 shifted, spacing 1.64 cm | Use the existing sieve flag; the arm name cannot identify the inserted sieve. |
| Centered sieve coordinates | x=(i−4)×2.54, y=(j−4)×1.524 | x=−12.5+2.5i, y=−8.2+1.64j cm | Physical column IDs must match the selected geometry. |
| Shifted SHMS sieve | No SHMS support | y=−8.2+0.82+1.64j cm for j=0…9 | Deferred. Explicit `SieveFlag=2` is rejected by centered-geometry stages. |
| In-plane vertex geometry | Positive foil-z sine term in the present HMS routine | Negative foil-z sine term with positive SHMS central angle | A distance-only switch is insufficient. |
| Horizontal-bender correction | None in the HMS truth routine | Momentum-dependent y-sieve displacement | Apply the same convention to truth and reconstructed-sieve checks; avoid double correction. |
| Mispointing | Existing angle-dependent HMS parameterizations | Survey note: xMP=−0.126 cm, yMP=−0.05 cm | Use a documented SHMS survey profile. Different analysis literals need attribution, not silent adoption. |
| Nominal momentum acceptance | Current workflow generally uses −10%<δ<10% | JLab specification: −10%<δ<22% | Preserve the distinction between nominal acceptance, analysis cuts and slice boundaries. |

SHMS geometry and sieve definitions are from
[Holly's configuration code](https://github.com/hszumila/SHMS_optics/blob/c7d70e88b0cf2d1d93359d941168ab8f35ab916e/source/src/myConfig.cpp).
The survey constants are also documented in section 5.4 of
[Holly's February 2019 optics note](https://hallcweb.jlab.org/DocDB/0010/001007/001/shmsNote.pdf).
The nominal acceptance is documented by
[Jefferson Lab](https://www.jlab.org/Hall-C/upgrade/HALLC_12GEV/shms_beam_envelope.html).

### SHMS truth equations to review with Holly and Mark

The event-truth calculation in
[SHMS_optics](https://github.com/hszumila/SHMS_optics/blob/c7d70e88b0cf2d1d93359d941168ab8f35ab916e/source/shms_optics.cpp#L800)
provides a coherent reference. With lengths in cm, theta in radians, and delta
in **percentage points**, its equations are:

```text
Xv = −react_y − xMP
Yv = −zfoil*sin(theta) + react_x*cos(theta) − yMP
Zv =  zfoil*cos(theta) + react_x*sin(theta)
C(delta) = −0.0398*delta + 0.000398*delta^2  [cm]
xptar_true = (xsieve_hole − Xv) / (253 − Zv)
yptar_true = (ysieve_hole − C(delta) − Yv) / (253 − Zv)
xtar_true = Xv − xptar_true*Zv
ytar_true = Yv − yptar_true*Zv
```

The correction expands the reference's coefficients
`−0.019δ + 0.00019δ² + 40(−0.00052δ + 0.0000052δ²)`.
The corresponding closure checks are
`xtar_true+253*xptar_true=xsieve_hole` and
`ytar_true+253*yptar_true+C(delta)=ysieve_hole`.
These equations are implemented in `spectrometer_config.h::targetTruth` and used
by `make_fit_ntuple_from_gmm.C`. Replay `extcor.ysieve` is read as supplied; the
HB correction is applied to the truth projection, not applied again to that branch.

Reference disagreements to document rather than copy indiscriminately:

- `SHMS_optics` uses xMP=−0.126, yMP=−0.05 cm; `cafe_optics` sets
  xMP=−0.126, yMP=−0.03 cm; `shms_optics_a1n` sets both to zero. The latter
  values are code settings, not evidence of replacement surveys.
- `cafe_optics` additionally reconstructs tracks from a separate matrix and
  uses both replay delta and recalculated delta in its truth expression.
  That complete analysis is not a drop-in replacement for this pipeline.
- Some older metadata use slice centers and widths, whereas this repository
  uses boundary arrays. Do not transfer those tables without conversion.
- The 2017 correction option in `SHMS_optics` changes foil-selection
  coordinates. Do not enable it solely because a run is from 2017; first
  establish what the supplied replay already reconstructs.

### What the cuts mean

There are distinct types of cuts: detector cuts select electron-like events;
foil cuts identify the target foil; sieve-band cuts identify the hole;
momentum cuts and delta slices define the portion of phase space fitted.
Changing one type is not equivalent to changing another.

The current ridge finder uses Cherenkov `npeSum>2`, and no calorimeter cut,
inside −10%<δ<10%. Most later HMS stages use `npeSum>6` and
`cal.etottracknorm>0.65`. Holly's `HMS_OPTICS_HIGH_MOM`, `deltaOpt` and later
SHMS angular routines include `npeSum>6` and `cal.etottracknorm>0.8`.
Her January 2019 setup instead uses `P.cal.etracknorm>0.8` and a
−12%<δ<20% window. `etracknorm` and `etottracknorm` are different branch
names and must not be treated as interchangeable without checking their
definitions. Some later angular selection code uses −15%<δ<24%.

The SHMS switch should support its nominal range automatically. Those other
windows are analysis selections and do not redefine the hardware acceptance.
The replay probe compares the current threshold, the stricter threshold, the
HMS-sized window, the nominal SHMS window and the wider historical window.
Event counts alone do not prove particle purity or determine the optimal cut.

## Implemented behavior and file-by-file changes

See [HMS_SHMS_CHANGES.md](HMS_SHMS_CHANGES.md) for the complete per-file record,
usage, target independence, and the separately approved diagnostic correction.
All root-level executable workflow scripts and relevant diagnostics now use
campaign-selected profiles. `analyze_fit_sample_balance.C` was inspected and
needs no changes: it groups generic TFit foil/slice/column IDs without assuming
nine columns. The matrix basis, numerical solver, existing PID thresholds,
and matrix-coefficient parser are unchanged.

The user explicitly approved assuming a **centered SHMS sieve** for now and
adding a position option later. No new spectrometer or sieve-position prompt
is introduced. A metadata flag that explicitly requests a shifted SHMS sieve
is rejected rather than interpreted as centered.

**Foil positions are independent of spectrometer.** The −10, 0, +10 cm
three-foil target is available to both HMS and SHMS; it was simply not used
in the NPS HMS campaigns. The arm profiles contain no target positions.
The optics metadata supplies foil coordinates and delta boundaries for each
run. SHMS ridge IDs account for the reversed ytar ordering at a positive
spectrometer angle, without changing those foil coordinates.

Historical campaign outputs, matrices and cuts are untouched. New
`SHMS_8p5695GeV` inputs and four appended optics metadata entries are
documented in that campaign's [README](../SHMS_8p5695GeV/README.md). Copied source
snapshots in `assistant/` are identified as historical; executable files and
this change record take precedence. The older Mac checkout's untracked
`angular_term_column_correlations.py` is not in the new checkout's GitHub
baseline and has not been silently imported.

## SHMS replay check and remaining run conditions

The four user-supplied paths are under
`/volatile/hallc/c-deuteron/gvill/ROOTfiles/prod/`:

```text
deut_replay_prod_3283_-1.root
deut_replay_prod_3284_-1.root
deut_replay_prod_3285_-1.root
deut_replay_prod_3286_-1.root
```

The supplied run description identifies 2017 sieve data on the −10, 0, +10 cm
three-foil target.
The user clarified that the target positions are −10, 0, +10 cm and that
this target is available for either spectrometer. The user subsequently
authorized assuming these coordinates for all four runs; initial processing
does not wait for another target confirmation.
Angles, momenta, beam energy, particle mass hypotheses and target mass agree across reports for runs 3283–3286, as documented below. The user also supplied an ifarm listing confirming that all four named ROOT
files exist (842M, 791M, 847M, 709M respectively; displayed July 7 timestamps).
Their contents have not been inspected here. Actual sieve position (centered
is the approved assumption), replay matrix and offset configuration remain
unverified. We have no
direct ifarm access. Do not construct invented rungroup rows or overwrite an
existing campaign from these assumptions. The new campaign uses reported
angles/momenta and the explicitly approved foil/sieve assumptions; the delta
boundaries are documented initial analysis bins.

### Runs 3283–3286 report evidence

Source: terminal output pasted by the user from
`/volatile/hallc/c-deuteron/gvill/REPORT_OUTPUT/prod/deut_prod_{3283,3284,3285,3286}_-1.report`.
The `SW_` fields agree across all four reports. For 3283, the repeated
general-run fields also agree with its `SW_` section. The reports have not
been accessed directly from this Mac. Every value in the table applies to
all four runs.

| Report field | Value | Interpretation |
|---|---|---|
| Beam energy | 10.600000 | 10.6 GeV replay setting |
| SHMS central momentum | 8.569500 | 8.5695 GeV/c replay setting |
| SHMS angle | +8.915000 | Positive angle, consistent with the implemented SHMS convention |
| SHMS particle mass | 0.000511 | Electron mass hypothesis |
| HMS central momentum | 2.194000 | 2.194 GeV/c replay setting |
| HMS angle | −54.960000 | Preserve the reported sign in the run record |
| HMS particle mass | 0.938272 | Proton mass hypothesis |
| Target mass | 12.010700 amu | Carbon mass setting; does not identify foil count or positions |

The particle mass entries describe replay configuration, not a PID-purity
measurement. In particular, this is not automatically a suitable electron-PID
validation sample for the HMS workflow. The broad 3283 report search did not return foil coordinates, sieve position,
matrix filename, or mispointing/offset settings. The supplied 3284–3286
extracts contain only the selected general-run fields. The initial campaign retains individual runs for diagnostics, with the same
user-approved three-foil/centered-sieve assumptions for all four. This choice
does not imply different target configurations and does not block initial tests.
Centered SHMS sieve remains the approved assumption; the independent target
coordinates remain user-supplied information. Four initial campaign rows and
matching shared optics metadata entries have been added. The reports do not
change PID thresholds, geometry constants, or solver settings. An appropriate
SHMS seed matrix is still needed before fitting; none has been invented.

Once the new probe files have been transferred with their directory layout,
run from the repository on ifarm:

```sh
bash diagnostics/run_replay_probe.sh > replay_probe.txt 2>&1
```

The report includes both schema status and threshold counts. A successful
schema check is not successful optics validation. The first 50,000 events can
include startup conditions; this is a bounded smoke test, not a full-run
efficiency estimate. The optional `REPLAY_INPUT_DIR` and `PROBE_MAX_EVENTS`
environment variables are for changed locations/sample sizes and local tests;
the default command requires no extra input.

### Validation record

- `python3 tests/test_replay_probe.py`: passed. Dual-arm isolation, exact cut
  counts, missing/nonfinite/empty inputs, invalid names, wrapper and unchanged
  input checksums.
- `python3 tests/test_spectrometer_switch.py`: shared C++/Python profiles,
  independent −10/0/+10 cm target on both arms, HMS truth regression against
  every TFit column from baseline `057b908`, SHMS geometry closure including
  displaced beam, GMM-mask intersection through TFit, column 10 and δ=16%,
  toy SVD event counts/output names, explicit shifted-sieve rejection, and
  the residual-plot boundary/interval regression.
- `python3 tests/test_switch_python_workflow.py`: passed. Both X/Y cleanup
  scripts with selected-arm raw-branch fallbacks, final physical column,
  independent target selection, SHMS inverse-vertex closure, and stability
  plots that include all three foils for both arms.
- `python3 tests/check_workflow_syntax.py`: all 17 ROOT macros loaded successfully in fresh ROOT interpreters. Python
  entry-point imports/`--help` and shell syntax checks passed. Test-only Python
  dependencies were installed in an isolated temporary environment. Shell/fixed-theta
  dry runs verified automatic routing and 24 commands (three foils × eight
  intervals) for the same target on each arm.

The four new SHMS campaign rows were also checked against both the Python
and C++ metadata readers: matching run IDs, +8.915° angle, −10/0/+10 cm
foils, centered flag, delta boundaries and exact replay paths. All pre-existing
shared metadata bytes remain unchanged.

These are software checks using synthetic events. They do not establish
SHMS optics accuracy, PID purity, peak/band choices, or replay conventions.
No real SHMS replay was read and no usable SHMS calibration matrix was made.
The one-term matrices in the integration test are temporary software fixtures.
The automated ridge/angle/relabel stages have been compiled, but their physics
behavior still needs the minimal real-replay pass and review with Holly/Mark.

## Reference revisions

All reference clones are read-only for this task.

| Repository | Inspected revision |
|---|---|
| `hszumila/HMS_OPTICS_HIGH_MOM` | `58d62dc0bf2f8b1604de4f9f38b31337c77e0396` |
| `hszumila/deltaOpt` | `1e3ecfce5071b4a13587dd75706c92f1091f84d3` |
| `hszumila/SHMS_optics` | `c7d70e88b0cf2d1d93359d941168ab8f35ab916e` |
| `hszumila/cafe_optics` | `4bdd2f949055f0a41f7cac81e8825dc75962f306` |
| `hszumila/shms_optics_a1n` | `cf7f20a8f838493051c4bd5ebc7c21f983a764fd` |
| `hszumila/optics-optimization-shms` | `f772b716d6278435ec8da4ce8a9aafd3d9d6d577` |
