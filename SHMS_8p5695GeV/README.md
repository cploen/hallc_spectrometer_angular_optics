# SHMS 8.5695 GeV/c campaign

The campaign name `SHMS_8p5695GeV` follows the existing spectrometer-and-momentum convention. Its `SHMS_` prefix selects the spectrometer. Runs 3283–3286 form one run group with shared reported settings and the assumed three-foil target, following the existing HMS campaign organization. This is a preparation and validation campaign, not a finished calibration.

## Established settings and approved assumptions

- All four replay reports agree: SHMS 8.5695 GeV/c, +8.915°, electron mass hypothesis; beam energy 10.6 GeV; carbon target mass setting 12.0107 amu.
- User instruction: assume the supplied foil coordinates **−10, 0, +10 cm for all four runs**. Do not require reconfirmation before initial processing. This target is independent of arm and is also available for HMS.
- User instruction: assume **centered SHMS sieve** (`SieveFlag=1`); position-option work remains deferred.
- User supplied an ifarm file listing confirming all four ROOT paths exist. Approximate displayed sizes: 3283 842M, 3284 791M, 3285 847M, 3286 709M. Displayed timestamps: July 7, respectively 15:24, 15:20, 15:17, 15:14; the listing did not show a year. Their contents have not been read from this Mac.
- Delta boundaries: **−10, −6, −2, +2, +6, +10, +22 percent** (six slices). These follow Holly's shms_optics_a1n run-group 10301 divisions within our unchanged −10 to +22% acceptance: omit her −12 to −10 slice and cap her +10 to +25 slice at +22. The combined and individual-run DAT entries use seven boundaries. Regenerate angle scans and downstream slice-dependent outputs after this change; old slice indices and tuned angles must not be reused as if the intervals were unchanged. Existing foil-ridge polygons do not depend on these slice boundaries.

## Inputs

`config/rungroups_8p5695_inputs.tsv` contains one group,
`rg01_theta8p915_foilpm10z0`, with optics ID **8569501** and all four run
numbers. Its corresponding entry is in `DATfiles/list_of_optics_run.dat`.
The individual IDs 3283–3286 remain available for diagnostics.

From the repository root on ifarm, prepare the combined ROOT input:

```sh
bash diagnostics/prepare_shms_8p5695_inputs.sh
```

This uses ROOT's `hadd` to combine the four original files into
`SHMS_8p5695GeV/inputs/shms_optics_8p5695_rg01_theta8p915_foilpm10z0.root`.
It checks source readability and refuses to overwrite an existing output.
Allow roughly the combined input size in additional storage. If a merge fails,
its partial output must be inspected and moved aside before retrying.
The table uses a repository-relative path: run workflow commands from the
repository root. Numbered output folders are created by the standard runners.

The next workflow stage is foil/ridge selection:

```sh
bash run_ytar_ridge_all_rungroups.sh SHMS_8p5695GeV
```

Inspect and tune those cuts before continuing through the standard candidate,
GMM, and fit stages in the repository README. This command is a normal workflow
run, not the bounded 50,000-event diagnostic.

`config/ztar_runlist.txt` is the existing standalone ztar validation format with the exact replay path appended. Foil index 0 means −10 cm, index 1 means 0 cm, index 2 means +10 cm. At positive SHMS angle the ridge ordering in ytar is reversed; the geometry-aware ridge assignment preserves these metadata IDs.

The original replay matrix/offset configuration remains unknown. No `config/oldfit.dat` seed is fabricated: matrix fitting waits for an appropriate SHMS seed. Existing PID thresholds are inherited initial settings; detector branches and cut populations must be checked on the replays.

After these files are installed on ifarm, the bounded read-only probe is:

```sh
bash diagnostics/run_replay_probe.sh SHMS_8p5695GeV > replay_probe.txt 2>&1
```

No X11 is needed for the probe. Reported tests read 50,000 events in each of
the four runs with zero read errors. The probe still flags 5–10 invalid or
nonscalar `P.extcor.ysieve` values per run; additional extreme finite values
occur before selection. For run 3283, none passed the tested PID and nominal
delta selection with the diagnostic `!(abs(P.extcor.ysieve)<1000)` condition.
That selected sample contains 16,686 events, visible sieve holes, and three
vertex peaks near −10, 0, +10 cm. A central-foil, narrower-delta view retains
2,427 events with separated holes. These are user-supplied ifarm results and
screenshots, not validation of the new fit geometry or the other runs' selected
distributions.

## Sharing

Share this campaign with the matching repository branch
`hms-shms-campaign-switch`, including its shared macros and DAT metadata.
The campaign folder alone is not a standalone program. ROOT inputs and generated
plots remain outside Git; collaborators with access to the original replay
directory can reproduce the combined input with the preparation command above.
No starting matrix is supplied yet: an appropriate SHMS `config/oldfit.dat`
is still needed before fitting.

See [the implementation record](../docs/HMS_SHMS_CHANGES.md) and [reference review](../docs/HMS_SHMS_REVIEW.md).
