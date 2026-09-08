# SHMS 8.5695 GeV/c campaign

The campaign name `SHMS_8p5695GeV` follows the existing spectrometer-and-momentum convention. Its `SHMS_` prefix selects the spectrometer; run numbers belong in the run metadata. These are initial inputs for individual-run checks, not a validated calibration or approved set of cuts.

## Established settings and approved assumptions

- All four replay reports agree: SHMS 8.5695 GeV/c, +8.915°, electron mass hypothesis; beam energy 10.6 GeV; carbon target mass setting 12.0107 amu.
- User instruction: assume the supplied foil coordinates **−10, 0, +10 cm for all four runs**. Do not require reconfirmation before initial processing. This target is independent of arm and is also available for HMS.
- User instruction: assume **centered SHMS sieve** (`SieveFlag=1`); position-option work remains deferred.
- User supplied an ifarm file listing confirming all four ROOT paths exist. Approximate displayed sizes: 3283 842M, 3284 791M, 3285 847M, 3286 709M. Displayed timestamps: July 7, respectively 15:24, 15:20, 15:17, 15:14; the listing did not show a year. Their contents have not been read from this Mac.
- Initial delta boundaries: −10, −8, −5, 0, 5, 10, 15, 20, 22 percent. These are chosen analysis bins covering the existing SHMS profile window, not values inferred from the reports or tuned using these events.

## Inputs

`config/rungroups_8p5695_inputs.tsv` retains one rungroup per run for the first checks, all with the same assumed three-foil target. This is an analysis organization choice, not a claim that the target differed among runs. Runs need not be combined to begin checking them. The corresponding shared metadata entries use optics IDs 3283–3286 in `DATfiles/list_of_optics_run.dat`.

`config/ztar_runlist.txt` is the existing standalone ztar validation format with the exact replay path appended. Foil index 0 means −10 cm, index 1 means 0 cm, index 2 means +10 cm. At positive SHMS angle the ridge ordering in ytar is reversed; the geometry-aware ridge assignment preserves these metadata IDs.

The original replay matrix/offset configuration remains unknown. No `config/oldfit.dat` seed is fabricated: matrix fitting waits for an appropriate SHMS seed. Existing PID thresholds are inherited initial settings; detector branches and cut populations must be checked on the replays.

After these files are installed on ifarm, the bounded read-only probe is:

```sh
bash diagnostics/run_replay_probe.sh SHMS_8p5695GeV > replay_probe.txt 2>&1
```

No X11 is needed. The probe samples at most 50,000 events per run. File existence and report metadata are established from user terminal output; branch/finite-value/cut-count validation remains pending.

See [the implementation record](../docs/HMS_SHMS_CHANGES.md) and [reference review](../docs/HMS_SHMS_REVIEW.md).
