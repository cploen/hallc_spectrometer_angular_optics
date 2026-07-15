# Hall C Spectrometer Optics Calibration

Semi-automated preparation of HMS and SHMS sieve data for optics-matrix fitting.

## 1. Software requirements

The workflow requires:

- HCANA
- ROOT with PyROOT
- Python 3
- NumPy
- pandas
- Matplotlib
- scikit-learn
- uproot

Check the environment:

```bash
command -v hcana
command -v root
python3 -c 'import ROOT, numpy, pandas, matplotlib, sklearn, uproot; print("Dependencies OK")'
```

## 2. Minimum campaign setup

Create a campaign directory containing one rungroup table and one starting optics matrix:

```text
HMS_<campaign>/
└── config/
    ├── rungroups_<campaign>_inputs.tsv
    └── oldfit.dat
```

Example:

```bash
CAMPAIGN=HMS_6p667GeV
mkdir -p "$CAMPAIGN/config"
cp /path/to/starting_matrix.dat "$CAMPAIGN/config/oldfit.dat"
```

The rungroup table must be tab-separated and contain:

```text
rungroup	optics_id	hms_angle_deg	foils	nruns	runs	rootfile
```

Example:

```text
rg01_theta12p490_foil0	666701	12.490	0	1	1544	/path/to/rg01.root
rg02_theta15p195_foil0	666702	15.195	0	1	5162	/path/to/rg02.root
rg03_theta12p490_foilpm8	666703	12.490	-8,8	1	1540	/path/to/rg03.root
```

Each `optics_id` must have a corresponding entry in:

```text
DATfiles/list_of_optics_run.dat
```

Each `rootfile` must exist before beginning the workflow.

The runner scripts create the numbered campaign output directories. Do not create them manually.

## 3. Workflow overview

1. Build and inspect the `ytar` ridge cuts.
2. Generate YFP/YPFP and XFP/XPFP band candidates.
3. Relabel the bands as physical `yscol` and `xscol`.
4. Build the event-level candidate trees.
5. Apply GMM cleanup.
6. Build the fit ntuple and optics matrix.
7. Replay and validate the result.

## 4. Projection tuning

Tune each spectrometer setting using its center-foil rungroup first. This
usually minimizes the adjustments needed for the other foils.

### Y projection

Choose the center-foil rungroup for one spectrometer setting and run the
initial Y angle scan:

```bash
CAMPAIGN=HMS_6p667GeV
REFERENCE_RUNGROUP=rg01_theta12p490_foil0
REFERENCE_FOIL=0

./run_y_angle_scan.sh \
  "$CAMPAIGN" \
  "$REFERENCE_RUNGROUP" \
  "$REFERENCE_FOIL"
```

Inspect the diagnostic PDFs. Close acceptable plots and leave only concerning
plots open in Evince.

Generate HCANA rerun commands for the open diagnostics:

```bash
python3 make_angle_rerun_commands.py "$CAMPAIGN"
```

Run the generated commands, inspect the new PDFs, and tune iteratively until
the center-foil projections are acceptable.

Apply the accepted center-foil angles to another rungroup at the same
spectrometer setting:

```bash
TARGET_RUNGROUP=rg03_theta12p490_foilpm8

python3 run_y_multifoils_fixedtheta.py \
  "$CAMPAIGN" \
  "$TARGET_RUNGROUP" \
  "$REFERENCE_RUNGROUP" \
  "$REFERENCE_FOIL"
```

Here:

- `REFERENCE_RUNGROUP` is the center-foil rungroup whose angles were tuned.
- `TARGET_RUNGROUP` is another foil rungroup at the same spectrometer setting.

Inspect the resulting PDFs and repeat for the other target rungroups at that
setting.

### X projection

Repeat the same center-foil-first procedure for X:

1. Run the X angle scan for the center-foil rungroup.
2. Close acceptable diagnostic PDFs and leave concerning PDFs open.
3. Generate and run the HCANA rerun commands.
4. Tune iteratively until the center-foil projections are acceptable.
5. Apply the accepted angles to each target rungroup:

```bash
python3 run_x_multifoils_fixedtheta.py \
  "$CAMPAIGN" \
  "$TARGET_RUNGROUP" \
  "$REFERENCE_RUNGROUP" \
  "$REFERENCE_FOIL"
```

Complete both Y and X projection validation before proceeding to relabeling.

## 5. Relabel the projection bands

After the Y and X projections are accepted, assign the physical sieve-column
labels to the selected bands.

### Y relabeling

Run Y relabeling for the full campaign:

```bash
./run_relabel_y_all.sh "$CAMPAIGN"
```

Inspect the PDFs in:

```text
CAMPAIGN/03a_relabel_y/plots/
```

Verify that the `yscol` labels follow the physical sieve columns consistently
across all foils and delta slices.

### X relabeling

Run X relabeling for the full campaign:

```bash
./run_relabel_x_all.sh "$CAMPAIGN"
```

Inspect the PDFs in:

```text
CAMPAIGN/03b_relabel_x/plots/
```

Verify that the `xscol` labels follow the physical sieve columns consistently
across both dynamic XFP/XPFP regions.

Do not proceed until both Y and X relabeling are acceptable.

The next step is to build the event-level candidate trees.

## 6. Build the event-level candidate trees

Before building the candidate trees, optionally veto any clearly bad relabeled
components.

### Optional pre-GMM veto file

Create the veto file only after inspecting the relabeling diagnostics:

```text
CAMPAIGN/config/pre_gmm_veto.tsv
```

Each line contains:

```text
rungroup axis exact_cut_name
```

Example:

```text
# rungroup axis exact_cut_name
rg03_theta12p490_foilpm8 y hYpFpYFp_cut_yscol_4_nfoil_0_ndel_2_part00_src_03
rg06_theta15p195_foilpm8 x hXpFpXFp_cut_xscol_7_nfoil_1_ndel_1_part00_zone_low_src_04
```

Use the exact cut name written by the relabeling step. The axis must be `x` or
`y`.

Vetoed events are excluded before GMM cleanup. They do not enter the candidate
tree, consume any part of the GMM removal budget, or reach the final fit.

If no vetoes are needed, do not create the file.

### Build the Y candidate trees

```bash
./run_make_yscol_candidates_all.sh "$CAMPAIGN"
```

Outputs are written under:

```text
CAMPAIGN/04a_candidate_trees_y/
```

Check the logs for missing veto names, overlapping assignments, or failed
rungroups.

### Build the X candidate trees

```bash
./run_make_xscol_candidates_all.sh "$CAMPAIGN"
```

Outputs are written under:

```text
CAMPAIGN/04b_candidate_trees_x/
```

Check the logs for missing veto names, overlapping assignments, or failed
rungroups.

The next step is GMM cleanup.

## 7. GMM cleanup

Run the GMM cleanup after both candidate-tree stages complete successfully.

Default settings:

```text
MIN_EVENTS=30
MAX_COMPONENTS=9
KEEP_FRAC=0.90
```

Run Y cleanup:

```bash
./run_gmm_cleanup_y_all.sh "$CAMPAIGN"
```

Run X cleanup:

```bash
./run_gmm_cleanup_x_all.sh "$CAMPAIGN"
```

The Y and X runners use separate inputs and outputs and may be run
simultaneously in separate terminals.

Outputs are written under:

```text
CAMPAIGN/05a_gmm_cleanup_y/
CAMPAIGN/05b_gmm_cleanup_x/
```

Inspect the diagnostic PDFs before proceeding. Check that the retained events
follow the intended sieve-hole population and that valid structure has not
been removed.

### Change a GMM setting

Override a default by setting its environment variable before the runner.

Example: retain 95% of the X events:

```bash
KEEP_FRAC=0.95 ./run_gmm_cleanup_x_all.sh "$CAMPAIGN"
```

Example: require at least 40 events for Y:

```bash
MIN_EVENTS=40 ./run_gmm_cleanup_y_all.sh "$CAMPAIGN"
```

Example: change multiple settings:

```bash
MIN_EVENTS=30 MAX_COMPONENTS=12 KEEP_FRAC=0.95 \
  ./run_gmm_cleanup_x_all.sh "$CAMPAIGN"
```

An override applies only to that command. Rerunning an axis replaces its
previous GMM outputs without changing the other axis.

The next step is to build the fit ntuple and optics matrix.

### Rerun one foil and delta slice

To rerun only one foil and one delta slice, call the underlying cleanup script
directly.

For X:

```bash
python3 gmm_cleanup_xscol_candidates_v2.py \
  --campaign "$CAMPAIGN" \
  --rungroup RUNGROUP \
  --foil FOIL_INDEX \
  --ndel DELTA_INDEX \
  --min-events 30 \
  --max-components 9 \
  --keep-frac 0.90
```

For Y:

```bash
python3 gmm_cleanup_yscol_candidates_v2.py \
  --campaign "$CAMPAIGN" \
  --rungroup RUNGROUP \
  --foil FOIL_INDEX \
  --ndel DELTA_INDEX \
  --min-events 30 \
  --max-components 9 \
  --keep-frac 0.90
```

Use a lower `--keep-frac` for more restrictive cleanup.

Example: rerun X for rungroup `rg05_theta12p495_foilpm3`, foil `1`,
and the −5% to 0% delta slice while retaining 90% of events:

```bash
python3 gmm_cleanup_xscol_candidates_v2.py \
  --campaign HMS_6p667GeV \
  --rungroup rg05_theta12p495_foilpm3 \
  --foil 1 \
  --ndel 2 \
  --min-events 30 \
  --max-components 9 \
  --keep-frac 0.90
```

This replaces the GMM outputs for that foil and delta slice only.

## 8. Build the fit ntuples and optics matrix

After the X and Y GMM cleanup is accepted, build the fit ntuples and run the
SVD optics fit.

Choose a descriptive output tag:

```bash
TAG=6p667_gmm_clean
```

First run in dry-run mode:

```bash
./run_build_gmm_fit_and_matrix.sh \
  "$CAMPAIGN" \
  "$TAG" \
  --dry-run
```

Check the printed inputs, output paths, and commands. Then run the full step:

```bash
./run_build_gmm_fit_and_matrix.sh \
  "$CAMPAIGN" \
  "$TAG"
```

The runner:

1. Builds and validates one fit ntuple for each rungroup.
2. Combines the accepted X and Y GMM selections.
3. Runs the SVD optics fit.
4. Writes the fitted matrix and QA outputs.

Outputs are written under:

```text
CAMPAIGN/06a_fit_ntuple/
CAMPAIGN/06b_svd_fit/
```

The fitted optics matrix is written as:

```text
CAMPAIGN/06b_svd_fit/matrices/nps_hms_newfit_<TAG>.dat
```

### Optional environment overrides

Default values include:

```text
FILEID=-1
NFIT_MAX=200000
NSETTINGS=-1
RUN_SVD=1
OVERWRITE=0
OLD_COEFFS=CAMPAIGN/config/oldfit.dat
CER_CUT=6.0
CAL_CUT=0.65
```

Override a value for one run by placing the environment variable before the
command.

Example: allow existing outputs to be replaced:

```bash
OVERWRITE=1 \
  ./run_build_gmm_fit_and_matrix.sh "$CAMPAIGN" "$TAG"
```

Example: limit the number of fit events:

```bash
NFIT_MAX=100000 \
  ./run_build_gmm_fit_and_matrix.sh "$CAMPAIGN" "$TAG"
```

Example: build the fit ntuples without running the SVD fit:

```bash
RUN_SVD=0 \
  ./run_build_gmm_fit_and_matrix.sh "$CAMPAIGN" "$TAG"
```

Inspect the fit-ntuple checks and SVD QA outputs before using the new matrix in
a replay.
