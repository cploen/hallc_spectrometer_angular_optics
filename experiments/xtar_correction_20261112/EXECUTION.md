# Remote execution

Started 12 November 2026, approximately 11:24 EST, ifarm2402.

- Source commit: `10f03bb1ecf116a53be8baf622083ee8deb517d4`
  on `codex/xtar-correction-20261112` in
  `cploen/hallc_spectrometer_angular_optics`.
- Source directory:
  `/home/cploen/nps_link/xtar_correction_2018_20261112`.
- Inputs:
  `/home/cploen/nps_link/ML_HMS_OPTICS_JOINT_SURFACE/HMS_LowMom_LowDelta_2018_Accepted`.
- Results directory: `results`; log: `run.log` within source directory.
- Source dependencies copied into `dependencies/`; source hashes recorded in
  `source_sha256.txt`, analysis checkout commit in `checkout_commit.txt`.
- Remote Python versions: numpy 2.0.2, scipy 1.13.1, scikit-learn 1.6.1.
- Seven targeted tests passed locally and on ifarm. Plot layout checked with
  synthetic data under local `layout_check/`; those plots are not real results.
- Initial preflight PID 1280242 stopped before loading data because the older
  campaign checkout's loader lacks the `xtar_coordinate` argument. No fit ran.
  Its output/log were preserved under `preflight_loader_failure/`.
- Actual source checkout is `/u/home/cploen/xtar_refit_2018_20261111`, commit
  `b184e5ad6c0d96c91135948b429aabae2a9c19a8`, whose loader checks geometric
  xtarT. The new correction runner remains the only experiment invoked.
- Restart at approximately 11:25 EST: PID 1283380. Completion not yet verified.

The README contains the foreground command. The actual launch used that same
command under `nohup`, with output redirected to `run.log`, and four threads.
