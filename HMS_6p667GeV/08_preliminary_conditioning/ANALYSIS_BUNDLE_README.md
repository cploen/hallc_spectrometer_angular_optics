# HMS 6.667 GeV preliminary angular-optics analysis bundle

This directory contains the data products, coefficient arrays, tables, plots,
summaries, and collaboration slide deck produced for the preliminary angular
optics conditioning and stability study.

The analysis makes no recommendation about committing a production matrix.

## Source inputs already tracked on this branch

- `HMS_6p667GeV/config/oldfit.dat`
- `HMS_6p667GeV/config/rungroups_6p667_inputs.tsv`
- `HMS_6p667GeV/config/pre_gmm_veto.tsv`
- The seven fit ntuples under `HMS_6p667GeV/06a_fit_ntuple/root/`
- `HMS_6p667GeV/06b_svd_fit/matrices/nps_hms_newfit_6p667_gmm_clean.dat`
- The fit logs and build summaries under `HMS_6p667GeV/06b_svd_fit/`

## Analysis code

- `diagnostics/conditioning/preliminary_angular_conditioning.py`
- `diagnostics/conditioning/angular_term_ladder.py`
- `diagnostics/conditioning/angular_solver_comparison.py`
- `diagnostics/conditioning/angular_per_foil_residual_ladder.py`
- `fit_opt_matrix_gmm.C`, including the strict matrix-parser check used during
  this study

## Results in this directory

- `REPRODUCTION.md`: chronological runbook, exact commands, sample definitions,
  solver definitions, validation landmarks, and retained ambiguities
- `reproduction_check/`: preserved uncorrected and strict-parser local rerun
  matrices and their QA files
- `angular_singular_values.tsv` and
  `leave_one_setting_out_conditioning.tsv`: initial conditioning diagnostics
- `term_ladder/`: current unscaled-XTX term ladder, solve-specific truncation
  margins, residuals, prediction steps, and saved coefficient arrays
- `solver_comparison/`: unscaled-XTX, direct-X, rescaled-direct-X, and
  rescaled-XTX comparisons, including training and held-out results
- `per_foil_residual_ladder/`: training and held-out residual ladders for each
  foil, target variable, and the two requested solve routes
- `COLLABORATION_SUMMARY.md`, `PRELIMINARY_FINDINGS.md`, and text summaries:
  written records of the preliminary findings and limitations
- `HMS_6p667_angular_stability_preliminary.pptx`: collaboration slide deck with
  speaker notes

The TSV files are the primary human-readable numerical record.  The NPZ files
retain fitted coefficient ladders so the solver results can be inspected
without rerunning the fits.  PDF and PNG files contain the same plotted results
in vector and raster formats.

Excluded from version control are Python caches, ROOT/ACLiC compilation
by-products, duplicate slide-render QA images, and temporary PowerPoint lock
files.  These are reproducible or transient and are not analysis inputs or
results.
