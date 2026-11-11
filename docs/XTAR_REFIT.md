# Conditional HMS x-target refit

This explicit experiment uses the existing frozen 2018 selections and full
geometric `xtarT` training coordinate. It changes only x-target-dependent
angular coefficients through a requested polynomial order. All non-x-target
terms, higher-order x-target terms, delta coefficients and matrix constants
are frozen. The zero-offset baseline must have zero matrix constants.

The first experiment allows the 56 existing x-target monomials through order
four, then selects smaller supports with training-only elastic net and beam
width 30. Coefficients in this free block are fitted in full; excluded terms
in the block are zero. The regularization reference is zero. Higher-order
simulated rows stay unchanged. No extra basis terms or fit intercept are added.

From an isolated source checkout on ifarm:

```bash
A=/u/group/nps/cploen/ML_HMS_OPTICS_JOINT_SURFACE
C="$A/HMS_LowMom_LowDelta_2018_Accepted"
B="$A/HMS_LowMom_2018_SievePreview_20261110/config/beam_zero_offsets.dat"
F="$C/06e_beam_search/beam30_inner20/tsv/folds.tsv"

bash run_xtar_refit.sh "$C" inner20 xtarT_order4_beam30 \
  --baseline "$B" --folds-from "$F" --max-order 4 --beam 30 \
  --steps 40 --threads 8 --check

bash run_xtar_refit.sh "$C" inner20 xtarT_order4_beam30 \
  --baseline "$B" --folds-from "$F" --max-order 4 --beam 30 \
  --steps 40 --threads 8
```

The baseline is the accepted 2018 Beam matrix with its complete constant row
removed, as used for the zero-offset replay experiment. Its expected SHA256 is
`f715879845f777f2da2aa62b7e66756907929bdf4b7b40809ac55370156f8801`.
The `--check` output prints the hash for comparison. Existing selection IDs,
80% core / 20% inner-shoulder mixture, metadata and saved training folds are
checked rather than regenerated. Missing or inconsistent geometric coordinates
stop the experiment. The old angular-fit entry points retain their defaults.

Outputs are immutable under `06g_xtar_refit/xtarT_order4_beam30/`:

- `seed.dat` and `manifest.json`: baseline, coordinate, free block, exact input
  hashes, configs, rank/conditioning and status.
- `matrices/{enet,full,start,beam}.dat`: complete transport matrices. Exports
  check frozen rows and every delta coefficient; no matrix is installed.
- `tsv/enet_*`, `tsv/search.tsv`, `tsv/conditioning.tsv`, `tsv/terms.tsv`:
  selection scores, sparsity, rank and coefficient records.
- `plots/cv.png`, `plots/terms.png`, `plots/development_cv_*`: training-only
  diagnostics. Predictions use fold fits, but adaptive basis selection and
  previously fitted frozen terms make these conditional development scores.

After completing the frozen fit, evaluate without fitting or changing terms:

```bash
bash run_xtar_refit.sh "$C" inner20 xtarT_order4_beam30 --evaluate
```

`evaluation/` contains reserved-ID residual heatmaps and geometric sieve
projections. Training/holdout ID overlap is rejected. Existing inspected
holdouts are development evidence, not a fresh independent test. The full
56-column reference can be rank-truncated; inspect the manifest. Beam supports
must have full numerical rank in every training fold.

Compare the candidate against `baseline_geometry`: both use `xtarT`.
`baseline_replay_coordinate` separately evaluates the old matrix with saved
replay `xtar`, isolating the coordinate substitution. All sieve projections
use the geometric `xtarT` anchor and are labeled as such; they are not replays.
Heatmaps include bias in RMS, show counts, and preserve sparse cells.
Sieve views cover [-8,10), [10,12), [12,14), and [-15,-8) percent when populated.

The user chose to skip afterburner closure as a prerequisite. These results do
not certify closure or deployable reconstruction. Review the matrix and plots
before a matched baseline/candidate HCANA replay. The separate coordinate
handoff records unresolved replay/HCANA convention questions for review.
