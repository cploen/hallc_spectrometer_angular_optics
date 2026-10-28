# Legacy residual diagnostics with a new HMS matrix

From the ifarm checkout on `codex/joint-surface-ifarm`, using the Python
environment already used for elastic net and beam search and a loaded ROOT
environment:

```bash
git pull --ff-only
campaign=HMS_6p667GeV_joint_min5
matrix="$campaign/06e_beam_search/beam10/matrices/beam.dat"

# Apply beam.dat to full replay events and inspect replay closure first.
./run_sieve.sh "$campaign" beam10_sieve --matrix "$matrix" --offsets 0 0 0 --check
./run_sieve.sh "$campaign" beam10_sieve --matrix "$matrix" --offsets 0 0 0

# Apply the same matrix to the saved training TFit events, then run all three
# legacy residual macros, including their foil/delta/sieve pages.
./run_legacy_diagnostics.sh "$campaign" beam10_legacy \
  --matrix "$matrix" --offsets 0 0 0 \
  --fit-dir "$campaign/06c_core_ntuple/equal15/fit/root"
```

The fitted `00000` coefficients are applied as part of `beam.dat`. Keep the
external additions at zero; do not copy those coefficients into `--offsets`.
The explicit offsets, when needed for other matrices, have units
`[xptar mrad, ytar cm, yptar mrad]`. Neither tool estimates offsets from data.
Choose a new output name to repeat a run; existing outputs are never replaced.

The sieve command uses the campaign replay files, electron/foil cuts and default
angle 12.490 degrees. Use `--angle` for another configured angle. It writes
`07_diagnostics/beam10_sieve/plots/`, `tsv/closure.tsv` and `manifest.json`.
Read any replay-closure warnings before interpreting the reconstruction.

The legacy command processes every nonempty configured rungroup in the specified
TFit directory, with its own metadata angle. These plots describe the saved
`equal15/fit` training sample, not an independent validation sample or the broader
full-replay selection used by the sieve command. It does no selection, fitting
or allocation. To inspect an already-built holdout/surplus export, supply that
export's `root/` directory explicitly and choose a distinct output name.
The original GMM-style `06a_fit_ntuple/root` directory can also be supplied.

Products under `07_diagnostics/beam10_legacy/`:

- `root/`: new `TFit` trees with the afterburner predictions.
- `residuals/<rungroup>/`: seven legacy PDFs and three histogram ROOT files
  per rungroup, for xptar, yptar and ytar.
- `logs/`: ROOT output, including fit warnings from sparse histograms.
- `manifest.json`, `matrix.dat`, `optics.dat` and the rungroup table: settings,
  exact input hashes, event counts and iteration-limit counts.

All residuals are reconstructed minus truth: angles in mrad, ytar in cm. The
previous +/-0.035 cm shifts in two ytar histograms have been removed.
`--dry-run` validates input paths/branches and prints the export/ROOT commands
without creating files. `--export-only` writes the trees without invoking ROOT;
they can then be passed directly to `diagnostics/fit/plot_*_residuals.C`.
`PYTHON` and `ROOT` can select the Python and ROOT executables.

The adapter uses the saved TFit `reactycalc` (= reaction y), `ztar` (= reaction z)
and focal-plane coordinates. It calls the same iterative `xtar` reconstruction
as the sieve tool, then replaces `xtar`, `xptar`, `ytar`, `yptar`, `xs` and `ys`.
Positions remain in cm and target slopes in rad in the exported tree. Original
predictions are retained in `replay_*` branches; event IDs, labels, all truth
values and other branches are retained. Old HMS ntuples with the uncorrected
`xptarT` definition are rejected with a rebuild instruction.

This is a fixed-vertex afterburner. `ztar` and `delta` remain the saved replay
values, including in the legacy auxiliary vertex histograms; they are not new
vertex or momentum reconstructions. No HCANA replay is required. The existing
`diagnostics/run_diagnostics.sh` still reports diagnostics from its supplied
original trees; use the new runner when applying another matrix first.
