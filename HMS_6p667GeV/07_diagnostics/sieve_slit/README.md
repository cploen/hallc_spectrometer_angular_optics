# Sieve slit scattering study: target-y distributions

This diagnostic compares saved replay `ytar` (cm) in dense sieve cores and the
out-of-core population. It is parallel to `geometry` and `compare` in
`07_diagnostics`; it does not refit a matrix or change selection cuts.

## Populations

- **Densest 50% of cores:** `quality == 2`, ranked by decreasing saved
  `core_score` separately in each run group, foil, delta slice, and X/Y hole.
  Keep `ceil(N/2)` events; break score ties by ascending original entry number.
  Thus odd-size strata contribute one extra event. This is half the accepted
  core events, not half the holes or a density-score threshold of 0.5.
- **Shoulders / between cores:** `quality == 1`, outside the accepted core of
  a supported hole. This is an operational noncore category; there is no
  additional geometric separation of shoulders from inter-hole valleys.
- **All out-of-core events:** `quality != 2`, including unsupported holes
  (`quality == 0`). Unsupported holes are tabulated separately and are not
  silently assigned to shoulders.
- The remaining half of accepted cores is tabulated as `outer_core`; it is
  not counted as out of core. Fit, surplus, and holdout events are all included.

Each run group and foil has three separate PNGs, with identical x limits and
100 common bins. All labeled holes and delta slices are pooled within that
foil. Counts are unweighted. `summary.tsv` includes N, mean, population standard
deviation, median, and 5th/95th percentiles for all five populations.
`histograms.tsv` stores every bin, including zero counts, for independent replotting.
No ytar range clipping is applied by this diagnostic.

## Findings for HMS 6.667 GeV

Across all 12 run-group/foil combinations, the supported noncore distribution
is broader than the densest core half: its standard deviation is about 11–21%
larger. For rg01 (12.490 degrees, foil 0 cm), the dense half contains 91,069
events with mean −0.174 cm and SD 0.122 cm; the shoulders contain 50,023 events
with mean −0.179 cm and SD 0.147 cm. Another 494 events are from unsupported
holes. The shoulder means remain within about 0.02 cm of the dense-core means
throughout this campaign.

These widths describe a selected, pooled population, not an intrinsic resolution
or a measured scattering fraction. Inputs contain unambiguous joint X/Y labels
and inherit the candidate producer's cuts, including target-y selection. Missing
or ambiguous candidates and events rejected upstream are unavailable here.
Those cuts can truncate tails. Different hole/delta mixtures can also affect
pooled widths. A broader noncore distribution is consistent with contamination,
but does not establish a sieve-scattering origin.

## Reproduce

From the repository root, with Python 3.11+ and NumPy, Matplotlib, and uproot:

```bash
./run_sieve_slit.sh HMS_6p667GeV --check
./run_sieve_slit.sh HMS_6p667GeV
# Optional: one run group, or an alternate configuration
./run_sieve_slit.sh HMS_6p667GeV --rungroup rg01_theta12p490_foil0
./run_sieve_slit.sh HMS_6p667GeV --config HMS_6p667GeV/config/sieve_slit.json
python3 -m unittest discover -s tests -p test_sieve_slit.py
```

The wrapper honors `PYTHON`. The actual source is
`diagnostics/validation/sieve_slit/plot_ytar.py`. Default configuration is
`<campaign>/config/sieve_slit.json`; every configured relative path is resolved
against the campaign directory, independent of the invoking working directory.
The rungroup table supplies names and optics IDs; `core_file` supplies the input
path template. All requested inputs are checked before writing plots.

For another campaign, create `HMS_<name>` or `SHMS_<name>` with a configuration
containing `rungroup_table`, `core_file`, `output`, and `bins`, then invoke the
same wrapper with that folder. Spectrometer identification uses the shared
`spectrometer_config.from_campaign`; no HMS branch names or geometry constants
are embedded in the diagnostic. Inputs must use the existing `CoreSample`
schema and quality definitions. HMS was run here; SHMS data were not available
for an end-to-end run.

The input files are the seven `05c_core_sample/min10/root/CoreSample_*.root`
files. Preserve these files for exact reproduction; the diagnostic does not
regenerate the frozen core selection. `manifest.json` records SHA-256 hashes
of every input, configuration, rungroup table, and plotting source, plus package
versions. Consult `CORE_SAMPLE.md` for upstream candidate/core production.
The committed PNGs and TSVs are compact analysis outputs; event ROOT files are
not duplicated in this study directory. Re-running overwrites outputs with the
same names. A partial `--rungroup` run updates the manifest/tables only for that
selection; use a separate configured output directory to preserve a full run.
