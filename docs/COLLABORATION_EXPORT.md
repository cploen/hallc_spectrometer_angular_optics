# Collaboration slide inputs

Run from the production project root on ifarm after the corrected standard
diagnostics have finished for all three campaigns:

```bash
python3 export_collaboration_update.py --make-sieve
git add -f collaboration_exports/nps_angular_update_20261103.zip
git commit -m "Export NPS angular fit results"
git push origin HEAD:codex/joint-surface-ifarm
```

The export uses `beam20_<tag>_standard_2024flags`, not the earlier zero-offset
comparison. It checks the pinned 2024 external offsets, zero external offsets
on modern matrices, matrix checksums, identical before/after replay inputs,
selection settings and foil event counts.

`--make-sieve` generates only missing representative sieve views, using the
existing matrices. It does not refit. The representative angles are 16.500°,
12.385° and 12.490° for 5.878, 6.117 and 6.667 GeV. Each view retains every
available foil at that exact angle. Nearby angle settings remain separate.
Both matrices use the same PID, delta and foil windows. Full replay sieve
histograms include events outside the training and protected samples.

The ZIP contains numeric heatmap tables, fold results, conditioning and term
tables, matrices, sieve histogram arrays, selected existing plots, and source
manifests. Its manifest records a SHA-256 checksum for each source file. It
contains no ROOT event files. These sources allow larger labels and a common
color scale when preparing the slides. `core_shoulder` means protected cores
plus **all** supported shoulders, with their original event counts.

If an existing sieve view has different settings or a different matrix, the
export stops instead of silently replacing it. An existing ZIP is also kept:
choose a new path with `--output` for a later export.
