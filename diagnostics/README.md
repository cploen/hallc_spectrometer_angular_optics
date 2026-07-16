# Fit diagnostics

Run from anywhere:

```bash
/path/to/hallc_spectrometer_angular_optics/diagnostics/run_diagnostics.sh CAMPAIGN
```

`CAMPAIGN` may be a directory path or a campaign directory name under the
repository root.

The runner locates the campaign configuration and fit trees, runs
`analyze_fit_sample_balance.C`, and produces xptar, yptar, and ytar residual
diagnostics for each configured fit tree. Products are written under:

```text
<campaign>/07_diagnostics/fit/
```

Use `--dry-run` to print commands without running ROOT.
