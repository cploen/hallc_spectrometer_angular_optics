# Single-pass X split preparation

The X fixed-angle runner previously read the complete replay once per foil
and delta slice, with all branches active. For three foils and six slices,
split preparation alone required 18 full reads before accounting for the
separate candidate-generation macros.

run_x_multifoils_fixedtheta.py now uses compute_splits to read once with only
the selected arm's Cherenkov, cal.etottracknorm, gtr.dp, gtr.y and dc.xp_fp
branches enabled. It collects each foil/slice sample and applies the existing
quantile function and minimum 1,000-event requirement. Independent foil
membership is preserved, including events in overlapping supplied cuts.

The PID thresholds, SHMS acceptance, half-open delta intervals, 5th/95th
percentile midpoint, saved reference angles, low/high gates, and downstream
macro commands are unchanged. No sampling or event limit is introduced.
Missing branches or failed entry reads stop with an explicit error.
Files are closed even if preparation fails.

Verbose progress now describes all foils/slices together, reports one pass
and five branches, and counts selected assignments (which can include an
event more than once if cuts overlap). All split samples are held in memory
together, increasing peak memory compared with one sample at a time.
The subsequent C++ macro runs still perform their own reads.
The existing --dry-run still computes splits and then prints macro commands.

tests/test_x_split_single_pass.py compares exact counts and split values with
the frozen pre-change reader on synthetic ROOT input for both arms, including
delta/PID boundaries and overlapping/empty foil cuts. It checks that input
files are unchanged. Real ifarm timing remains to be measured.

Run from the repository root:

```sh
python3 -u run_x_multifoils_fixedtheta.py SHMS_8p5695GeV \
  rg01_theta8p915_foilpm10z0 rg01_theta8p915_foilpm10z0 1 --verbose
```
