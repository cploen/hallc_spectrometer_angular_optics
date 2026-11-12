# xtar tolerance scan

Rescore the saved candidate set, with unchanged data, run-group folds,
geometry, constants, EN/ridge penalties, six-term cap and 0.2 cm correction cap.
Compare cell/overall downstream SD limits of 2%/2%, 5%/2%, and 10%/5%.
The signed-bias guard remains 2% of the corresponding baseline width.
Populated cells have N>=100; core, shoulder and union are guarded separately.

```bash
E=/home/cploen/nps_link/xtar_correction_2018_20261112
cd "$E"
python3 -B -u xtar_tolerance_scan.py --experiment "$E" \
  --checkout /u/home/cploen/xtar_refit_2018_20261111 \
  --out "$E/tolerance_scan_5_2_and_10_5"
```

The original 2%/2% eligibility and validation scores must reproduce before any
new reserved evaluation. Shared fold solves are cached and use one BLAS thread;
all unpenalized/regularized numerical conventions remain unchanged.

For each tolerance, report the strongest eligible EN and beam candidate and
the original simplicity-selected candidate separately. Freeze each strongest
model before reserved evaluation. Report xtar width gain, all downstream overall
and worst-cell SD changes, separate biases, and whether the reserved guards pass.
No candidate is installed or adopted automatically. These previously inspected
reserved events remain development evidence.

This scan evaluates the original saved candidate set; it does not claim an
exhaustive search under the relaxed constraints. New relaxed-gate beam searches
could explore additional bases if needed.
