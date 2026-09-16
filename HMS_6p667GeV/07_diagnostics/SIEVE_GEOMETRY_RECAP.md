# HMS 6.667 — sieve afterburner and unresolved geometry checks

Conversation recap, 15 September 2026. Method:
[SIEVE_AFTERBURNER.md](../../docs/SIEVE_AFTERBURNER.md).
Equation definitions, informal proofs and pinned external references:
[HMS_GEOMETRY_CONSISTENCY.md](../../docs/HMS_GEOMETRY_CONSISTENCY.md).

## What was run and why

The user requested an independent visual view of the reconstruction: GMM versus
beam, at exactly 12.490°, for physical foils 0 and ±8 cm. The tool uses rg01/run
1544 and rg03/run 1540, not the nearby 12.495° groups. It applies the same original
good-electron selection to both matrices: Cherenkov >2, calorimeter E/p >0.65,
−10<delta<10%, and original reaction-z within ±2 cm of the nominal foil.
There are **no GMM or core event-selection cuts** in these plots.

GMM denotes the historical fitted matrix, not GMM-selected events refitted on
beam terms. Beam is `beam10/matrices/beam.dat`. The fixed-vertex afterburner
reevaluates the polynomial and iterates xtar, retaining the original measured
track and reaction vertex. It does not rerun tracking, recompute the reaction
vertex, or refit momentum. Six PNGs share binning, axes and a logarithmic scale.
The side-by-side review suggested modest shape changes, not dramatic universal
sharpening; the unresolved closure check limits interpretation.

## Actual replay-check result reported by the user

The user ran the tool on ifarm and pasted `sieve/tsv/closure.tsv`. These numbers
are a conversation record; the output manifest/tables are not in this local
checkout at the audit date.

| Group | Selected N | xsieve RMS difference (cm) | ysieve RMS difference (cm) |
|---|---:|---:|---:|
| rg01, foil 0 | 265,877 | 0.047633 (pass) | 0.028387 (pass) |
| rg03, foils −8/+8 combined | 98,382 | 0.080956 (fail) | 0.035244 (pass) |

Tolerance is 0.05 cm = 0.5 mm per coordinate. This checks whether applying the
**presumed original matrix** reproduces the original saved sieve coordinates;
it does not compare either new fit with physical truth. One failed coordinate
marks every figure “replay check needs review.” In rg03, 78 reference events
reached the five-iteration limit; whether they explain the failure is untested.

The reference is `config/oldfit.dat` with zero external angular additions.
Source directory naming suggests a zero-offset replay, but **the actual replay
matrix, hmsflags and HCANA version have not been matched to those settings**.
The parser ignores commented constants and includes active `00000` rows. A
matrix/constant/flag or parser convention mismatch remains a hypothesis, not a
diagnosed cause. Establish that provenance before attributing the difference
to geometry or iteration. No offset or parser adjustment was made to force closure.

## Geometry work and limits

The proof note identifies two algebraic target-projection mismatches: the
vertical slope omits vertex terms retained in xtar, and the horizontal slope
and intercept differ in their beam-position contribution to the longitudinal
vertex. It separately discusses saved-versus-recomputed xtar and finite replay
iteration. The laboratory vertical `react.y` belongs in the spectrometer x
equations; that coordinate mapping is not itself an error.

A freely fitted constant may absorb a nearly uniform offset; it does not prove
the target equations or replay conventions agree. Physical consequences remain
conditional on convention/provenance checks. **No target equations, truth values,
or reconstruction code were changed by the proof.**

Code: `sieve_afterburner.py`, `spectrometer_config.h`, `elastic_diagnostics.py`.
Three development tests in `tests/test_sieve.py` covered polynomial/projection
arithmetic, cuts and synthetic six-plot closure. The proof identities were checked
by substitution. Neither check establishes closure for actual farm data.

Retain `sieve/{manifest.json,tsv/,plots/,histograms.npz}` and the exact replay
configuration. Deferred substantive studies are weighted SVD using all development
cores, Huber loss and hole-centroid displacement; xtar-term changes and elastic
delta scans were deferred too. Nearby momenta must not be pooled under an assumed
identical transport map, and runs 1541–1543 were short false starts per the user.
