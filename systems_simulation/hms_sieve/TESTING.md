# Testing and reproduction

[Overview and file guide](README.md) · [iFarm](IFARM.md) · [Next steps](NEXT_STEPS.md)

Use the **full repository checkout**, not this directory alone. The generator
records hashes of parent optics sources, and one test compiles the actual
`spectrometer_config.h`. The coordinate checks require Python **3.7 or later**
(standard library only) and a C++17 compiler. They do not require ROOT, NumPy,
Geant4, a matrix, or experimental data.

Last local verification (2026-09-22): 11 tests passed with no skips, including
compiled fit-source parity. Regeneration into a separate directory was
byte-identical to the committed baseline artifacts. The user subsequently
reported 11 tests passing with no skips on iFarm at `7fff412`, plus a successful
native Geant4 build, VRML export and Qt/OpenGL display on 2026-09-24. The exact
reported environment and limits are in [IFARM_VALIDATION.md](IFARM_VALIDATION.md).
The layered update passes 15 local tests on 2026-09-24: the original 11 plus
four annotation checks. Native validation of the new layers is still pending.

## Run the coordinate tests

From the repository root:

```sh
python3 --version
c++ --version
python3 -m unittest discover -s systems_simulation/hms_sieve/tests -v
```

Expected for this revision: **15 tests, OK, no skips** (11 coordinate checks and
four visualization-metadata checks). The source parity test skips if no compiler
or parent optics header is found. A skipped
test is incomplete source validation, even if unittest prints OK. A compiler
failure is a failure, not a skip. `CXX` can name a different compiler executable
(for example `CXX=g++`); do not put flags or a shell command in that variable.
The local macOS run used `CXX=/opt/homebrew/opt/llvm/bin/clang++` because the
system compiler could not locate its C++ standard headers.

Coverage: frame signs, handedness and round trips; distinct lab/HMS reference
planes; central/off-axis holes; foil Z=-8,0,+8 cm; deliberately displaced beam
positions; the unchanged fit construction; HCANA old/new slope interface;
one fixed offset per setting; explicit per-setting translations; and refusal
to invent a material sieve solid.

## Reproduce without overwriting the committed results

Run from the repository root (POSIX shell):

```sh
hms_check_dir=$(mktemp -d)
python3 systems_simulation/hms_sieve/build_geometry.py --output "$hms_check_dir"
diff -ru systems_simulation/hms_sieve/generated "$hms_check_dir"
```

For an unchanged checkout and the same parameters, expect no substantive
differences. Tiny floating-point formatting differences can occur between
platforms; inspect them rather than silently replacing the baseline. Source or
configuration edits deliberately change the hashes in `validation.json`.

The baseline contains **168 reference rays across two settings** (504 rows for
A/B/C). Maximum endpoint closure is approximately **1.78e-15 cm**; the tests
require closure below **1e-12 cm**. Nonzero B/C residuals are expected and
reported, not forced to zero. C is the bounded coordinate-interface trace
described in the README, not a complete HCANA replay.

`closure.tsv` retains raw and offset-subtracted residuals in cm/mm, target-plane
intercepts, slopes and angular differences. `scene.json`, `scene.tsv` and
`coordinate_views.svg` show the default setting. To display another setting:

```sh
python3 systems_simulation/hms_sieve/build_geometry.py \
  --setting HMS_15p195_code_defaults --output "$hms_check_dir/setting_15p195"
```

`--setting` selects the scene; the numerical comparison still covers **all**
settings in the selected configuration. Use `--config /path/to/geometry.json`
for a deliberate parameter experiment. Keep its outputs separate from the
baseline and record the source of any alignment adjustment.

## What passing means

The equations and declared geometric inputs are internally checked, including
agreement with the existing C++ fit construction. This does not certify a
surveyed apparatus, physical hole dimensions, historical replay parameters,
or a Geant4 build. Native compilation and visual inspection are separate checks
in [IFARM.md](IFARM.md). They passed for the reported `7fff412` baseline; the
new layers must be checked independently using [VIEWER.md](VIEWER.md).
