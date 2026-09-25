# iFarm validation record

[Overview](README.md) · [iFarm commands](IFARM.md) · [Viewer controls](VIEWER.md)

## Baseline: user-reported success, 2026-09-24

This record transcribes the user's iFarm handoff and attached Qt/OpenGL
screenshot. These commands were not rerun by the local coding agent. The
tested commit was **`7fff412`** on `optics-validation-studies`, titled
`Add Huber expansion and equal-foil studies with sieve comparison plots`.
The checkout was:

```text
/work/hallc/nps/cploen/analysis/hallc_spectrometer_angular_optics
```

The screenshot identifies the display host as `ifarm2402.jlab.org`.

| Component | Reported version or result |
|---|---|
| Python | 3.9.25 |
| GCC/C++ | 11.5.0 |
| Geant4 module | `geant4/11.2.1`; `geant4-config --version` returned `11.2.1` |
| CLHEP | 2.4.6.4 |
| Qt | 5.15.13 |
| Deterministic tests | 11 tests in 0.517 s; `OK`; no skips |
| Native build | `hms_geometry_viewer` built successfully |
| File export | `g4_00.wrl`, reported size 42K; 114 coordinate primitives |
| Interactive display | Qt/OpenGL succeeded; axes, planes, grid and A/B/C were visible |

The test command, from the repository root, was:

```sh
python3 -m unittest discover -s systems_simulation/hms_sieve/tests -v
```

The reported noninteractive build and export used:

```sh
module load geant4/11.2.1
geant4-config --version
hms_repo_dir=$PWD
hms_build_dir=$(mktemp -d)
cmake -S "$hms_repo_dir/systems_simulation/hms_sieve" -B "$hms_build_dir" \
  -DHMS_WITH_INTERACTIVE=OFF
cmake --build "$hms_build_dir" --parallel 2
cd "$hms_build_dir"
./hms_geometry_viewer \
  "$hms_repo_dir/systems_simulation/hms_sieve/generated/scene.tsv" \
  "$hms_build_dir/export_vrml.mac"
```

CMake automatically found
`/u/group/halla/apps/Geant4/11.2.1/el9/lib64/cmake/Geant4/Geant4Config.cmake`.
The handoff reports a subsequent successful interactive launch and supplies
its screenshot; it does not supply that launch's exact build/command log.

`VRML2FILE` warned that `G4Text` is not implemented. Geometry export succeeded,
but text labels are unsupported in that format. The Qt/OpenGL image showed
heavy label overlap, which motivates the annotation-layer update.

## Scope and current pending check

The baseline confirms deterministic coordinate/source checks and native
geometry rendering. It does **not** establish surveyed sieve dimensions,
historical replay geometry, scattering, acceptance, electron-cut survival or
material transport. No experimental sieve events enter the closure tests.

**The new layered viewer has not yet been validated natively on iFarm.**
Its geometry inputs and equations are unchanged. Follow [VIEWER.md](VIEWER.md)
to check independent frame toggles, naming layers, diagnostic rays and clean
overview, and record the new commit and result here. The baseline's 114
primitives and 42K export describe that tested version, not mandatory output
counts or file sizes for later visualization versions.

Generated numerical reports retain their original local-runtime status;
this dated record supplies the later external validation without regenerating
or changing the geometric benchmark.
