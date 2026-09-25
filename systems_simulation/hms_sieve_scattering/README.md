# Provisional HMS sieve-scattering run

This is now a separate Geant4 application layered on the closed-tested coordinate
model in `../hms_sieve`. It places a provisional Carbon 0.5% foil in the LAB
frame, the HMS sieve at the established 168 cm HMS position, and a sensitive
Q1-entrance scoring aperture. It does not edit or regenerate the reference
coordinate configuration. `case.py` checks its SHA-256, imports its frame
transform, and writes the flat runtime input consumed by the C++ application.

## First run setup

The primary electron starts at `LAB (x,y,z)=(x_beam,y_beam,-5 cm)` with direction
`+LAB Z`, 10.6 GeV energy, and target-plane centroid `(+0.22,-0.31) mm`.
Independent uniform raster offsets of ±1 mm are added in LAB X and Y. The
intrinsic spot width is provisionally zero because no pass-matched harp sigma has
been supplied. The source emits 10,000 primaries. Current does not change a
single-particle Geant4 track; the first run uses unit event weights so timing
and process distributions can be checked before applying a current-based rate
normalization.

The carbon slab is 0.032 in thick, a provisional guess from the thinner dummy
foil dimension, with effective density chosen to preserve the target report's
`0.1749 g/cm²` areal density. Its 4 by 4 cm transport patch only contains the
beam envelope; it is not the target's asserted physical outline. The target is
modeled as pure carbon at the derived density; its unspecified 0.05% remainder
is omitted.

The sieve plate is centered at HMS `z=168 cm`, with the drawing's thickness,
outer dimensions, 9 by 9 bore pattern, alloy recipe, reduced center bore, and
two blocked positions. Mapping drawing-horizontal to HMS x and drawing-vertical
to HMS y is an explicit provisional orientation. Mispointing is the existing
HMS code-default translation at 12.5 degrees; it remains distinct from the
beam centroid.

The Hall C Standard Equipment Manual gives Q1's inner pole radius as 25 cm,
actual length as 2.34 m, and vacuum-vessel inner radius as 20.05 cm. This first
run uses the documented 20.05 cm vacuum-vessel radius for a thin sensitive
scoring disk just 0.1 cm downstream of the provisional sieve's downstream
face. That disk is a Q1 entrance proxy because the manual does not set its
position relative to this sieve configuration. Q1 is not otherwise modeled:
there are no magnetic fields, Q1 material, or detector response. The radius is
a geometric scoring cut, not a simulated HMS acceptance.

The hit file reports momentum and
`delta = (p_at_scoring_plane - 6.667 GeV/c) / 6.667 GeV/c`, with a flag for
`|delta| <= 0.20`. This is the HMS momentum window, not a spread on the 10.6 GeV
incident beam. Geant4 process names in `material_steps.tsv` identify processes
that acted within carbon or sieve material. They are not the experiment's
analysis `procID` codes; a mapping would need an explicit definition.

## Build and run on iFarm

From the repository root:

```sh
module load geant4/11.2.1
hms_repo_dir=$PWD
hms_case_dir="$hms_repo_dir/systems_simulation/hms_sieve_scattering"
hms_run_dir=$(mktemp -d)

python3 "$hms_case_dir/case.py" \
  --write-g4-config "$hms_run_dir/geant4_input.tsv"
cmake -S "$hms_case_dir" -B "$hms_run_dir/build"
cmake --build "$hms_run_dir/build" --parallel 2

cd "$hms_run_dir"
/usr/bin/time ./build/hms_sieve_scattering geant4_input.tsv
```

The 10,000-primary count comes from the case file. Output is written in the run
directory:

- `material_steps.tsv`: track IDs, particle codes, process names, energies,
  and positions for steps in the target and sieve.
- `q1_entrance_hits.tsv`: tracks crossing the sensitive disk, with creator
  process, momentum, HMS delta, and LAB position/direction.

An empty sieve section or no Q1 hits is a possible result for 10,000 incident
beam primaries: the undeflected beam is along LAB Z, while the HMS is at 12.5
degrees. At the 168 cm HMS plane, the nominal LAB beam ray crosses about 37.2 cm
from the HMS origin vertically, outside the sieve's 12.7 cm half-height, so it
misses the plate unless target interactions redirect it. Large-angle target
scattering is rare. Such a result tests build
and runtime, but is not evidence that sieve scattering cannot occur. A later
high-statistics study may need importance sampling or a second, explicitly
outgoing-electron source mode to sample target-to-sieve rays efficiently.

## Provenance and unresolved items

| Input | First-run value | Status |
|---|---:|---|
| Coordinate transform and mispoint | Frozen `hms_sieve` transform; 12.5° HMS code default | Reused and hash-checked; not survey alignment |
| Incident beam | 10.6 GeV, +LAB Z, source at LAB z = −5 cm | User-selected provisional energy/source distance |
| Beam centroid and raster | (+0.22, −0.31) mm; ±1 mm/axis | From supplied email; run applicability TBD |
| Carbon 0.5% thickness | 0.032 in | Provisional dummy-foil guess; areal density preserved |
| Sieve center | HMS z = 168 cm | Provisional placement requested for first run |
| Drawing orientation | Horizontal→HMS x; vertical→HMS y | Provisional; installed orientation still to verify |
| Q1 sensitive aperture | Radius 20.05 cm | Manual vacuum-vessel inner radius; simplified scorer |
| Q1 longitudinal position | 0.1 cm after sieve downstream face | Provisional Q1 entrance proxy |
| Magnetic transport and cuts | Not modeled | Outside this run's scope |
| Beam-current rate scale | Not applied | Add after selecting the run current and normalization |

The [Hall C Standard Equipment Manual](https://hallcweb.jlab.org/safety-docs/March2026/Standard-Equipment-Manual-2026.pdf)
describes the HMS Q1 geometry and asymmetry. The corresponding [JLab HMS Q1
parameters page](https://www.jlab.org/Hall-C/VRML/hmsqm.html) reports the 2.34 m
actual length, 25.0 cm inner pole radius, and 20.05 cm Q1 vacuum-vessel inner
radius used above. Geant4's [EM reference](https://geant4.web.cern.ch/documentation/pipelines/master/prm_html/PhysicsReferenceManual/electromagnetic/elastic_scattering/singlescat.html)
describes its single Coulomb-scattering process; `G4EmStandardPhysics_option4`
is selected for this first mechanism-presence run.

The first pass is deliberately a mechanism and runtime check, not a yield
prediction. It does not establish transport through the HMS magnets, detector
cuts, or the actual experimental `procID` classification.
