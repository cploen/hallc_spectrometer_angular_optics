# HMS sieve-slit study: geometry and coordinate validation

This is a systems simulation, separate from `07_diagnostics`. It traces
deterministic geometric constructions from the optics and analyzer conventions.
**No sieve data, fitted matrix, particle transport, or scattering events are
used as closure tests.** The optics fitting/reconstruction code is unchanged.

The source audit, parameter model, ray comparisons, and visualization input
are implemented. Eleven automated tests pass, including a compiled check
against the actual `spectrometer_config.h`. Across two nominal settings,
168 endpoint-defined rays close within **1.78e-15 cm (1.78e-14 mm)**.
**Native Geant4 compilation and display remain unverified:** local CMake stops
at the missing `Geant4Config.cmake`. This stage therefore has not met the full
Geant4 display requirement.

## File guide

| Need | Start here |
|---|---|
| Understand the model, frames and three constructions | This README; [quantity/source audit](GEOMETRY_AUDIT.md) |
| Run tests and reproduce the committed artifacts | [TESTING.md](TESTING.md) |
| Transfer to iFarm, build, and display/export | [IFARM.md](IFARM.md) |
| Continue the work and identify the remaining gates | [NEXT_STEPS.md](NEXT_STEPS.md) |
| Distinguish established inputs from unresolved geometry | [OPEN_GEOMETRY.md](OPEN_GEOMETRY.md); [source inspection status](SOURCES.md) |
| Change a setting's angle or mispointing | [config/geometry.json](config/geometry.json); parameter section below |
| Inspect numerical results and the figure | [generated/RESULTS.md](generated/RESULTS.md); [coordinate_views.svg](generated/coordinate_views.svg) |
| Feed a future CAD renderer | [generated/scene.json](generated/scene.json), produced by [build_geometry.py](build_geometry.py) |

## What the model represents

The scene contains the nominal laboratory beam line, target center and center
foil marker, HMS central axis, laboratory/HMS/sieve coordinate triads, nominal
hole centers, and two separately named target reference planes:

- **Lab `Z=0`:** perpendicular to the incident beam through the nominal target.
- **HMS `z=0`:** perpendicular to the spectrometer central axis through its
  declared object origin. This is the target-intercept reference plane.

The sieve **projection** plane is HMS `z=168 cm`. Its association with a physical
plate face is unresolved. The equipment reference places the assembly before
the first quadrupole, but the reported 166.4 cm position lacks an established
face reference. No physical sieve or foil solid is invented. Plane outlines,
markers, and displayed axes have drawing sizes, not material dimensions.

All Python/JSON vectors are **`(x,y,z)`**, and laboratory vectors are `(X,Y,Z)`.
The user's longitudinal-first notation `(z_v,x_v,y_v)` is consequently stored
as `(x_v,y_v,z_v)`. Lengths are cm; the Geant4 reader multiplies by `cm` once.
Slopes are dimensionless; reported angular differences use `atan(slope)` in mrad.

## Three constructions and the offset test

**A** joins a declared lab vertex to a nominal sieve-hole center after an
explicit rigid lab-to-HMS transformation. Both the slope/intercept form and
an independent lab line/plane intersection must close numerically.

**B** evaluates the unchanged HMS `targetTruth` equations. Synthetic horizontal
displacement uses the declared test hypothesis `X=-xbpm`; this does not assert
the sign of an actual raster provider. It uses the source's angle-dependent
mispointing formulas even if a different geometric alignment is selected.

**C** traces the HCANA extended-target **coordinate stage** recorded in the
previous source audit and directly readable local mirror. It receives analytic
`p_old=p_new=p_A`, `q_A`, and `yt_A/100` at the matrix interface. Its x intercept
uses `-Y-mx-p_old*Z*cos(theta)` and projection uses `p_new`. No matrix is
evaluated and no reaction vertex is reconstructed. The y result tests unit
conversion and projection only; it is not an independent HCANA y-vertex
construction. A separate test prescribes unequal old/new slopes to check that
interface explicitly. Full HCANA/ReactionPoint behavior remains OPEN.

For each setting/construction, one constant sieve-plane registration offset is
defined by the central foil, central hole, zero-displacement case. The same
offset is subtracted from every other case. Nothing is refitted for individual
holes, foils, or beam displacements. Raw results are preserved alongside the
offset-subtracted values. Target intercepts and angular differences are also
reported: agreement up to an offset at one plane does not prove full ray
equivalence.

For the declared in-plane, transport-pointing benchmark, substitution gives

```
B: dx = -Y - mx;     dy = X*q_B*sin(theta)
C: dx = -X*p_A*sin(theta);  dy = 0 (A supplied the y interface)
```

Thus at zero transverse beam displacement B matches the grid up to the setting
offset `(-mx,0)`. The deliberate displacement tests expose terms that a single
setting offset cannot remove. These are geometric identities of the stated
constructions, not conclusions about measured scattering populations.

See [numerical results](generated/RESULTS.md),
[all individual rays](generated/closure.tsv), and
[coordinate figure](generated/coordinate_views.svg).

## Parameters and reproduction

[config/geometry.json](config/geometry.json) is the parameter source. Each entry
under `settings` has its own central angle, out-of-plane test angle, translation
model, and translation-frame interpretation. The included 12.490 and 15.195
degree angles are nominal metadata settings, not surveyed alignments.

To test a setting-specific geometric adjustment, set that entry's
`mispointing_model` to `explicit_translation` and supply
`pointing_translation_cm: [mx,my,0]`, with its source/status. Choose
`pointing_mode` explicitly: `transport_components` or `podd_lab_components`.
These represent the unresolved code-path interpretations; they must not be
mixed silently. A nonzero longitudinal translation is rejected because it is
not supported by the inspected HCANA setup. Changes affect the reference
geometry and displayed frame, not the existing fit equations.

Use [TESTING.md](TESTING.md) for commands, expected results, compiler selection,
and reproduction into a temporary directory. It also explains how to select
another setting without overwriting the committed baseline.

The generator writes `scene.json`, containing frame origins, basis vectors,
planes, and named primitives in laboratory coordinates. Both the SVG and
Geant4's `scene.tsv` derive from that same representation. A future CAD renderer
should consume this scene or the validated parameter module; do not duplicate
angles, translations, or hole coordinates manually. `validation.json` records
configuration/source hashes and numerical results.

Use [IFARM.md](IFARM.md) for transfer, environment requirements, native build,
file export and interactive viewing. The viewer draws coordinate primitives
and runs no events; native execution remains unverified.

## Gate before material scattering or acceptance claims

[GEOMETRY_AUDIT.md](GEOMETRY_AUDIT.md) contains the quantity/convention/source
table. [SOURCES.md](SOURCES.md) distinguishes freshly inspected upstream code,
inherited audit evidence, local implementation, and unavailable drawings.
[OPEN_GEOMETRY.md](OPEN_GEOMETRY.md) lists the unresolved values and the evidence
needed to resolve them.

Before a physical scattering model: run and inspect the native viewer; obtain
the mechanical drawing/survey defining plate faces, bores and installed
orientation; establish foil dimensions and placement; identify the actual
HCANA/Podd revisions, pointing meanings, beam-provider conventions and relevant
parameter overrides for each setting. Geometry closure alone cannot establish
those inputs, HMS acceptance, electron-cut survival, or scattering yields.
