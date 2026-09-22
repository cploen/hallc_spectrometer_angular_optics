# Next steps and completion gates

[Overview and file guide](README.md) · [Testing](TESTING.md) · [iFarm](IFARM.md)

The deterministic coordinate model, source audit, per-setting translations,
A/B/C comparisons, numerical outputs and SVG exist. Eleven coordinate/source
tests pass. **The overall geometry milestone remains incomplete because native
Geant4 build/display validation is outstanding.** Do not treat compilation alone
as evidence of physical alignment or scattering.

1. **Run and inspect the native viewer.** Follow IFARM.md. Record versions and
   commands, preserve the exported scene/image, and check axis directions,
   origins, both target reference planes, hole labels and displayed ray
   intersections. Resolve any implementation mismatch before physics work.
2. **Trace remaining coordinate code with explicit inputs.** Recover the exact
   HCANA/Podd sources needed to verify the inherited ExtTarCor/ReactionPoint
   statements and pointing-vector call path. Extend synthetic tests only from
   those equations. Continue using prescribed rays, not sieve events, as
   geometric closure tests; C currently tests only the stated coordinate stage.
3. **Establish each setting's alignment convention.** Keep translations and
   angular offsets distinct. Record whether a translation is an exploratory
   adjustment, code default, replay parameter or surveyed value. One offset is
   held fixed across the complete hole/foil/beam-displacement scan for that
   setting. Preserve raw residuals and intercept/slope differences.
4. **Resolve physical geometry before introducing material.** Obtain the
   drawing/survey establishing sieve face location, bores, plate orientation,
   blocked-hole placement and nearby structure, plus foil dimensions and
   placement. Do not identify 166.4 cm or 168 cm with a physical face by inference.
   [OPEN_GEOMETRY.md](OPEN_GEOMETRY.md) lists the missing evidence item by item.
5. **Only then define the next physics stage.** A material-scattering model
   requires those physical inputs and a separate acceptance/transport plan.
   The current model establishes neither survival of electron cuts nor
   transmission to the focal plane. Detailed fields, detector volumes and
   production physics lists remain outside this milestone.

For any later CAD view, consume `generated/scene.json` or the validated geometry
module. Do not maintain a second set of geometry constants. Document the scope
of each completed check and update the status in the README and numerical report
generator when the corresponding gate is actually satisfied.
