# Next steps and completion gates

[Overview and file guide](README.md) · [Testing](TESTING.md) · [iFarm](IFARM.md)

The deterministic coordinate model, source audit, per-setting translations,
A/B/C comparisons, numerical outputs and SVG exist. Eleven coordinate/source
tests passed on the baseline, including the user-reported iFarm run at
`7fff412`. Native build, Qt/OpenGL display and VRML export also succeeded there;
see [IFARM_VALIDATION.md](IFARM_VALIDATION.md). Do not treat those checks as
evidence of physical alignment or scattering.

1. **Validate the layered viewer update on iFarm.** Follow [IFARM.md](IFARM.md)
   and [VIEWER.md](VIEWER.md). Record the new commit and environment. Inspect
   the clean default view, each independently selected frame, variable/branch
   naming views and optional A/B/C layers. Check that switching annotations
   leaves origins, planes, holes and rays fixed. Inspect labels in Qt/OpenGL;
   VRML export does not support `G4Text`. The prior native validation applies
   to the baseline viewer, not automatically to these new controls.
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
and commit of each completed check in the validation record. Keep historical
generated numerical artifacts distinct from later viewer-validation evidence.
