# Replay beam-position compatibility

The supplied runs 3283–3286 were reported by the probe to lack
`P.rb.raster.fr_xbpm_tar` and `P.rb.raster.fr_ybpm_tar`. The branch listing
for run 3283 contains `P.react.x`, `P.react.y`, and `P.react.z`.
This is treated as a replay output/configuration difference, not evidence of
a fundamental HMS/SHMS difference. The precise replay configuration is unverified.

Holly's SHMS_optics source/src/myEvent.cpp binds P.react.x/y to xVer/yVer;
source/shms_optics.cpp uses them in the physical vertex equations.
Our SHMS geometry already follows those expressions. No geometry, signs,
foil positions, or cuts are changed by this compatibility adjustment.

- make_fit_ntuple_from_gmm.C: retains commented SHMS BPM bindings, uses the
  existing required react bindings, and removes the redundant SHMS BPM requirement.
- make_yscol_candidate_tree.C: same input adjustment; checks availability of
  optional raster diagnostics and omits unfilled SHMS BPM output placeholders.
- diagnostics/validation/replay/probe_replay_inputs.C: retains commented SHMS
  BPM checks; required react fields continue to be validated.
- tests/test_replay_probe.py: synthetic SHMS input now omits BPM branches while
  retaining react fields, matching the observed schema.

The established HMS calculation still reads BPM fields. This preserves its
existing behavior; it is not a claim that other HMS replays cannot use react
fields. The X candidate reader already treats BPM branches as optional.
Real replay validation must still be rerun on ifarm.
