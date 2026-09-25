# HMS sieve scattering case: provisional 12.5° setup

This directory starts the material-geometry/scattering study without altering
`../hms_sieve/config/geometry.json` or its closed-tested coordinate equations.
The case is an overlay: `case.py` imports the existing frame transform and
checks the reference configuration's SHA-256 before creating the 12.5° frame.
If the reference changes, the overlay refuses to load until it is deliberately
reviewed and retested.

## Provisional target thickness

The selected target is Carbon 0.5% at nominal lab `Z=0`. The target report gives
99.95% carbon and `0.1749 ± 0.00035 g/cm²` areal density. It does not state the
physical thickness. At Christine's request, this first geometry guess borrows
the thinner dummy-foil certificate dimension cited in that report: `0.032 in`
(`0.08128 cm`). The material's effective density is derived as
`0.1749 / 0.08128 = 2.15182 g/cm³`, so the model preserves the report's measured
areal density. This is a provisional geometric slab, not a claim that the
Carbon 0.5% target was measured to be 0.032 in thick. The other cited dummy
foil dimension, `0.050 in`, is retained as an alternative for a later
thickness-sensitivity comparison. The target's transverse outline and tilt
remain unresolved.

## Mispointing and beam input

The frame uses the existing transport-component interpretation and the HMS
source-code default pointing translation evaluated at `12.5°`:
`(0.14825, 0.06825, 0) cm`. This preserves and exercises the existing
mispointing convention. It is not asserted to be the selected run's active
replay alignment; a per-run additional translation remains unset until replay
parameters or the planned ROOT/BPM comparison establish it.

The supplied email's nominal target centroid is recorded separately as
`(+0.22, -0.31) mm` in LAB `X,Y`, with the BPM x sign already flipped as Dave
described. The raster is a static uniform rectangle of ±1 mm per axis. Harp
sigmas remain null until a pass-matched scan is selected. These beam inputs are
not merged into the spectrometer pointing translation.

## Sieve dimensions and limits

Drawing `6732-C-5681` supplies a 10.00 by 8.25 in plate envelope, 1.25 in
thickness, 9 by 9 nominal pattern, 1.000 in vertical and 0.600 in horizontal
pitch, 78 nominal 0.20 in holes, one 0.10 in center hole, and two positions
without holes. The pattern center is 0.25 in above the plate's geometric
center in the drawing's front view. Hole-index orientation in HMS coordinates
and the absolute physical front/center/back location remain separate from this
drawing and are therefore explicitly unresolved in the case file. No material
solid is placed at 166.4 or 168 cm by assumption.

The intended first transport endpoint is a scoring surface just downstream of
the sieve at the Q1 entrance. It can establish that a sieve interaction can
produce an electron with specified position, direction, and momentum there.
Without HMS field transport and detector response, it does not establish full
spectrometer acceptance or the survival of detector-based electron cuts.

## Checks

Run the focused overlay checks from the repository root:

```sh
python3 -m unittest discover \
  -s systems_simulation/hms_sieve_scattering/tests -v
```

These tests verify the reference hash, 12.5° source-default mispointing, target
areal-density preservation, and explicit TBD values. They do not build or run
the eventual Geant4 material-scattering transport.
