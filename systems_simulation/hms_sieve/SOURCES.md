# HMS geometry sources

This milestone validates coordinates and nominal positions. It is a marker-only
geometry, not a surveyed material model or a scattering simulation. Dimensions
listed below are source evidence; their listing does not imply that a physical
plate or foil has been constructed.

## Analyzer conventions and campaign inputs

- [HCANA `THcHallCSpectrometer.cxx`](https://github.com/JeffersonLab/hcana/blob/master/src/THcHallCSpectrometer.cxx):
  transport-coordinate conventions, spectrometer rotations, and pointing offsets.
  The accessible `master` source was inspected during this work. This moving
  reference is not a verified checkout of the campaign replay version.
- [Podd `THaSpectrometer.cxx`](https://github.com/JeffersonLab/analyzer/blob/master/Podd/THaSpectrometer.cxx):
  `SetCentralAngles` and `LabToTransport` were inspected in accessible upstream
  `master`. The analyzer revision used by the historical replay is OPEN.
- [HCANA `THcExtTarCor.cxx`](https://github.com/JeffersonLab/hcana/blob/master/src/THcExtTarCor.cxx):
  extended-target correction and projection to the sieve plane. This file and
  `THcReactionPoint.cxx` were **not freshly retrieved** in this task; their
  descriptions are inherited from the existing local geometry derivation below.
- [Existing geometry derivation](../../docs/HMS_GEOMETRY_CONSISTENCY.md):
  the campaign's equation-level discussion and references to HCANA revision
  `fcef8e23e16228fbb930078b9c44ef5031bed4f9`. That pinned upstream revision could
  not be fetched in this session; its equivalence to the inspected `master`
  must not be assumed.
- [Fit ntuple construction](../../make_fit_ntuple_from_gmm.C), lines 471–490 and
  591–602, directly establishes the saved-branch inputs and constructed targets.
  [Current Python fit](../../fit_elastic.py), lines 42–55, uses saved `xtar` in
  the matrix inputs and `xptarT`, `ytarT`, `yptarT` as fit targets.
  [Existing geometry diagnostic](../../diagnostics/validation/geometry/core_geometry.C),
  lines 271–277, evaluates `targetTruth` and the HMS x closure residual.
  These local sources were read directly; no fitting code was changed.
- [Nominal spectrometer profiles](../../spectrometer_profiles.def): HMS has a
  9 by 9 grid, nominal sieve plane at 168 cm, transport-x pitch 2.54 cm, and
  transport-y pitch 1.524 cm. Nominal centers are
  `x = (i - 4) * 2.54 cm`, `y = (j - 4) * 1.524 cm`, for indices 0 through 8.
  These are reconstruction conventions, not an independent mechanical survey.
- [Campaign sieve mask](../../HMS_6p667GeV/config/sieve_mask.json):
  blocked labels `[[2, 3], [5, 5]]`. This is an analysis mask. The physical
  drawing orientation has not been checked against these labels.

Use the configuration and generated provenance of this geometry for the selected
angle, beam position, target reference position, and pointing offsets. The
campaign's current choices are not universal HMS constants.

## Physical sieve references

[Hall C wiki, HMS Optics / HMS Collimator/Sieve](https://hallcweb.jlab.org/wiki/index.php/Main_Page#HMS_Collimator/Sieve)
was inspected as web text. It states:

| Quantity | Published value |
|---|---|
| Assembly location | 166.4 cm from the target, before the first quadrupole |
| Material | Mi-Tech HD-17 tungsten alloy |
| Density | 17 g/cm3 |
| Composition | 90% W, 6% Ni, 4% Cu |
| Sieve thickness | 1.25 in = 3.175 cm |

The location statement does not specify a front face, midplane, or another
assembly reference. In particular, it does not establish that the nominal
168 cm reconstruction plane is the physical plate center.

[Hall C Standard Equipment Manual](https://hallcweb.jlab.org/safety-docs/current/Standard-Equipment-Manual.pdf),
section 4.7.1, **HMS**, printed page 83, was inspected through the indexed text.
It gives outer sieve dimensions of 10.00 in vertically by 8.25 in horizontally
(25.4 by 20.955 cm), thickness 1.25 in, a reduced central aperture, and two
blocked holes. It describes the bottom ladder installation and the role of
shielding above the sieve in distinguishing the upper hole row. These details
must be considered before extending this model to material scattering.

[Hall C DocDB 1104, version 2, HMS sieve drawing](https://hallcweb.jlab.org/DocDB/0011/001104/002/HMS-Sieve.pdf)
is the intended mechanical authority. Retrieval was unsuccessful and the
drawing was **not visually inspected**. No bore dimensions, bore angles,
blocked-hole orientation, or plate-to-pattern offset in this model are claimed
to be verified from that drawing. See [open geometry questions](OPEN_GEOMETRY.md).

## Geant4 viewer API

The native coordinate renderer follows the official
[Geant4 11.4 application guide, Visualization User Actions and Standalone Visualization](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/Visualization/compiledcontrol.html#visualization-user-actions).
`G4Polyline`, `G4Circle`, `G4Text`, named user actions and scene extents were
checked against that guide and the corresponding visualization source. API
inspection does not substitute for compiling/running this viewer; that check
remains open because the local Geant4 development package is unavailable.
