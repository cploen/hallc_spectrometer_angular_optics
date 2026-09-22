# Open physical geometry questions

The delivered model represents the incident beam line, interaction reference,
center-foil reference, nominal sieve plane, and nominal hole centers. Markers
have no claimed physical thickness, material, or aperture. This is sufficient
to check coordinate transformations and conventions; it cannot predict material
interactions, acceptance, or transmission to the focal plane.

| Item | Established | Unresolved before a material simulation |
|---|---|---|
| Longitudinal sieve reference | Analysis profile uses 168 cm; wiki places assembly at 166.4 cm | Which physical face/reference each distance denotes; surveyed plate location |
| Nominal hole pattern | Analysis uses 9 by 9 centers with x/y pitches 2.54/1.524 cm | Mechanical verification of pattern and drawing view orientation |
| Blocked positions | Campaign label mask is `[[2,3],[5,5]]`; manual describes two blocked holes | Correspondence of labels to the installed plate's transport coordinates |
| Bore sizes | Manual identifies a smaller central hole | Ordinary and reduced diameters, tolerances, and any additional special bores |
| Bore directions | Not verified | Parallel versus angled bores; entry/exit center coordinates; chamfers or tapers |
| Plate extent | Manual gives 25.4 cm vertical by 20.955 cm horizontal, 3.175 cm thick | Outer-plate origin relative to the central hole, edge details, and mounting features |
| Neighboring shielding | Manual describes the ladder and shielding above the sieve | Relevant installed solids and gaps, particularly near the top row |
| Material | Wiki gives HD-17, density 17 g/cm3, 90% W / 6% Ni / 4% Cu | As-installed provenance and whether a more detailed material specification is needed |
| Center foil | Run 1544 metadata labels nominal Z=0 and “Carbon 0.5%”; only its center is represented | Physical thickness/areal density interpretation, transverse dimensions, tilt, and surveyed position |
| Beam | Configured incident line and vertex can be transformed | Whether settings reproduce the selected run's measured beam position, slopes, and raster |
| Pointing | Source defaults are known, but the inspected Podd transform and inherited Hall C extended-target convention imply different component frames | Actual replay call path, lab versus TRANSPORT meaning of the stored vector, active parameter overrides, and any independent sieve survey offset |

Required evidence is specific to each question: the installed plate drawing
and survey must establish distances, pattern orientation, bore geometry,
outline and mounting/shielding; target records must establish foil properties
and placement. The actual replay job's parameter files, HCANA and analyzer
revisions, optics matrix, and configured beam provider must establish beam
coordinates and pointing conventions. Equipment design values or upstream
defaults alone cannot close these run-specific questions. See the equations
and source-path conflict in [GEOMETRY_AUDIT.md](GEOMETRY_AUDIT.md).

Do not resolve the 166.4 cm versus 168 cm difference by silently assigning the
former to the entrance face: adding half of the published thickness gives
167.9875 cm, which is suggestive but does not establish that interpretation.
Likewise, do not turn the analysis mask into a surveyed blocked-hole placement
without checking the drawing's viewing direction and the installed orientation.

The center foil is a reference marker only. No guessed carbon slab dimensions
or thickness are introduced. Beam rays are geometric references and do not
constitute a model of electron scattering at the target.

The local environment has no usable GEANT4 runtime. The deterministic coordinate
tests pass; CMake stops at the missing Geant4 package. Native compilation and
inspection of the generated coordinate primitives remain to be performed in a
GEANT4-enabled environment such as iFarm. No GDML or material apertures are
constructed at this stage. No physics or yield validation is claimed.

The geometry supports a separate explicit translation for each setting. This
is an adjustment parameter, with its frame and provenance recorded, not evidence
that a fitted constant sieve-plane residual equals a physical mispointing.
The reported registration offset is held fixed across all synthetic rays of a
setting; target-intercept and slope differences are retained separately.

The source inspection status, including the unavailable mechanical drawing and
unverified pinned HCANA revision, is recorded in [SOURCES.md](SOURCES.md).
