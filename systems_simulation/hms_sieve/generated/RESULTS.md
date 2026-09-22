# Deterministic geometry results

All inputs are declared geometric vertices, hole labels and code parameters. No sieve data, fitted matrix or simulated scattering events enter these checks.

One constant sieve-plane registration offset is fixed by the central foil / central hole / zero-displacement case for each construction and setting. It is then held fixed across the scan. Offsets are reported, not applied to the geometry or source equations.

| Setting / construction | Offset x,y (mm) | Max raw x,y (mm) | Max after offset x,y (mm) | Max angle difference p,q (mrad) |
|---|---:|---:|---:|---:|
| HMS_12p490_code_defaults / A_endpoint | 0, 0 | 1.776357e-14, 1.776357e-14 | 1.776357e-14, 1.776357e-14 | 0, 0 |
| HMS_12p490_code_defaults / B_fit_target | -1.48306, 0 | 3.48306, 0.02147502 | 2, 0.02147502 | 2.174336, 0.01336946 |
| HMS_12p490_code_defaults / C_HCANA_coordinate_step | 0, 0 | 0.02782667, 1.776357e-14 | 0.02782667, 1.776357e-14 | 0, 0 |
| HMS_15p195_code_defaults / A_endpoint | 0, 0 | 1.776357e-14, 1.776357e-14 | 1.776357e-14, 1.776357e-14 | 0, 0 |
| HMS_15p195_code_defaults / B_fit_target | -1.340296, 0 | 3.340296, 0.02716525 | 2, 0.02716525 | 2.084039, 0.01689775 |
| HMS_15p195_code_defaults / C_HCANA_coordinate_step | 0, 0 | 0.03365659, 1.776357e-14 | 0.03365659, 1.776357e-14 | 0, 0 |

The raw TSV retains every intercept, slope, cm/mm residual and angular difference. `D(L)-D(0)=L*(delta_p,delta_q)` is checked independently; an offset at one plane does not prove full ray equivalence.

C is only the HCANA coordinate stage traced with analytic inputs: A supplies p_old=p_new, q and ytar/100. Its y closure is a unit/projection check, not an independent HCANA y-vertex reconstruction. A separate automated test preserves the distinction between old and updated slopes.

The source default translations are setting dependent. They are not measured alignments. See the per-setting configuration and OPEN_GEOMETRY.md before assigning a physical offset.

Native Geant4 compilation/display: NOT VERIFIED (local Geant4 package unavailable). Physical survey, active replay configuration, ReactionPoint reconstruction and matrix iteration: OPEN / not exercised.
