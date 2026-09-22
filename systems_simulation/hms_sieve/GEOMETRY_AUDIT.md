# HMS beam, target, and sieve geometry audit

Source audit: 2026-09-21; deterministic comparison updated 2026-09-22. This directory is a systems simulation, separate from
the optics diagnostics. This note records the conventions needed before a
physical GEANT4 geometry can be identified with the historical replay.

**Status: OPEN.** A geometric construction and closure test can be evaluated
under explicit conventions. They cannot yet be certified as the geometry of
the actual replay. In particular, the inherited analyzer transform and the
Hall C extended-target convention do not interpret the stored pointing vector
the same way in the sources available to this review. Do not silently resolve
that conflict or change the existing fit equations to make a comparison close.

## Evidence and provenance

- **H1 — fresh upstream source review:** official
  [HCANA `master`, THcHallCSpectrometer.cxx](https://github.com/JeffersonLab/hcana/blob/master/src/THcHallCSpectrometer.cxx),
  read through the web tool by the task coordinator on 2026-09-21. The fetched
  result was cached approximately three weeks earlier. This is a moving
  upstream source, **not a verified historical replay revision**. Relevant
  ranges in that fetched version: 286–289 (HMS angle sign), 295–325 (default
  mispointing), 367–379 (central angles and stored pointing vector), and
  508–541 (matrix inputs and output offsets). Line numbers may move on master.
- **P1 — fresh upstream dependency review:** official
  [analyzer `master`, Podd/THaSpectrometer.cxx](https://github.com/JeffersonLab/analyzer/blob/master/Podd/THaSpectrometer.cxx),
  methods `SetCentralAngles` and `LabToTransport`, read by the coordinator in
  this task. The dependency revision linked to the actual replay is unknown.
- **L1 — existing local audit, indirect evidence for pinned HCANA:**
  [HMS_GEOMETRY_CONSISTENCY.md](../../docs/HMS_GEOMETRY_CONSISTENCY.md),
  dated 2026-09-14. It records an earlier review of HCANA
  `fcef8e23e16228fbb930078b9c44ef5031bed4f9`. The pinned files could not be
  retrieved in this task. Statements below about `THcExtTarCor` and
  `THcReactionPoint` are inherited from this local audit, not a fresh reading.
- **F1 — current fit implementation, directly read locally:**
  [spectrometer_config.h](../../spectrometer_config.h), especially lines
  23–38 and 62–85, and
  [spectrometer_profiles.def](../../spectrometer_profiles.def).
  The existing code has not been changed by this audit.
- **M1 — campaign metadata:**
  [DATfiles/list_of_optics_run.dat](../../DATfiles/list_of_optics_run.dat)
  and campaign configuration. Metadata labels nominal settings; it is not a
  substitute for survey records or the active replay parameter files.
- **F2 — saved branches, fit inputs and diagnostic, directly read locally:**
  [make_fit_ntuple_from_gmm.C](../../make_fit_ntuple_from_gmm.C), lines 471–490
  and 591–602;
  [fit_elastic.py](../../fit_elastic.py), lines 42–55;
  [fit_opt_matrix_gmm.C](../../fit_opt_matrix_gmm.C), lines 575–577; and
  [core_geometry.C](../../diagnostics/validation/geometry/core_geometry.C),
  lines 271–277. These establish the distinction between saved reconstruction
  quantities and constructed optimization targets.

## Convention and provenance table

All axes below form right-handed systems. Capital letters denote laboratory
coordinates; lower-case letters denote HMS TRANSPORT coordinates. Numerical
lengths in the analysis convention are cm; convert once to GEANT4 units.

| Quantity | Meaning | Frame | Origin/reference plane | Axis/sign convention | Units | Value/formula | Source | Confidence/conflict |
|---|---|---|---|---|---|---|---|---|
| Laboratory axes | Beam and hall coordinates | Lab | Nominal target/object reference point | +Y upward; +Z downstream along incident beam; +X beam left | cm | Right-handed Cartesian basis | H1 coordinate comments; L1 | Upstream convention established; actual beam-provider branch mapping OPEN |
| HMS TRANSPORT axes | Spectrometer ray coordinates | HMS | Spectrometer object reference point | +x downward at zero out-of-plane angle; +z along central ray; +y completes right-handed basis | cm | Basis vectors in the next section | H1; P1 | Rotation established for inspected upstream sources; actual replay dependency OPEN |
| Laboratory target plane | Plane perpendicular to nominal beam | Lab | Lab Z=0 through nominal object point | Normal +Z | cm | Z=0; physical central foil is placed here only by nominal assignment | Geometry definition; M1 | Distinct from TRANSPORT target plane |
| TRANSPORT target plane | Track-intercept reference plane | HMS | HMS z=0 through object origin | Normal along HMS central ray | cm | z=0; xtar and ytar are intercepts here | P1; L1 | Not the laboratory Z=0 plane at nonzero central angle |
| Object origin | Reference point of the unshifted spectrometer coordinates | Lab/HMS | Nominal object point, before a specified mispoint translation | Origin and translation must be defined independently of beam position | cm | O_lab=(0,0,0) is a coordinate choice, not a survey measurement | Geometry definition | Survey relationship to target ladder OPEN |
| Signed central angle | Horizontal direction of central ray | Lab | Object origin | HMS angle negative in inspected HCANA; positive campaign angle is its magnitude | degrees | alpha=-theta for zero corrections; HCANA first enforces HMS sign, then adds fThetaCentralOffset in degrees | H1 286–289, 367–379; F1 | Actual angle and offset parameter files OPEN; do not reuse unsigned fit angle as signed lab angle |
| Out-of-plane central angle | Vertical direction of central ray | Lab | Object origin | P1 central unit vector is (sin(alpha)cos(phi), sin(phi), cos(alpha)cos(phi)) | degrees at SetCentralAngles call | H1 forms ph=fPhi_lab+fPhiOffset*RadToDeg and calls SetCentralAngles(...,false) | H1 367–379; P1 | H1 loads fOopCentralOffset earlier but this call uses fPhiOffset; preserve this distinction; actual values OPEN |
| Foil positions | Physical planes of target material | Lab | Nominal target reference along beam | Positive downstream Z; normal to beam only if assigned so explicitly | cm | Center foil nominal Z=0; other zfoil values are lab Z coordinates | M1; F1 | Center assignment allowed; dimensions, thickness, material, tilt and survey offsets must be separately sourced |
| Sieve reconstruction plane | Plane used for straight-line replay projection | HMS | Distance along HMS z from its object reference | Positive toward sieve/magnets | cm | L=168; xs=xtar+L*xptar; ys=ytar+L*yptar when both lengths are cm | L1 citing pinned ExtTarCor; F1 profile | Reconstruction convention established indirectly; physical face interpretation OPEN |
| Physical sieve faces | Upstream/downstream surfaces of material | Sieve/HMS | Must be referred to the same object origin as L | Thickness along physical plate normal | cm | No certified face location in this audit | No verified survey/drawing here | OPEN: 166.4 cm must not be equated with a front face or 168 cm with mid-plane without a thickness/drawing reference |
| Hole coordinates | Nominal grid centers in analysis | HMS sieve plane | Nominal centered sieve grid | x row increases downward; y column follows HMS +y | cm | x_i=(i-4)*2.54; y_j=(j-4)*0.6*2.54; i,j=0..8 | F1 23–26; HMS profile | Analysis grid known; drilled positions, missing/blocked holes, diameters, bores and plate survey OPEN in this coordinate audit |
| Beam/raster position | Incident beam position at a stated lab Z | Lab/provider | Beam-provider reference plane must be identified | Physical beam X,Y versus branch signs not assumed | cm in analysis | Fit uses B=-H.rb.raster.fr_xbpm_tar; vertical input react.y; fit HMS branch does not use react.x | F1 62–85; L1 | OPEN: B=X_lab is not established; raster/BPM and reaction coordinates need actual provider configuration |
| Beam direction | Incident beam slopes | Lab/provider | Same reference plane as position | dX/dZ and dY/dZ | dimensionless | Straight reference beam may be explicitly assigned (0,0,1); measured slopes require provider values | Geometry definition | Zero slope is a nominal input, not a recovered run measurement |
| Pointing translations | Position displacement of spectrometer axis/origin | Frame is disputed between code paths | Object plane or lab origin depending code path | H1 stores (mx,my,0); P1 subtracts it before rotation; Hall C vertical convention subtracts mx after rotation | cm | See conflict and default formulas below | H1 295–325, 367–379; P1; L1 | OPEN CONFLICT; do not treat stored components as both lab and transport displacements |
| Angular reconstruction offsets | Additions to reconstructed slopes | Matrix target outputs | No physical position origin | xptar gets fPhiOffset; yptar gets fThetaOffset | dimensionless, conventionally radians | xptar=sum[0]+PhiOffset; yptar=sum[2]+ThetaOffset | H1 508–541 | Separate from pointing translations; H1 also uses PhiOffset in central ph, so source pathways must be retained |
| Central-angle offsets | Changes to coordinate-frame angles | Lab/HMS frame rotation | Object origin | fThetaCentralOffset changes central in-plane angle | radians internally, converted to degrees for call | theta_lab += ThetaCentralOffset*RadToDeg | H1 367–379 | Do not substitute htheta_offset, a matrix-slope offset, for this parameter |
| Matrix length units | Polynomial inputs/output | Focal plane and target | Matrix reference planes | hut x,y and xtar divided by 100 | meters internally | ytar=sum[1] is meters; replay sieve projection uses 100*ytar | H1 508–541; L1 | cm/m conversion explicit; no physical length inferred from coefficient names |
| Reconstructed momentum offset | Addition to matrix momentum output | Reconstruction | Central momentum reference | delta=sum[3]+DeltaOffset in matrix evaluator | Matrix convention | Must track conversion to published percent delta separately | H1 508–541 | Actual central momentum, active offset and percent conversion not certified by this audit |
| Historical replay configuration | Executable and parameters that generated the tree | All | Actual production job | Overrides control active values | revision/parameter provenance | Not yet recovered | L1 and present review | OPEN; upstream defaults are not evidence of active run settings |

## Rigid geometry and translation conflict

At zero out-of-plane angle, P1 gives these basis vectors expressed in lab
coordinates, for the signed central angle alpha:

```
ex = (0, -1, 0)
ey = (cos(alpha), 0, -sin(alpha))
ez = (sin(alpha), 0,  cos(alpha))
R  = [ex ey ez]             # columns are transport basis vectors in lab
```

For an explicitly lab-defined origin O, the physical rigid transform is

```
r_transport = transpose(R) * (r_lab - O_lab)
r_lab       = O_lab + R * r_transport
```

Vectors/directions use the rotation only. Translations are not applied to
momenta, slopes or unit directions. With alpha=-theta, this means

```
x = -(Y-OY)
y = (X-OX)*cos(theta) + (Z-OZ)*sin(theta)
z = -(X-OX)*sin(theta) + (Z-OZ)*cos(theta)
```

Two different interpretations of H1's stored `(mx,my,0)` must remain separate:

1. **P1 `LabToTransport`:** subtract the stored vector before rotating,
   `transpose(R)*(r_lab-fPointingOffset)`. This treats its components as lab
   displacement coordinates. For `(mx,my,0)` the vertical coordinate is
   `x=-Y+my`, rather than `-Y-mx`.
2. **Hall C convention recorded by L1:** the extended-target code uses
   `x_v=-vertex.y-pointingOffset.x`. This treats the x component as a
   transport-vertical translation. Consistent rigid geometry for a full
   transport translation m would be `transpose(R)*r_lab-m`, equivalent to
   the lab origin `O_lab=R*m`.

This is a source-path interpretation conflict; it is **not** proof that the
historical HMS event reconstruction calls the inherited method in this way.
Establish the actual call path and parameter meanings before choosing a
physical alignment from these values. For a purely transverse transport
translation `(mx,my,0)`, the translated origin still lies in the same geometric
transport z=0 plane. A general lab translation can move that plane.

The angle-dependent default formulas recorded by H1 and used explicitly by F1
are, with the positive angle magnitude in degrees:

```
ax = min(abs(angle_deg), 50)
ay = min(abs(angle_deg), 40)
mx_cm = 0.1 * (2.37 - 0.086*ax + 0.0012*ax*ax)
my_cm = 0.1 * (0.52 - 0.012*ay + 0.002*ay*ay)
```

They are defaults, not a survey of the current apparatus. F1 evaluates them
from the campaign angle; replay can use parameter overrides and a different
sequence of angular corrections.

## Independent closure and comparison with existing fit/replay

For one selected, explicit rigid convention, transform a lab vertex to
`(xv,yv,zv)` and define a hole point `(xh,yh,L)`. The straight ray is exactly

```
xp = (xh-xv)/(L-zv)       yp = (yh-yv)/(L-zv)
xt = xv-xp*zv            yt = yv-yp*zv
xs = xt+L*xp             ys = yt+L*yp
```

Require `xs-xh=0` and `ys-yh=0` to floating-point tolerance. Round-trip
lab/transport coordinates, right-handedness, orthonormality and unchanged
vector norm are separate checks. These tests validate the stated geometric
model, not whether its input alignment matches the apparatus.

Keep the current F1 HMS fit equations as a separate, unchanged comparison.
With `s=sin(theta)`, `c=cos(theta)`, `B=-xbpm`, `z0=zfoil*c`, they are

```
yc   = zfoil*s + B*c - my
yp_B = (yh-yc)/(L-z0)
yt_B = zfoil*(s-yp_B*c) + B*(c+yp_B*s) - my
xp_B = xh/(L-z0)
xt_B = -reacty-mx-xp_B*z0
```

Direct substitution gives the independent, nonzero residuals

```
(xt_B + L*xp_B) - xh = -reacty-mx
(yt_B + L*yp_B) - yh = B*yp_B*sin(theta)
```

The first combines a slope numerator with no vertical vertex displacement
and an intercept that includes it. The second uses `z0` in its slope but
`z0-B*sin(theta)` in its intercept. Neither identity depends on a fitted
matrix. They do not by themselves quantify an error in the production replay.

For the **replay-style projection comparison**, L1 records pinned
`THcExtTarCor` using `x_tg=-vertex.y-pointingOffset.x`, updating xtar with
the reconstructed lab vertex Z and current xptar, reevaluating the matrix,
then projecting `xs=xtar+168*xptar`, `ys=100*ytar_matrix+168*yptar`.
That iteration uses a reconstructed vertex and a matrix, not a known physical
beam-to-hole ray. The available local afterburner mirrors a five-update limit
and 2 mrad stopping criterion. It must not be presented as exact rigid geometry.

Source links retained from L1:
[pinned THcExtTarCor](https://github.com/JeffersonLab/hcana/blob/fcef8e23e16228fbb930078b9c44ef5031bed4f9/src/THcExtTarCor.cxx#L112),
[pinned THcReactionPoint](https://github.com/JeffersonLab/hcana/blob/fcef8e23e16228fbb930078b9c44ef5031bed4f9/src/THcReactionPoint.cxx#L56).
Those files were not freshly retrieved in this audit.

## Saved reconstruction quantities versus optimization targets

The `T` suffix identifies constructed fit targets in the ntuple; it does not
certify a physical event origin. The following mappings are explicit in F2.
Lengths are cm unless a conversion is stated. Slopes are dimensionless,
conventionally reported using the small-angle approximation in radians.

| Quantity | Input or construction | Reference and use | Exact local source |
|---|---|---|---|
| Saved `xtar`, `ytar` | `H.gtr.x`, `H.gtr.y` | Final saved target-plane intercepts; distinct from `xtarT`, `ytarT` | `make_fit_ntuple_from_gmm.C:471–472` |
| Saved `xptar`, `yptar` | `H.gtr.th`, `H.gtr.ph` | Final saved TRANSPORT slopes | `make_fit_ntuple_from_gmm.C:477–478` |
| `reactx`, `reacty`, `reactz` | `H.react.x`, `H.react.y`, `H.react.z` | Saved reaction coordinates; the active reaction/beam-provider implementation remains to be identified | `make_fit_ntuple_from_gmm.C:473–475` |
| `xbpm_tar`, `ybpm_tar` | `H.rb.raster.fr_xbpm_tar`, `H.rb.raster.fr_ybpm_tar` | Raw raster/BPM branches; the HMS target construction uses `-xbpm_tar`, not `reactx` | `make_fit_ntuple_from_gmm.C:489–490`; `spectrometer_config.h:80–85` |
| `xs`, `ys` / local `xsieve`, `ysieve` | `H.extcor.xsieve`, `H.extcor.ysieve` | Saved replay sieve projections, distinct from assigned nominal hole centers | `make_fit_ntuple_from_gmm.C:483–484,511–514` |
| `xsT`, `ysT` | Centers selected by the event's assigned hole labels | Nominal analysis grid at the reconstruction sieve plane; not an independently known scattering origin | `make_fit_ntuple_from_gmm.C:599–600` |
| `ztarT` | Nominal foil coordinate `zf=info.zfoil[foilT]` | Laboratory longitudinal foil label; not HMS longitudinal coordinate `z_v` | `make_fit_ntuple_from_gmm.C:591,601` |
| Saved `ztar` | `reactz` | Reconstructed reaction coordinate; not the nominal foil label | `make_fit_ntuple_from_gmm.C:602` |
| `xtarT`, `ytarT`, `xptarT`, `yptarT` | `targetTruth` from nominal foil/grid, reaction coordinates and raw `xbpm_tar` | Construction B above, including the identified nonclosing terms | `make_fit_ntuple_from_gmm.C:592–595`; `spectrometer_config.h:80–85` |
| Fifth polynomial input | Saved `xtar/100`, together with focal-plane coordinates | The Python fit uses saved `xtar`, not constructed `xtarT`, including terms held fixed in this input | `fit_elastic.py:42–55` |
| Fit outputs | `xptarT`, `ytarT/100`, `yptarT` | Target y length is converted to meters for fitting | `fit_elastic.py:55`; `fit_opt_matrix_gmm.C:575–577` |

The existing geometry diagnostic also evaluates the unmodified target
construction and compares its x closure with `-reacty-xMis(angle)`
(`core_geometry.C:271–277`). Agreement with that identity tests the source
transcription; it does not establish physical vertex-to-hole closure.

Saved quantities are documented for provenance only; no saved events enter
this milestone. Construction C below receives prescribed analytic values at
the HCANA matrix interface. A full HCANA iteration comparison would require
the actual matrix, parameters, revisions and beam-provider configuration.

## Closure gate before claiming run-specific geometry

### Current deterministic test protocol

No sieve events are validation inputs. The model traces A's exact endpoint ray,
B's unchanged fit construction, and C's coordinate interface with synthetic
inputs. Each setting has its own declared translation model. One constant
sieve-plane registration offset is anchored to the central zero-displacement
case, then held fixed across holes, foils, and transverse displacements.
Raw and offset-subtracted results are both retained. Target-plane differences
and `atan` slope differences are retained because a constant residual at one
plane does not establish full-ray equivalence.

C supplies exact A values at the unevaluated matrix interface. The directly
readable [local mirror](../../sieve_afterburner.py), lines 47–57, distinguishes
the previous slope used for the intercept from the updated projection slope:

```
x0 = -Y - mx
xt_C = x0 - p_old * Z_lab * cos(theta)
xs_C = xt_C + L * p_new
ys_C = 100 * y_new_matrix_m + L * q_new
```

In the coordinate benchmark `p_old=p_new=p_A`, `q_new=q_A`, and
`y_new_matrix_m=yt_A/100`. The known synthetic lab Z replaces the saved
reaction Z at this interface; Z is not reconstructed. For transport-component
pointing at zero out-of-plane angle, `xs_C-xh=-X*p_A*sin(theta)`. The y agreement only checks units and
projection. An additional test injects unequal old/new slopes to retain the
term `Z*cos(theta)*(p_new-p_old)`. This is not a matrix iteration or a complete
HCANA/ReactionPoint test. The upstream provenance remains L1, with the local
mirror independently readable; no unavailable ReactionPoint equation is
invented.

### Evidence still required

- Recover the actual HCANA and analyzer revisions, replay parameter files and
  beam-provider configuration; identify the coordinate-conversion call path.
- Resolve whether active pointing components denote transport translations or
  a lab displacement, and whether they represent measured alignment or
  effective reconstruction corrections.
- Resolve beam/raster branch signs, units, slopes and reference planes against
  the lab vertex definition. Do not infer `B=X_lab` from the fit formula alone.
- Establish the physical plate front, back, normal, thickness, hole sizes and
  axes, exceptional holes and offsets from a drawing/survey. Keep those
  separate from the 168 cm reconstruction plane.
- Record foil material, dimensions, thickness, orientation and survey offsets;
  a nominal center-foil placement is not a complete material specification.
- Evaluate independent rigid-ray closure, the unchanged fit comparison and the
  replay-style projection separately. Report residuals and open conflicts;
  passing the first test must not be used to declare the others equivalent.

Until these gates are closed, any geometry is an explicitly labeled reference
model. It does not establish focal-plane transmission, detector cuts or the
yield of the observed inter-hole population.
