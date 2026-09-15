# HMS target and replay geometry: an informal consistency proof

**Review note — 14 September 2026. No reconstruction or fit code has been changed.**

This note asks a narrow question: **Do the target positions and slopes supplied to the HMS matrix fit describe a ray that passes through the nominal sieve hole, in the coordinate convention used to display the replay?**

It builds the geometry first, then substitutes the implemented equations. There are two algebraic target-equation mismatches, followed by two separate approximation issues in the fit/replay chain. None of these findings establishes the magnitude of an error in the actual 6.667 data without checking the historical replay configuration.

Repository inspected: `cploen/hallc_spectrometer_angular_optics`, commit `b6b1fcf36b62d3f40d43b3ae59f7d5af70145035`. HCANA reference inspected: commit `fcef8e23e16228fbb930078b9c44ef5031bed4f9`. The latter is a pinned upstream reference, **not a verified identification of the HCANA binary used to create the campaign files**.

## 1. Coordinate systems and definitions

A particle starts at a beam–foil intersection and travels through a sieve hole before entering the spectrometer magnets. For the reference geometry, consider the straight line joining that vertex to the hole. Scattering and finite hole size produce a distribution around this reference ray; they do not change its defining line equations.

The laboratory and spectrometer coordinate systems use different axis names:

- Laboratory **Y** is vertically upward; laboratory **Z** follows the beam downstream. Laboratory **X** is the horizontal transverse direction, conventionally toward beam left in the HCANA description.
- Spectrometer TRANSPORT **x** is vertically downward; **z** follows the central spectrometer ray toward the magnets; **y** completes the right-handed system.
- Thus a laboratory vertical position called `react.y` belongs in the spectrometer **x** equations. Its appearance there is not an accidental exchange of axes.
- The target reference plane is the plane **z = 0 in spectrometer coordinates**. It is not the physical location of every target foil.

See the [HMS coordinate description, section II A of Puckett et al.](https://misportal.jlab.org/sti/publications/14949/attachments/440/1707.07750.pdf#page=2) and the [HCANA coordinate comments](https://github.com/JeffersonLab/hcana/blob/fcef8e23e16228fbb930078b9c44ef5031bed4f9/src/THcHallCSpectrometer.cxx#L55). Beam-provider signs and actual spectrometer-angle settings still require confirmation for a particular replay.

All lengths below are **cm**. Slopes are dimensionless, conventionally reported in radians or mrad under the small-angle approximation. A slope is dx/dz or dy/dz; the exact geometric angle would be its arctangent.

| Symbol | Definition and connection to the code |
|---|---|
| `R_y` | Laboratory vertical reaction position, read from `H.react.y`. |
| `R_x` | Laboratory horizontal reaction position, read from `H.react.x`. The current HMS target routine receives this value but does not use it in its HMS branch. |
| `z_f` | Assigned nominal/surveyed foil position along the laboratory beam direction; stored as `ztarT`. |
| `z_r` | Reaction z reconstructed by the original replay; read from `H.react.z`. This is distinct from `z_f`. |
| `theta` | Positive HMS angle magnitude used by the campaign target routine; here 12.490 degrees in the numerical example. |
| `s`, `c` | `sin(theta)` and `cos(theta)`, with degrees converted to radians before evaluation. |
| `m_x`, `m_y` | Spectrometer x and y mispointing translations, in cm. The target routine uses the repository's angle-dependent HMS formulas. Replay may instead have explicit parameter values. |
| `x_v` | Vertical vertex coordinate relative to the spectrometer axis: `x_v = -R_y - m_x`, following the implemented convention. |
| `B` | The horizontal beam quantity used by the HMS target code: **`B = -xbpm_tar`**, where `xbpm_tar` is read from `H.rb.raster.fr_xbpm_tar`. Whether this equals the laboratory horizontal position used by the replay must be verified; it is not assumed here. |
| `y_v` | Horizontal vertex coordinate implied by the target routine: `y_v = z_f s + B c - m_y`. |
| `z_0` | Approximate spectrometer-longitudinal vertex coordinate `z_f c`, omitting the horizontal beam-position contribution. |
| `z_v` | Spectrometer-longitudinal vertex coordinate of a fully specified reference ray. A common value must be used in both its x and y projections. |
| `L` | Spectrometer-longitudinal position of the sieve plane. Both the current HMS profile and referenced replay use **168 cm**. This note does not independently verify the surveyed distance. |
| `x_h`, `y_h` | Nominal sieve-hole center coordinates relative to the sieve coordinate origin, supplied by the fixed hole grid. They are stored as `xsT`, `ysT`. They are not event-by-event reconstructed hit positions. |
| `x_t`, `y_t` | Track intercepts at the spectrometer target reference plane, z = 0; `xtar`, `ytar` in geometric units. Neither is generally equal to the vertex coordinate at a displaced foil. |
| `p`, `q` | Track slopes: `p = xptar = dx/dz`, `q = yptar = dy/dz`. |
| subscript `T` | A target value constructed for fitting; for example `p_T = xptarT`. “True” in a branch name does not guarantee that its defining geometry is correct. |
| subscript `R` | A value saved by the original replay, for example `x_R = saved xtar`. |
| `u` | Four measured focal-plane variables: xfp, xpfp, yfp, ypfp. |
| `F_M(u,x_t)` | The polynomial target reconstruction using matrix M and the supplied xtar input. Its xptar component is written `F_M^x`. |

**Unit caution:** transport matrices use length variables in meters, while the event tree and the geometry above use cm. The Python fit converts xfp, yfp and xtar by dividing by 100; matrix ytar output is converted back by multiplying by 100. The mismatches derived below persist after those conversions and are not a missing factor of 100.

## 2. Establish the ray identities before looking at code

Project the three-dimensional ray into the spectrometer x–z plane. The two endpoints are:

- vertex V: `(z_v, x_v)`;
- hole H: `(L, x_h)`.

The horizontal run in this drawing is `L - z_v`; the vertical rise is `x_h - x_v`. The extension of the line to z = 0 determines the intercept x_t. If the foil is downstream of z = 0, that extension runs backward from the vertex; if upstream, the reference plane is between vertex and sieve.

Therefore,

$$p=\frac{x_h-x_v}{L-z_v}. \tag{1}$$

The same line must satisfy

$$x_t=x_v-pz_v. \tag{2}$$

Projecting it to the sieve gives

$$x_t+Lp=x_v+p(L-z_v)=x_h. \tag{3}$$

Exactly the same argument in the y–z projection gives

$$q=\frac{y_h-y_v}{L-z_v},\qquad y_t=y_v-qz_v,\qquad y_t+Lq=y_h. \tag{4}$$

These are line identities, not a proposed fitting method. Any slope/intercept pair claimed to describe this reference ray must obey them in one common coordinate system. No matrix, SVD or regularization enters the proof.

## 3. First finding: the vertical slope and intercept use different vertex assumptions

The current [HMS target routine](https://github.com/cploen/hallc_spectrometer_angular_optics/blob/b6b1fcf36b62d3f40d43b3ae59f7d5af70145035/spectrometer_config.h#L78) uses

$$p_T=\frac{x_h}{L-z_0}, \tag{5}$$

$$x_{t,T}=x_v-p_Tz_0,\qquad x_v=-R_y-m_x. \tag{6}$$

Equation (5) treats the vertical vertex displacement as zero. Equation (6) retains it. To determine whether these can nevertheless describe the same ray to the nominal hole, substitute both into the sieve identity:

$$
\begin{aligned}
x_{t,T}+Lp_T
 &=x_v-p_Tz_0+Lp_T\\
 &=x_v+p_T(L-z_0)\\
 &=x_v+x_h.
\end{aligned} \tag{7}
$$

The mismatch is therefore

$$\boxed{\epsilon_x\equiv x_{t,T}+Lp_T-x_h=-R_y-m_x.} \tag{8}$$

**The contradiction:** if x_h is the absolute nominal hole position in the same coordinates used by the projection, and x_v is nonzero, equations (5) and (6) cannot both describe a ray through x_h. They describe a ray through **x_h + x_v** instead.

This is about **both the missing `react.y` contribution and missing x mispointing in the slope numerator**. They are included in xtarT, so it is not accurate to say they are absent from the whole calculation.

Under the same approximate longitudinal geometry, a consistent slope would be

$$p_{\mathrm{ray}}=\frac{x_h-x_v}{L-z_0}=\frac{x_h+R_y+m_x}{L-z_0}. \tag{9}$$

Equation (9) is the consequence of the stated assumptions, **not an instruction to patch the fit immediately**. A different historical coordinate origin or an intended intermediate, beam-relative slope would have to be documented and followed by the corresponding conversion to the final physical slope.

### Does the code actually pass a fixed nominal hole coordinate?

Yes. The [builder creates the hole grid](https://github.com/cploen/hallc_spectrometer_angular_optics/blob/b6b1fcf36b62d3f40d43b3ae59f7d5af70145035/make_fit_ntuple_from_gmm.C#L534), passes the grid coordinate to `targetTruth`, and then stores that same grid coordinate as `xsT`. There is no event-dependent subtraction of x_v at that call site. See the [target construction and storage](https://github.com/cploen/hallc_spectrometer_angular_optics/blob/b6b1fcf36b62d3f40d43b3ae59f7d5af70145035/make_fit_ntuple_from_gmm.C#L591).

That strengthens the coordinate concern. It does not prove that every historical matrix was intended to use these target values as final physical angles; that interpretation still needs checking.

### Numerical counterexample

Choose the central foil and central hole: z_f = 0, x_h = 0. Let R_y = 0. The current HMS mispointing formula at 12.490 degrees gives

$$m_x=0.1[2.37-0.086(12.490)+0.0012(12.490)^2]=0.148306\ \mathrm{cm}.$$

The implemented targets are p_T = 0 and x_t,T = −0.148306 cm. They project to −0.148306 cm at the sieve, not to the hole at zero.

The reference ray to zero instead has

$$p_{\mathrm{ray}}=0.148306/168=0.000882774\approx0.883\ \mathrm{mrad}.$$

This is a counterexample to the simultaneous definitions under the stated conventions. **It is not a measurement of a 0.883 mrad error in this campaign.** Beam position, actual mispointing settings, fitted offsets and replay xtar can change the observed result.

## 4. Second finding: the horizontal target pair mixes orders in the beam-position correction

Use the code's definition B = −xbpm_tar without assuming its laboratory sign is verified. Write

$$y_v=z_fs+Bc-m_y,\qquad z_0=z_fc.$$

The current HMS equations are

$$q_T=\frac{y_h-y_v}{L-z_0}, \tag{10}$$

$$y_{t,T}=z_f(s-q_Tc)+B(c+q_Ts)-m_y. \tag{11}$$

Expand equation (11):

$$y_{t,T}=y_v-q_Tz_0+Bq_Ts. \tag{12}$$

Now project to the sieve:

$$
\begin{aligned}
y_{t,T}+Lq_T
 &=y_v+q_T(L-z_0)+Bq_Ts\\
 &=y_h+Bq_Ts.
\end{aligned} \tag{13}
$$

Thus

$$\boxed{\epsilon_y\equiv y_{t,T}+Lq_T-y_h=Bq_T\sin\theta.} \tag{14}$$

It vanishes if B, q_T or sin(theta) is zero. Otherwise the pair does not exactly project to y_h.

The source of the inconsistency becomes visible by rewriting equation (12) as

$$y_{t,T}=y_v-q_T(z_0-Bs).$$

That intercept uses the longitudinal vertex coordinate **z_v = z_0 − B s**, while the slope denominator in equation (10) uses **z_0**. Consistency would require either retaining that longitudinal correction in both equations or omitting it in both, with the approximation documented. With z_v = z_0 − B s, the ray slope would be

$$q_{\mathrm{ray}}=\frac{y_h-y_v}{L-z_0+Bs}. \tag{15}$$

This is separate from the missing laboratory vertical position in the x-angle numerator. It is a horizontal beam-position/longitudinal-distance term and can be much smaller.

For illustration only, B = 0.1 cm, q_T = 0.03 and theta = 12.490 degrees give epsilon_y ≈ 0.000649 cm = **0.00649 mm**. These are example inputs, not measured beam conditions. The size scales with B and q_T.

A unified three-dimensional derivation should use one z_v for both projections. The x proof above deliberately used z_0 throughout so that its larger numerator mismatch does not depend on resolving this smaller longitudinal correction.

## 5. Why a free constant does not prove the geometry is repaired

For one event, equation (9) differs from the current target by

$$p_{\mathrm{ray}}-p_T=-\frac{x_v}{L-z_0}.$$

A constant coefficient could supply that correction if it were the same for all relevant events. In general it changes with beam position and foil location.

There is a more basic point: the optimizer is given p_T as its desired answer. A zero residual against p_T is rewarded whether or not p_T describes the intended physical ray. A free constant can remove an average fit residual without correcting the target definition. Other polynomial terms may partly compensate, but that does not establish the correctness of the supplied geometry.

There are also two different meanings of “offset” that must remain separate:

- m_x and m_y are **position translations**, in cm, used in geometry.
- A matrix constant or external hphi_offset/htheta_offset is a **slope addition**, in radians, except for the ytar constant, which has length units.

An external angular correction and an active constant coefficient both add to the polynomial output. They must not be counted twice. The reference [HCANA matrix evaluator](https://github.com/JeffersonLab/hcana/blob/fcef8e23e16228fbb930078b9c44ef5031bed4f9/src/THcHallCSpectrometer.cxx#L573) adds hphi_offset to xptar and htheta_offset to yptar. That naming is counterintuitive but independent of the two projection identities above.

## 6. What the replay correction does—and what it does not establish

The inspected [HCANA extended-target correction](https://github.com/JeffersonLab/hcana/blob/fcef8e23e16228fbb930078b9c44ef5031bed4f9/src/THcExtTarCor.cxx#L112) constructs x_tg = −vertex.y − pointingOffset.x, updates xtar using the reconstructed vertex z and current xptar, reevaluates the matrix, and projects to xsieve = xtar + 168 xptar and ysieve = 100 ytar + 168 yptar. The factor 100 converts matrix ytar from meters to cm.

Thus the referenced replay explicitly includes the vertical vertex translation when constructing xtar. The missing-vertex concern identified here is in the fit's x-angle target equation, not an omission of vertex.y from this replay correction.

However, matching this upstream implementation is not proof of matching the historical job. We must check the deployed HCANA revision, central-angle conventions, active mispointing values, beam-provider coordinates, external offsets, and any focal-plane transformations. The reference [reaction-point code](https://github.com/JeffersonLab/hcana/blob/fcef8e23e16228fbb930078b9c44ef5031bed4f9/src/THcReactionPoint.cxx#L56) obtains horizontal and vertical coordinates from its configured beam object and calculates z from the track. This does not alone establish the relationship between that beam object's horizontal coordinate and `-H.rb.raster.fr_xbpm_tar`.

## 7. A separate limitation: the fit uses saved replay xtar, not xtarT

The [fit-ntuple builder](https://github.com/cploen/hallc_spectrometer_angular_optics/blob/b6b1fcf36b62d3f40d43b3ae59f7d5af70145035/make_fit_ntuple_from_gmm.C#L471) reads xtar from H.gtr.x and stores it separately from xtarT. The [modern fit](https://github.com/cploen/hallc_spectrometer_angular_optics/blob/b6b1fcf36b62d3f40d43b3ae59f7d5af70145035/fit_elastic.py#L42) evaluates the fixed xtar-dependent terms using that saved xtar.

Consequently, a perfect fit residual means

$$F_M^x(u,x_R)=p_T.$$

It does **not** mean that the fitted model has solved the joint geometry equations for a new x_t. If the perfect target slope were projected using the saved x_R, equation (8) gives

$$x_R+Lp_T-x_h=(x_R-x_{t,T})+x_v. \tag{16}$$

The two terms could add or partially cancel. This is why equation (8) alone cannot predict the actual sieve displacement of a fitted matrix in the frozen-input diagnostic.

A new replay generally supplies a different xtar input after evaluating the new matrix and reaction geometry. For a small input change, the corresponding prediction change is approximately

$$\Delta p\approx\frac{\partial F_M^x}{\partial x_t}\,\Delta x_t.$$

No value of that derivative is assumed here. This is a **fit-to-replay consistency limitation**, not a second proof that the x target algebra is wrong. Retaining saved inputs was an explicit choice that made the modern matrix comparisons controlled; it limits what those comparisons establish about a full new replay.

## 8. A further approximation: replay stops the xtar iteration at finite tolerance

The referenced extended-target routine stops when the xptar change is at most **2 mrad**, or after five updates. In an update, it constructs an intercept from the previous slope p_old and then calculates p_new from that intercept:

$$x_{t,\mathrm{stored}}=x_v-p_{\mathrm{old}}z_r c.$$

Relative to an exactly self-consistent pair using p_new, the remaining vertex-line mismatch is

$$x_{t,\mathrm{stored}}+p_{\mathrm{new}}z_r c-x_v=z_r c(p_{\mathrm{new}}-p_{\mathrm{old}}). \tag{17}$$

If convergence stops at a 2 mrad change and |z_r| = 8 cm, this is at most approximately **0.16 mm**. The bound does not apply when the five-update limit is reached while the change is still larger. Actual residual changes may be far smaller than the stopping threshold.

This is an identifiable numerical approximation, **not evidence that all replay events violate geometry by 0.16 mm**. The afterburner currently mirrors this stopping rule. A comparison with a tighter iteration tolerance would isolate its effect without changing the matrix or target definitions. It still would not address all differences between a saved-vertex afterburner and a complete replay.

## 9. Why an apparently good ztar diagnostic cannot settle these questions

The current [ztar residual calculation](https://github.com/cploen/hallc_spectrometer_angular_optics/blob/b6b1fcf36b62d3f40d43b3ae59f7d5af70145035/elastic_diagnostics.py#L12) recovers the horizontal beam quantity algebraically from the exported ytarT, yptarT and ztarT. It then uses the inverse of that same target relation to calculate reconstructed z.

For supplied predictions exactly equal to those targets, it necessarily returns z_f, provided the denominators are nonsingular. That is a useful internal consistency property. It is **not an independent verification of the beam-coordinate convention or the nominal hole projection**. In particular, the horizontal targets can satisfy that z equation while equation (14) still gives a nonzero sieve mismatch.

This does not invalidate comparisons of candidate matrices under the same definitions. It explains why those comparisons cannot by themselves certify the physical definitions.

## 10. What is proved, and what remains open

| Finding | What is established | What is not established |
|---|---|---|
| Vertical target mismatch | The implemented targets give x_t,T + L p_T − x_h = −R_y − m_x exactly. | Whether an intended historical coordinate conversion elsewhere compensates it; the actual event-level error or resolution impact. |
| Horizontal target mismatch | The implemented targets give y_t,T + L q_T − y_h = B q_T sin(theta) exactly. | Whether this retained/omitted term is a deliberate acceptable approximation; its size on actual data; the beam-provider sign mapping. |
| Saved xtar in fitting | The fit uses x_R, whereas a new reconstruction can produce a different xtar input. | Whether the resulting prediction change is material for this campaign. |
| Finite replay iteration | The source uses previous and updated slopes with a finite stopping tolerance. | The actual residual iteration error on the selected events. |
| ztar self-consistency | Its beam recovery and inverse use the same target relation. | Independent validation of that geometry or of a full replay. |

**Confidence:** high in the algebra and identification of the inspected code paths; conditional in the physical interpretation until coordinate origins and the historical replay settings are reconciled. No claim is made that correcting these items will necessarily improve a measured resolution.

## 11. A small, decisive next check

Before modifying any matrix or regenerating targets:

1. Establish a short convention sheet from the actual replay: laboratory/transport signs, foil positions, sieve origin and distance, active mispointing and angular offsets, and the relationship between the beam object's coordinates and the exported BPM/react branches.
2. Evaluate the line identities for a central-hole ray at z_f = 0 and at the outer foils. Repeat with a deliberate nonzero vertical displacement and then a horizontal displacement. A nominal central-hole, zero-displacement-only test can miss the issue.
3. On existing exported events where the necessary quantities are available, calculate x_t,T + L p_T − xsT and y_t,T + L q_T − ysT. Compare against equations (8) and (14). Do not substitute saved replay xtar for xtarT in this particular test.
4. Separately reproduce the old replay sieve coordinates using its original matrix/settings, then compare frozen-xtar evaluation with the iterative afterburner. This measures a different effect; it is not a substitute for step 3.
5. Only after the convention sheet identifies the intended reference ray should corrected targets be considered. Any fit comparison must then state that its target definitions changed; do not mix old-target and corrected-target residuals as though they shared the same truth.

## External references

The links below are useful for an independent discussion; none identifies the historical ifarm binary by itself.

- [Puckett et al., Technical Supplement, sections II A and II D](https://misportal.jlab.org/sti/publications/14949/attachments/440/1707.07750.pdf): transport coordinates and HMS optical reconstruction context. The line-identity proofs in this note are direct derivations, not claims that the paper endorses a code correction.
- [Pinned HCANA THcExtTarCor](https://github.com/JeffersonLab/hcana/blob/fcef8e23e16228fbb930078b9c44ef5031bed4f9/src/THcExtTarCor.cxx): extended-target iteration and sieve projection.
- [Pinned HCANA THcReactionPoint](https://github.com/JeffersonLab/hcana/blob/fcef8e23e16228fbb930078b9c44ef5031bed4f9/src/THcReactionPoint.cxx): beam-object coordinates and reaction-z calculation.
- [Pinned HCANA THcHallCSpectrometer](https://github.com/JeffersonLab/hcana/blob/fcef8e23e16228fbb930078b9c44ef5031bed4f9/src/THcHallCSpectrometer.cxx): coordinate comments, mispointing defaults, polynomial evaluation and external offsets.
- [Pinned campaign target routine](https://github.com/cploen/hallc_spectrometer_angular_optics/blob/b6b1fcf36b62d3f40d43b3ae59f7d5af70145035/spectrometer_config.h#L62): equations under review.
- [Pinned fit-ntuple builder](https://github.com/cploen/hallc_spectrometer_angular_optics/blob/b6b1fcf36b62d3f40d43b3ae59f7d5af70145035/make_fit_ntuple_from_gmm.C): source branches, nominal hole grid and target storage.
