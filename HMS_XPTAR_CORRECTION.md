# HMS xptarT correction

The fit reference slope now traces the ray from the vertical interaction vertex
to the nominal sieve-hole center:

```text
xptarT = (xsieveT + reacty + x_mis) / (168 - zfoil*cos(CentAngle))
xtarT  = -reacty - x_mis - xptarT*zfoil*cos(CentAngle)
```

Positions are in cm and the angle in this expression is in radians. The existing
HMS longitudinal approximation and mispointing convention are retained.
The pair satisfies `xtarT + 168*xptarT = xsieveT` for center and displaced foils.
The fitted polynomial continues to use replay `xtar` consistently. SHMS retains
its own vertex and longitudinal expressions.

## Pulling this update

Rebuild TFit ntuples before fitting. Existing ROOT files, coefficients, residual
archives, plots, and frozen study code were not regenerated or rewritten.
The SVD and elastic-net/beam adapters reject old HMS ntuples rather than silently
fit the old reference slope. Residual plotting from TFit must use rebuilt files.

Where present, elastic net, pooled/refit/convergence studies, beam search, matrix
comparisons, and conditioning diagnostics consume the same rebuilt TFit values.
Huber expansion constructs corrected slopes directly from saved `reacty` and
requires a beam reference produced with the corrected definition. Equal-foil
Huber requires a corrected unweighted Huber reference. Run their dependencies in
that order in a fresh campaign/output series; historical directories are immutable.

The focal-plane 3D sheet geometry and physical labels do not need a new definition.
Its downstream ntuple builder uses the corrected shared HMS geometry when the
angular fit is run. There is no second vertex correction to add to its residuals.

In the current geometry diagnostic, the reference-ray contribution to the residual
slope is now zero. The former +5.952 mrad/cm contribution arose from the omitted
vertex term, not a prediction for the remaining reconstruction error.
