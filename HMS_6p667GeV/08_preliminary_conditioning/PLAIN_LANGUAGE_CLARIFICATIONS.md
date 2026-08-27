# Plain-language clarifications for the 6.667 GeV angular study

This file collects the main clarifications from the preliminary angular-optics study. It is meant for collaboration discussion. It does not recommend which matrix should be used for production.

## 1. Matrix rows, fit terms, coefficients, and SVD modes

These words do not mean the same thing.

- A **matrix row** is one polynomial pattern, such as a particular set of powers of the focal-plane variables.
- A **fit term** is one polynomial pattern that is allowed to change in this angular fit.
- A **coefficient** is the number multiplying a fit term for one reconstructed quantity.
- An **SVD mode** is an independent combination of fit terms that the data may constrain well or poorly. A mode is usually not one individual term.

The seed matrix has 461 valid rows:

- 209 rows without explicit `xtar` dependence are refitted.
- 252 rows with explicit `xtar` dependence are carried through without being refitted.
- The fit also adds a constant term. This gives 210 fitted terms.

Each of the 210 fitted terms has a separate adjustable coefficient for `xptar`, `ytar`, and `yptar`. The angular fit therefore adjusts 630 numbers in total: 210 terms times 3 reconstructed quantities.

The SVD has at most 210 modes because the design matrix has 210 fitted-term columns. The same 210-column basis is used for each of the three reconstructed quantities. The 461 matrix rows therefore do **not** mean 461 fitted modes.

The matrix file also contains a delta coefficient column. This angular fit does not adjust that column.

## 2. What one point on the term ladder means

Every value of `N` is a separate fit.

For example, the `N=5` point fits the data using only the first five allowed terms. It is not the result of fitting all 210 terms and then inspecting the first five coefficients.

When `N` increases, every allowed coefficient is solved again. Earlier coefficients may change because the new terms give the fit more freedom.

The ladder ran quickly because the events and polynomial values were loaded once. The full `X^T X` matrix and `X^T y` vectors were also formed once. For each `N`, the code selected the leading `N` rows and columns and performed a new solve. This is a valid nested refit.

One limit is important: the ladder sums were formed with NumPy, while the production ROOT macro forms the sums in its event loop. Both represent the same equations, but floating-point rounding may differ slightly. This matters when the equations are extremely poorly conditioned.

## 3. Complete polynomial degrees

With four focal-plane input variables, the complete polynomial endpoints are:

| N | Included polynomial set |
|---:|---|
| 5 | Complete through linear order |
| 15 | Complete through quadratic order |
| 35 | Complete through cubic order |
| 70 | Complete through quartic order |
| 126 | Complete through fifth order |
| 210 | Complete through sixth order |

`N=35` is the last complete cubic fit. It contains every term through cubic order.

`N=39` contains all terms through cubic order plus only the first four quartic terms in the file order. It is an incomplete quartic prefix.

An incomplete order is allowed, but its result depends on which terms happen to come first. It should not be described as a general quartic fit. This matters if we compare `N=39` with another matrix that uses a different term order.

## 4. The two main ways of solving the fit

Let `X` be the table of polynomial values for the selected events, `c` the coefficient vector, and `y` the target values.

The direct problem is:

`X c = y`

The current ROOT method first forms the normal equations:

`(X^T X) c = X^T y`

`X^T X` is also called the Gram matrix. In these notes, **production Gram** means the unscaled `X^T X` matrix used by the current fitting method. It does not mean that the matrix has been approved for production.

The **right-hand side** is `X^T y`. There is one right-hand-side vector for each output: `xptar`, `ytar`, and `yptar`.

To **accumulate** `X^T X` and `X^T y` means to add the contribution from every training event:

- A Gram-matrix value is a sum of products of two term values.
- A right-hand-side value is a sum of a term value times the target value.

Changing the order or software used for these sums changes the last floating-point digits. A stable solve should barely respond. A poorly conditioned solve can respond much more.

## 5. Scaling the columns

Different polynomial columns of `X` can have very different numerical sizes. Column scaling changes the units of those columns before solving so that they have comparable sizes. The coefficients are converted back to the original units afterward.

Scaling does not add data, remove data, or change the polynomial model. It makes the numerical problem easier for the solver.

For the full 210-term fit:

| Method | Approximate condition number | Modes retained |
|---|---:|---:|
| Unscaled `X^T X`, current ROOT method | `5.7 x 10^26` | 111 of 210 |
| Direct solve on unscaled `X` | `3.5 x 10^13` | 210 of 210 |
| Scaled `X^T X` | `1.2 x 10^11` | 210 of 210 |
| Direct solve on scaled `X` | `3.5 x 10^5` | 210 of 210 |

Solving directly on unscaled `X` improves the condition number by about 13 powers of ten compared with solving unscaled `X^T X`. Scaling `X` improves it by about another eight powers of ten.

In this dataset, scaling also allowed the `X^T X` solve to retain all 210 modes. The direct scaled-`X` solve was still much better conditioned.

## 6. ROOT's truncation threshold

SVD separates the fit into modes. ROOT compares each mode strength with a numerical cutoff. Modes below that cutoff are set to zero in the pseudoinverse.

The cutoff is recalculated for every solve because it depends on the strongest mode in that solve. There is no single fixed cutoff line for the whole term ladder.

The plots therefore show the distance between the weakest mode and that solve's own cutoff. A margin of zero means the weakest mode is at the cutoff.

A **decade** means a factor of ten:

- `+1 decade` means 10 times above the cutoff.
- `+2 decades` means 100 times above the cutoff.
- `-1 decade` means 10 times below the cutoff.

For the current unscaled `X^T X` solve:

- Through `N=39`, every mode is retained with a useful margin above the cutoff.
- From `N=40` through `N=53`, every mode is still retained, but the weakest mode is only about 1.36 to 1.48 times above the cutoff.
- At `N=54`, the first mode is discarded.
- At `N=210`, ROOT retains 111 modes and discards 99 modes.

The first discarded mode at `N=54` does **not** mean that every later coefficient is suppressed. As more terms are added, the retained rank continues to rise. It reaches 111 at `N=176` and stays at 111 through `N=210`.

ROOT suppresses weak **combinations of terms**, not necessarily the newly added term or one named coefficient. We therefore cannot say, “all coefficients after term 54 are suppressed.”

No mode was discarded in either direct-`X` ladder through `N=210` under the tested NumPy/LAPACK cutoff. The numerical rank reported by that solve was 210. This does not mean all coefficient directions are equally precise; the unscaled direct solve still has a large condition number.

## 7. What “rank deficient” means here

The current ROOT `X^T X` solve is numerically rank deficient under ROOT's cutoff: it uses only 111 of the 210 possible modes at `N=210`.

This is not proof that the original design matrix `X` has an exact algebraic rank of only 111. Direct solves on `X` retained all 210 modes. Much of the loss comes from forming `X^T X`, which makes the numerical conditioning much worse.

The practical risk is that the data do not determine some combinations of coefficients reliably under the current calculation. Small changes in rounding, parsing, or solver method can move individual coefficients while leaving the training-event predictions nearly unchanged. Predictions away from the strongly sampled calibration region can change more.

The phrase “not an unconstrained 210-parameter determination” means that ROOT did not independently determine all 210 coefficient combinations for each output. It set 99 weak combinations to zero according to its cutoff. The result is therefore partly defined by the cutoff rule.

In plain language: the high-order current solution is not determined only by the calibration data. It is also shaped by which weak directions ROOT decides are too small to use.

## 8. Would the scaled solution be very different?

There are two separate questions:

1. Are the fitted coefficients similar?
2. Are the reconstructed predictions similar on training and held-out events?

Poor conditioning can make coefficient vectors look quite different even when predictions on the training sample are close. A better-conditioned method can therefore give a different matrix without a large training-residual change.

For this study, the full scaled direct-`X` fit had slightly smaller training residuals than the current unscaled `X^T X` fit:

| Output | Current `X^T X` training RMS | Scaled direct-`X` training RMS |
|---|---:|---:|
| `xptar` | 1.978 mrad | 1.906 mrad |
| `ytar` | 0.1478 cm | 0.1449 cm |
| `yptar` | 1.126 mrad | 1.117 mrad |

“Slightly smaller” applies only to the calibration events used in the fit. It does not mean the full scaled fit performed better overall.

On the eligible events not selected for fitting, the full scaled 210-term fit was worse:

| Output | Current `X^T X` held-out RMS | Scaled direct-`X` held-out RMS |
|---|---:|---:|
| `xptar` | 1.874 mrad | 3.049 mrad |
| `ytar` | 0.1455 cm | 0.5141 cm |
| `yptar` | 1.069 mrad | 3.838 mrad |

This is overfit-like behavior: the added freedom improves the training sample slightly but gives poorer results on held-out events. Much of the held-out change is associated with the `+8 cm` foil sample. Scaling fixes a numerical problem; it does not by itself prevent overfitting.

## 9. Training and held-out samples

The study used:

- 166,013 selected training events.
- 278,684 eligible events that were not selected for the fit.

The held-out events come from the same runs and cuts. They are useful for detecting poor generalization inside this dataset, but they are not a fully independent run-level validation sample.

The residual ladders were made separately for each foil and for both solve methods:

- Current unscaled `X^T X` solve.
- Scaled direct-`X` solve.

Each `N` point uses coefficients refitted at that `N`. The plots do not take one 210-term fit and simply hide later terms.

## 10. Parser correction

The supplied `oldfit.dat` contains a blank line after a separator. The old parser did not check whether all nine expected values were read before it stored a row.

As a result, the prior ifarm run included one extra machine-dependent fixed-`xtar` row: 253 fixed rows instead of 252. The malformed row was not one of the 210 adjustable terms. It still changed the fixed contribution subtracted from the targets, so it could change the fitted angular coefficients.

The strict parser stores a row only when all nine expected values were read.

The corrected valid coefficient vectors differ from the prior ifarm output by 2.3% to 3.0% in relative L2 norm. This is the size of the whole-vector difference; it does not mean every coefficient changed by that percentage.

The training-event prediction changes between those two matrices were small:

- `xptar`: 0.00455 mrad RMS.
- `ytar`: 0.000781 cm RMS.
- `yptar`: 0.00264 mrad RMS.

The coefficient change, together with the cutoff and summation sensitivity, is evidence that individual coefficients are numerically sensitive. The small training-prediction change shows that several coefficient combinations can describe the calibration sample almost equally well.

## 11. What “GMM-clean matrix” means

The prior ifarm file is named `nps_hms_newfit_6p667_gmm_clean.dat`.

Here, **clean** refers to the GMM-cleaned event sample. It does not mean that the parser was correct, that the numerical solve was stable, or that the matrix was approved for production.

To avoid confusion, this study calls it the **prior ifarm GMM-clean matrix**. The corrected local output is called the **strict-parser local rerun**.

## 12. Leave-one-setting-out calculation

A leave-one-setting-out condition-number table was produced during the preliminary work. It removed one setting at a time from `X` and recomputed the condition number. It did not refit and validate a new optics matrix for each removal.

This was not requested before it was run. It should not be used as support for removing any setting, and no setting was removed from the main fits or residual plots. All seven settings were used for the main results.

If setting-removal studies are ever wanted, they need a separate approved plan with complete refits and training and held-out checks.

## 13. What the preliminary results do and do not establish

The results support these limited statements:

- The current unscaled `X^T X` method becomes close to ROOT's cutoff at `N=40` and first discards a mode at `N=54`.
- At `N=210`, that method retains 111 of 210 modes.
- Directly solving on `X` avoids the large extra conditioning penalty caused by forming `X^T X`.
- Scaling the columns improves the numerical conditioning further.
- The full scaled 210-term fit has lower training residuals but much worse held-out residuals in this preliminary split.
- The parser correction is necessary for a reproducible input matrix.

The results do not yet establish:

- Which term count should be used for production.
- Which solve method should be used for production.
- That retaining all 210 modes gives better physics results.
- That one data setting should be removed.
- Performance on independent runs or all production kinematics.

The safe plain-language conclusion is: the current high-order solve is numerically sensitive, while the full scaled solve shows overfit-like held-out behavior. Numerical stability and held-out physics performance must both be checked. Improving one does not guarantee the other.
