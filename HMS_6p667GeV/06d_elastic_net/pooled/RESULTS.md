# Reduced-basis seed study: pooled training validation

12/12 target/penalty cases completed every fold. Incomplete cases are excluded from pooled scores and plots.

Each balanced training event supplies exactly one out-of-fold residual per complete case. Every fold independently selects terms, learns scaling, and refits scaled X directly by unpenalized SVD on the other folds. No protected events are used.

Read plots/plateau.png first: pooled median absolute residual, mean foil/delta MSE, P90 absolute residual, and the range of fitting-matrix condition numbers. Faint points show separate folds; horizontal spans show their term-count range. The condition-number point is the median of fold condition numbers, not a condition number computed from pooled residuals. The constant is counted as a term but handled separately from the centered/scaled nonconstant columns used for conditioning.

Median and P90 pool individual events. cell_mse averages MSE equally over populated physical foil/delta cells, retaining the earlier scoring convention; mse in the tables is the ordinary event-weighted value. Residual units are mrad for angles and cm for ytar. Quantiles are computed from pooled residuals, never by averaging fold quantiles.

plots/foil_delta.png checks mean-squared residuals throughout acceptance. foil_delta.tsv also supplies median and P90 by cell. coverage.tsv gives pooled and individual-fold counts down to setting/hole, with qa_low as a reporting flag only. Pooling does not create extra events in sparse holes.

Use modest, consistent changes across nearby term counts to nominate several seeds. These four penalty settings can suggest a stable range, not establish a finely resolved plateau. There is no automatic winning model or requirement to beat full-basis SVD. Basis membership may differ between folds; terms.tsv and coefficients.npz preserve each actual seed. Do not treat the average term count as a single fitted basis. Beam search may reintroduce omitted terms and refit without penalties.

fold_scores.tsv preserves individual-fold metrics; pooled.tsv preserves exact pooled metrics. residuals.npz stores complete-case residual pairs in [full_svd, refit_svd] order with event IDs and folds. seed.dat, folds.tsv, code/, checkpoints.tsv and manifest.json make the study reproducible. No final replay matrix is exported.
