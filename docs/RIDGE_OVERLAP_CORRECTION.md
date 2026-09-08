# Foil-envelope overlap correction

The overview PDF was drawn before the final overlap guard. Its envelopes
therefore differed from the saved TCutG polygons. Boundary and width graphs
also retained values from before protection.

The guard additionally paired neighboring foils by array index. Foils can
skip different delta rows, so equal indices need not describe equal delta.

In ytar_ridge_cut.C, ProtectRidgeEnvelopes now aligns piecewise-linear
envelopes on the union of retained delta knots, without extending any foil's
delta support. Every overlapping pair is protected in geometric ytar order,
including outer foils where an intermediate foil has no support. The existing
0.04 cm gap is retained. The common knots protect the intervening linear
polygon segments as well as the sampled points.

Inserted rows interpolate diagnostic quantities; nrow=-1 identifies them as
interpolated rather than measured rows. Existing measured row counts remain.
Crossed or insufficiently separated centers, or inverted envelopes, stop
processing before output rather than writing an invalid new cut.

Protection and polygon/graph/count updates now occur before every output.
CSV boundaries, final widths, plots and ROOT objects share the final geometry.
The minRowEntries default remains 40; this correction is separate from the
temporary 100-event threshold experiment.

Targeted test: tests/test_ridge_overlap.C exercises unequal grids, reversed
foil IDs, unequal support, intervening-segment separation, disjoint supports,
and rejection of crossed/too-close centers. Real replay plots must be
regenerated on ifarm. Existing output files are not automatically updated.
