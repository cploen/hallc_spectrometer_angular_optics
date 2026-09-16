# HMS 6.667 — exact-membership fit-tree export

Conversation recap, 15 September 2026. Method/commands:
[CORE_SAMPLE.md](../../CORE_SAMPLE.md#6-build-compatible-hms-fit-trees) and
[ELASTIC_NET.md](../../docs/ELASTIC_NET.md).

Export `equal15` separately as `fit`, `holdout`, then `surplus`. Reuse the frozen
allocation metadata rather than current mutable run metadata (`50ecfe9`). The
builder joins requested replay entries, constructs the existing HMS target
quantities, and verifies that precisely those IDs appear in TFit. Filenames
retain `_fit_tree_gmm.root` for compatibility; the containing directory and
membership, not that suffix, identify the new core sample.

The farm checks reported **68,535 fit** and **101,258 holdout** exact-membership
events. The later surplus build reported **241,978 events**. “Output exists”
means that an immutable export directory exists; it does not prove a fit ran.
`run_build_core_fit.sh` never solves SVD. Validate/reuse the existing export;
do not delete it just to silence that message.

The original ROOT solver could otherwise impose its old admission caps again.
`preallocated_svd.py` provides an acceptance-only check bypassing them. Its
`--fit` mode is **not the modern direct-X reference used in this study**.
The elastic-net/beam Python path consumes all verified fit IDs and performs
the centered/scaled direct-X solves. No protected or surplus IDs enter those fits.

The fitting inputs retain saved replay xtar. The seed's xtar-dependent terms
and delta column are fixed; the angular/ytar constants are freely fitted.
Exported truth uses `spectrometer_config.h` via `make_fit_ntuple_from_gmm.C`.
Geometry questions documented later were **not corrected in these exports**.

Code: `build_core_fit.py`, `make_fit_ntuple_from_gmm.C`, `preallocated_svd.py`.
Tests: `tests/test_frozen_reallocation.py`, `tests/test_preallocated_svd.py`,
and `tests/test_elastic.py` cover export/admission and downstream separation.
ROOT-dependent fixture tests are identified in the balance verification note.

For exact continuation retain the complete `equal15` parent tag and
`06c_core_ntuple/equal15/{fit,holdout,surplus}/`, including `build.json`, frozen
metadata/seed, membership files and ROOT outputs. These exports are not tracked
in GitHub at this audit. Rebuilding also needs the original replay trees from
the campaign rungroup table. Plots alone cannot replace these inputs.
