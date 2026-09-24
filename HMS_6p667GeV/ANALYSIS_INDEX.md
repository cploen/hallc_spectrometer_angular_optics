# HMS 6.667 — analysis index for this Codex conversation

Audit: 15 September 2026, through code commit `2634d39`, branch
`core-sample-campaign`. Scope: core selection, historical sample audit, foil
balancing, elastic net, beam search, comparison and geometry discussion.
Earlier labeling/GMM development is context, not a new task in this record.

Paths and commands below are relative to the **repository root**. Set once:

```bash
C=HMS_6p667GeV
B="$C/06b_svd_fit/diagnostics/balance"
python3 -m pip install -r requirements_core.txt
```

These are reproduction commands, **not instructions to rerun completed work**.
Reuse verified existing tags; most production commands refuse overwrite.
For a changed experiment use new tags and update dependent source names. Exact
historical reproduction also requires the saved code/configuration and inputs.

## Sequence and documentation inventory

Latest addition (24 September): [smooth Huber and sample expansion results](../docs/HUBER_RESULTS_20260924.md)
and [reproduction](../docs/HUBER_EXPANSION.md). Actual frozen-data fits compare
68,535 balanced cores with up to 396,320 supported development events. Pooled
gains do not establish improvement across the foil/delta acceptance.

Follow-up: [equal-foil weighting results](../docs/HUBER_FOIL_RESULTS_20260924.md)
and [reproduction](../docs/HUBER_FOIL_WEIGHTS.md). Four fits retain all eligible
cores or cores plus shoulders with 20% base weight per foil. Core-only gains
are small; shoulder-inclusive fits retain an angular core/shoulder tradeoff.

1. **Select local cores; recover sparse edge holes (`core` → `min10`).**
   [Design/commands](../CORE_SAMPLE.md) · [decisions/QA](05c_core_sample/CORE_SAMPLE_RECAP.md) ·
   [code](../core_sample.py).
   `./run_core_sample.sh "$C" min10 --check`, then
   `./run_core_sample.sh "$C" min10`.
   Original `core` used min_events=20; current campaign configuration uses 10.

2. **Inspect hole counts; choose 1.5×P25 proposal and log-circle maps.**
   [Definitions](../CORE_SAMPLE.md#hole-count-histograms) ·
   [histogram code](../plot_core_counts.py) · [map code](../plot_core_map.py).
   `python3 plot_core_counts.py "$C/05c_core_sample/min10" rg01`;
   `python3 plot_core_map.py "$C/05c_core_sample/min10" rg01 --delta 0 1 2 3 4 --cap 1.5 --log`.

3. **Reconstruct actual historical SVD membership and three-stage rg01 clouds.**
   [Counts/provenance](06b_svd_fit/diagnostics/balance/HMS_6p667_SVD_SIEVE_BALANCE.md) ·
   [clouds](06b_svd_fit/diagnostics/balance/HMS_6p667_RG01_EVENT_CLOUDS.md) ·
   [extraction code](../diagnostics/conditioning/svd_counts.py) · [cloud code](../plot_svd_clouds.py).
   Commands:
   ```bash
   python3 diagnostics/conditioning/svd_counts.py "$C" "$C/06b_svd_fit/logs/fit_opt_matrix_6p667_gmm_clean_20260716_125023.log" "$B" --metadata "$B/optics.dat"
   python3 plot_svd_counts.py "$C" "$B"
   python3 plot_svd_clouds.py "$C" rg01
   ```

4. **Freeze strict physical-foil equality plus delta/hole caps (`equal15`).**
   [Exact design](../docs/CORE_BALANCE_DESIGN.md) · [verification](../docs/CORE_BALANCE_VERIFICATION.md) ·
   [campaign recap](05c_core_sample/FOIL_BALANCE_RECAP.md) ·
   [all-foil/all-delta results](05c_core_sample/equal15/HMS_6p667_CORE_BALANCE.md) ·
   [quota code](../core_balance.py) / [allocation code](../reallocate_core.py).
   `./run_core_balance.sh "$C" min10 equal15_review --preview`, then
   `./run_core_balance.sh "$C" min10 equal15`.
   Policy: [config/core_balance.json](config/core_balance.json).

5. **Export exact fit and protected memberships.**
   [Recap/input requirements](06c_core_ntuple/CORE_EXPORT_RECAP.md) · [builder](../build_core_fit.py).
   `./run_build_core_fit.sh "$C" equal15 fit`;
   `./run_build_core_fit.sh "$C" equal15 holdout`.
   Optional admission-only check: `python3 preallocated_svd.py "$C" equal15 "$C/06c_core_ntuple/equal15/acceptance"`.
   This is not the modern SVD solve.

6. **Initial elastic net (`enet`).**
   [Runbook](../docs/ELASTIC_NET.md) · [study decisions](06d_elastic_net/ELASTIC_NET_RECAP.md) ·
   [driver](../fit_elastic.py) / [solver](../elastic_net.py).
   `./run_elastic.sh "$C" equal15 enet --check`, then
   `./run_elastic.sh "$C" equal15 enet`.

7. **Resolve convergence; refit selected terms without penalties (`conv`, `refit`).**
   [Convergence code](../elastic_convergence.py) / [refit code](../elastic_refit.py) ·
   [results](06d_elastic_net/refit/RESULTS.md).
   `./run_elastic.sh "$C" equal15 conv --convergence`;
   `./run_elastic.sh "$C" equal15 refit --refit`.

8. **Pool all three validation folds; choose compact starts (`pooled`).**
   [Results/metric definitions](06d_elastic_net/pooled/RESULTS.md) · [code](../elastic_pool.py).
   `./run_elastic.sh "$C" equal15 pooled --pooled`.
   Review `plots/plateau.png` and `plots/foil_delta.png` in that output.

9. **Beam search with SVD refits: width 4, then width 10.**
   [Algorithm](../docs/BEAM_SEARCH.md) · [decisions/results](06e_beam_search/BEAM_SEARCH_RECAP.md) ·
   [driver](../fit_beam.py) / [search](../beam_search.py).
   `./run_beam.sh "$C" equal15 beam --check`;
   `./run_beam.sh "$C" equal15 beam --beam 4 --steps 20 --threads 8`;
   `./run_beam.sh "$C" equal15 beam10 --beam 10 --steps 40 --threads 8`.

10. **Export surplus; compare five frozen matrices on identical events.**
    [Runbook](../docs/MATRIX_COMPARISON.md) · [identities/findings/limits](07_diagnostics/MATRIX_COMPARISON_RECAP.md) ·
    [numerical report](07_diagnostics/compare/MATRIX_COMPARISON.md) · [code](../compare_matrices.py).
    `./run_build_core_fit.sh "$C" equal15 surplus`;
    `./run_compare.sh "$C" equal15 compare --source beam10 --check`;
    `./run_compare.sh "$C" equal15 compare --source beam10 --threads 8`.

11. **Full-electron-sample GMM/beam sieve afterburner, 12.490°, three foils.**
    [Runbook](../docs/SIEVE_AFTERBURNER.md) · [actual closure warning](07_diagnostics/SIEVE_GEOMETRY_RECAP.md) ·
    [code](../sieve_afterburner.py).
    `./run_sieve.sh "$C" --check`; `./run_sieve.sh "$C"`.
    Already run on ifarm; closure remains unresolved.

12. **Geometry consistency proof; no reconstruction change.**
    [Definitions/equations/references](../docs/HMS_GEOMETRY_CONSISTENCY.md) ·
    [target code](../spectrometer_config.h) / [diagnostic equations](../elastic_diagnostics.py).
    No production command; this was an inspection and algebraic check.

13. **Fix comparison plots without reevaluation.**
    [Plot guide/command](../docs/MATRIX_COMPARISON.md#refresh-existing-plots-quickly) ·
    [plot code](../comparison_plots.py).
    `./run_compare.sh "$C" equal15 compare --replot`.
    Code pushed in `2634d39`; farm refresh not yet confirmed here.

## Is the record sufficient?

**Methods and choices: yes with the linked recaps. A GitHub-only exact rerun: no.**
The guides already specify the algorithms; the recaps add the accepted/rejected
choices, chronology, test evidence and unresolved qualifications. Earlier full-SVD
conditioning is documented separately in
[08_preliminary_conditioning/REPRODUCTION.md](08_preliminary_conditioning/REPRODUCTION.md)
and [PLAIN_LANGUAGE_CLARIFICATIONS.md](08_preliminary_conditioning/PLAIN_LANGUAGE_CLARIFICATIONS.md).

| Required record | Availability at this audit |
|---|---|
| Frozen `05c_core_sample/min10/` | Seven event masks, models, manifest and supporting records tracked; do not re-estimate regions for reallocation. |
| `05c_core_sample/equal15/` | Review reports/maps/tables tracked; original event-level manifest and masks remain on ifarm. |
| `06c_core_ntuple/equal15/{fit,holdout,surplus}/` | Original build manifests, metadata/seed and TFit ROOT files required; not tracked. |
| `06d_elastic_net/`, `06e_beam_search/` | Saved study reports, manifests and numerical products tracked; verify each manifest's inputs before rerunning. |
| `07_diagnostics/compare/residuals.npz` | On ifarm; required for the fast replot. Report/tables/manifest are tracked. |
| Original replay trees; `07_diagnostics/sieve/` | Farm data/products; sieve screenshots and closure numbers were supplied in chat, not a complete tracked output. Exact original replay matrix/flags/version remain to be verified. |

Keep original manifests/builds with their data: regenerated files can have new
hashes even with equivalent event membership. Do not bypass downstream checksum
checks to mix old and new results. No expensive reruns were performed for this
documentation audit; the recaps identify the existing targeted tests and records.
