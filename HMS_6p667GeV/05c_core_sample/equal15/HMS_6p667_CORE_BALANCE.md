# HMS 6.667 core balance

**Frozen core reallocation; no SVD fit run.**

Fixed hole/delta multipliers 1.5/1.5; qa_min=10 is reporting only.

Exact campaign budget: **68,535 = 5 × 13,707**.

| zfoil | available | capacity | quota | actual_selected | hole_loss | delta_loss | budget_loss | surplus | limiting |
|---|---|---|---|---|---|---|---|---|---|
| -8 | 54299 | 33294 | 13707 | 13707 | 21005 | 0 | 19587 | 40592 | False |
| -3 | 23339 | 14475 | 13707 | 13707 | 8864 | 0 | 768 | 9632 | False |
| 0 | 160951 | 92188 | 13707 | 13707 | 68763 | 0 | 78481 | 147244 | False |
| 3 | 25088 | 13707 | 13707 | 13707 | 11381 | 0 | 0 | 11381 | True |
| 8 | 46836 | 24862 | 13707 | 13707 | 21974 | 0 | 11155 | 33129 | False |

- [Delta 0: [-10, -8)%, width 2 pp](pages/delta0.md)
- [Delta 1: [-8, -5)%, width 3 pp](pages/delta1.md)
- [Delta 2: [-5, 0)%, width 5 pp](pages/delta2.md)
- [Delta 3: [0, 5)%, width 5 pp](pages/delta3.md)
- [Delta 4: [5, 10)%, width 5 pp](pages/delta4.md)

[Exact missing inputs](INPUT_AVAILABILITY.md) · [All leaves](tsv/leaves.tsv) · [Setting totals](tsv/settings.tsv) · [QA: populated holes below 10](tsv/qa_low.tsv) · [Delta capacities](tsv/deltas.tsv)

P25 references use original positive eligible counts. Geometry exclusions precede references. No automatic relaxation, reserve, weights, or legacy 80%/400 caps. Surplus includes all unselected eligible core events.

Equal foil totals do not enforce matching delta cells or equal setting contributions. Count equality does not establish equal leverage, conditioning, or matrix stability.

Historical comparison uses different populations: historical GMM → SVD versus frozen accepted development cores. Historical membership was reconstructed and verified against the saved log; original solver ID lists were not saved.

[Actual historical foil/delta totals alongside new quotas](tsv/historical_comparison.tsv)

- [Preserved historical context: HMS_6p667_SVD_SIEVE_BALANCE.md](../../06b_svd_fit/diagnostics/balance/HMS_6p667_SVD_SIEVE_BALANCE.md)
- [Preserved historical context: HMS_6p667_RG01_EVENT_CLOUDS.md](../../06b_svd_fit/diagnostics/balance/HMS_6p667_RG01_EVENT_CLOUDS.md)
