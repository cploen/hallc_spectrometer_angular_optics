# HMS 6.667 GeV — rg01 event clouds

Original joint-labeled events → GMM survivors → events admitted to the saved July 16, 2026 SVD fit. All five figures use the same axes, binning and log color scale. No density boundaries or new core-selection cuts are applied.

All 207,731 GMM survivors match the original candidates by event ID, sieve coordinates, foil, delta slice, and X/Y labels. The original joint set contains 232,528 events; 48 candidate entry IDs were excluded by the existing unique-joint-label checks. SVD membership is reconstructed from checksum-matched inputs and the retained selection rules; the campaign total and all logged Y-column counts match. No original SVD event-ID list was saved.

| Delta slice | Original labeled | After GMM | Used by SVD | Plot |
|---|---:|---:|---:|---|
| 0 | 22,100 | 19,460 | 1,468 | [Open slice 0](clouds/clouds_rg01_theta12p490_foil0_foil0_ndel0.png) |
| 1 | 38,901 | 34,543 | 2,572 | [Open slice 1](clouds/clouds_rg01_theta12p490_foil0_foil0_ndel1.png) |
| 2 | 67,428 | 60,399 | 4,095 | [Open slice 2](clouds/clouds_rg01_theta12p490_foil0_foil0_ndel2.png) |
| 3 | 60,696 | 54,563 | 3,852 | [Open slice 3](clouds/clouds_rg01_theta12p490_foil0_foil0_ndel3.png) |
| 4 | 43,403 | 38,766 | 3,013 | [Open slice 4](clouds/clouds_rg01_theta12p490_foil0_foil0_ndel4.png) |

## Reproduce

From the repository root, using the core-sample Python dependencies:

```bash
python3 plot_svd_clouds.py HMS_6p667GeV rg01
```

This reads the two rg01 candidate ROOT files, the original `Optics_666701_-1_fit_tree_gmm.root`, and the verified `selected.npz` described in [HMS_6p667_SVD_SIEVE_BALANCE.md](HMS_6p667_SVD_SIEVE_BALANCE.md). Rerunning replaces only the cloud figures and their provenance JSON. Event counts include any points outside the displayed frame; those counts appear above each panel.
