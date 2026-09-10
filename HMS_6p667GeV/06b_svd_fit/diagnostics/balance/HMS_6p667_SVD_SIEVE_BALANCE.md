# Historical SVD sieve counts

Saved fit: July 16, 2026, `6p667_gmm_clean`.

The right-hand maps show the reconstructed events admitted to SVD after its limits. The left-hand maps show all labeled events in the original GMM-cleaned TFit inputs. These are not the new core samples.

Verified 166,013 selected events from 444,697 input events, matching all 540 Y-column counts in the original log. All seven input trees and the original metadata match the SHA-256 checksums in the earlier reproduction record.

Selection reads TFit entries in order, retaining at most 1,000 per rungroup/foil/delta/Y column and 15,000 per rungroup/foil, with a 200,000 campaign maximum. There is no random draw in the SVD selection loop. The original event-ID list was not saved; membership is reconstructed from the unchanged inputs and retained selection code. `selected.npz` saves reconstructed TFit indices and replay-entry IDs by optics ID.

The current metadata file includes later changes, so `optics.dat` preserves the original checksum-matched metadata from input commit 5cb66ee. No fitting or matrix changes were performed.

Known blocked positions are marked, but any historical events assigned to them are still displayed and included in totals: 126 available and 45 selected campaign-wide. This prevents the diagnostic from silently correcting the historical sample.

## Reproduce

From the repository root (PyROOT and NumPy for extraction; matplotlib and NumPy for plotting):

```bash
python3 diagnostics/conditioning/svd_counts.py HMS_6p667GeV HMS_6p667GeV/06b_svd_fit/logs/fit_opt_matrix_6p667_gmm_clean_20260716_125023.log HMS_6p667GeV/06b_svd_fit/diagnostics/balance --metadata HMS_6p667GeV/06b_svd_fit/diagnostics/balance/optics.dat
python3 plot_svd_counts.py HMS_6p667GeV HMS_6p667GeV/06b_svd_fit/diagnostics/balance
```

## Maps

### rg01_theta12p490_foil0

Foil 0: [delta 0](plots/svd_rg01_theta12p490_foil0_foil0_ndel0.png) · [delta 1](plots/svd_rg01_theta12p490_foil0_foil0_ndel1.png) · [delta 2](plots/svd_rg01_theta12p490_foil0_foil0_ndel2.png) · [delta 3](plots/svd_rg01_theta12p490_foil0_foil0_ndel3.png) · [delta 4](plots/svd_rg01_theta12p490_foil0_foil0_ndel4.png)

### rg02_theta15p195_foil0

Foil 0: [delta 0](plots/svd_rg02_theta15p195_foil0_foil0_ndel0.png) · [delta 1](plots/svd_rg02_theta15p195_foil0_foil0_ndel1.png) · [delta 2](plots/svd_rg02_theta15p195_foil0_foil0_ndel2.png) · [delta 3](plots/svd_rg02_theta15p195_foil0_foil0_ndel3.png) · [delta 4](plots/svd_rg02_theta15p195_foil0_foil0_ndel4.png)

### rg03_theta12p490_foilpm8

Foil 0: [delta 0](plots/svd_rg03_theta12p490_foilpm8_foil0_ndel0.png) · [delta 1](plots/svd_rg03_theta12p490_foilpm8_foil0_ndel1.png) · [delta 2](plots/svd_rg03_theta12p490_foilpm8_foil0_ndel2.png) · [delta 3](plots/svd_rg03_theta12p490_foilpm8_foil0_ndel3.png) · [delta 4](plots/svd_rg03_theta12p490_foilpm8_foil0_ndel4.png)

Foil 1: [delta 0](plots/svd_rg03_theta12p490_foilpm8_foil1_ndel0.png) · [delta 1](plots/svd_rg03_theta12p490_foilpm8_foil1_ndel1.png) · [delta 2](plots/svd_rg03_theta12p490_foilpm8_foil1_ndel2.png) · [delta 3](plots/svd_rg03_theta12p490_foilpm8_foil1_ndel3.png) · [delta 4](plots/svd_rg03_theta12p490_foilpm8_foil1_ndel4.png)

### rg04_theta15p195_foilpm3

Foil 0: [delta 0](plots/svd_rg04_theta15p195_foilpm3_foil0_ndel0.png) · [delta 1](plots/svd_rg04_theta15p195_foilpm3_foil0_ndel1.png) · [delta 2](plots/svd_rg04_theta15p195_foilpm3_foil0_ndel2.png) · [delta 3](plots/svd_rg04_theta15p195_foilpm3_foil0_ndel3.png) · [delta 4](plots/svd_rg04_theta15p195_foilpm3_foil0_ndel4.png)

Foil 1: [delta 0](plots/svd_rg04_theta15p195_foilpm3_foil1_ndel0.png) · [delta 1](plots/svd_rg04_theta15p195_foilpm3_foil1_ndel1.png) · [delta 2](plots/svd_rg04_theta15p195_foilpm3_foil1_ndel2.png) · [delta 3](plots/svd_rg04_theta15p195_foilpm3_foil1_ndel3.png) · [delta 4](plots/svd_rg04_theta15p195_foilpm3_foil1_ndel4.png)

### rg05_theta12p495_foilpm3

Foil 0: [delta 0](plots/svd_rg05_theta12p495_foilpm3_foil0_ndel0.png) · [delta 1](plots/svd_rg05_theta12p495_foilpm3_foil0_ndel1.png) · [delta 2](plots/svd_rg05_theta12p495_foilpm3_foil0_ndel2.png) · [delta 3](plots/svd_rg05_theta12p495_foilpm3_foil0_ndel3.png) · [delta 4](plots/svd_rg05_theta12p495_foilpm3_foil0_ndel4.png)

Foil 1: [delta 0](plots/svd_rg05_theta12p495_foilpm3_foil1_ndel0.png) · [delta 1](plots/svd_rg05_theta12p495_foilpm3_foil1_ndel1.png) · [delta 2](plots/svd_rg05_theta12p495_foilpm3_foil1_ndel2.png) · [delta 3](plots/svd_rg05_theta12p495_foilpm3_foil1_ndel3.png) · [delta 4](plots/svd_rg05_theta12p495_foilpm3_foil1_ndel4.png)

### rg06_theta15p195_foilpm8

Foil 0: [delta 0](plots/svd_rg06_theta15p195_foilpm8_foil0_ndel0.png) · [delta 1](plots/svd_rg06_theta15p195_foilpm8_foil0_ndel1.png) · [delta 2](plots/svd_rg06_theta15p195_foilpm8_foil0_ndel2.png) · [delta 3](plots/svd_rg06_theta15p195_foilpm8_foil0_ndel3.png) · [delta 4](plots/svd_rg06_theta15p195_foilpm8_foil0_ndel4.png)

Foil 1: [delta 0](plots/svd_rg06_theta15p195_foilpm8_foil1_ndel0.png) · [delta 1](plots/svd_rg06_theta15p195_foilpm8_foil1_ndel1.png) · [delta 2](plots/svd_rg06_theta15p195_foilpm8_foil1_ndel2.png) · [delta 3](plots/svd_rg06_theta15p195_foilpm8_foil1_ndel3.png) · [delta 4](plots/svd_rg06_theta15p195_foilpm8_foil1_ndel4.png)

### rg07_theta12p495_foilpm8

Foil 0: [delta 0](plots/svd_rg07_theta12p495_foilpm8_foil0_ndel0.png) · [delta 1](plots/svd_rg07_theta12p495_foilpm8_foil0_ndel1.png) · [delta 2](plots/svd_rg07_theta12p495_foilpm8_foil0_ndel2.png) · [delta 3](plots/svd_rg07_theta12p495_foilpm8_foil0_ndel3.png) · [delta 4](plots/svd_rg07_theta12p495_foilpm8_foil0_ndel4.png)

Foil 1: [delta 0](plots/svd_rg07_theta12p495_foilpm8_foil1_ndel0.png) · [delta 1](plots/svd_rg07_theta12p495_foilpm8_foil1_ndel1.png) · [delta 2](plots/svd_rg07_theta12p495_foilpm8_foil1_ndel2.png) · [delta 3](plots/svd_rg07_theta12p495_foilpm8_foil1_ndel3.png) · [delta 4](plots/svd_rg07_theta12p495_foilpm8_foil1_ndel4.png)


## Event-cloud comparison

[rg01 three-stage event-cloud plots](HMS_6p667_RG01_EVENT_CLOUDS.md)
