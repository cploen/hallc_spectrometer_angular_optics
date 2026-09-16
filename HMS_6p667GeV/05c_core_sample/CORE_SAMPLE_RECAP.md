# HMS 6.667 — core-selection decisions and review

Conversation recap, 15 September 2026. Complements the method and commands in
[CORE_SAMPLE.md](../../CORE_SAMPLE.md); it does not redefine the frozen samples.

## Choices made

- Start from **joint X/Y candidate labels before GMM**, retaining exact
  `(rungroup, entry)` identities. Assume the inherited labels are meaningful;
  the selector does not relabel events or prove label correctness.
- Reserve 20% within each labeled hole/foil/delta/rungroup stratum **before**
  learning density regions. Use deterministic hashing, seed 667. The same learned
  regions classify the reserve afterward. The reserve therefore contains core,
  noncore, and unsupported events.
- Choose label-seeded smoothed density and watershed regions rather than a
  global density cut or HDBSCAN. Each hole uses its own peak/background contrast,
  allowing sparse edge holes to survive alongside densely populated holes.
- The initial `core` tag used `min_events=20`; the accepted sparse-edge trial
  `min10` uses 10 **development events needed to attempt a region**, not ten fit
  events per hole. Other support/contrast/stability tests remain in force.
- Leave unusual contour shapes alone: recovering useful edge populations was
  the priority. A contour only selects events already carrying its label. It
  cannot itself change that label, but it also cannot certify inherited labels.
- Keep unused eligible core **events** as surplus. One hole has one accepted
  core region, potentially containing many events. Surplus is development data
  because it helped establish the region; it is not the protected reserve.

## Reviews and adjustments

The rg01 slices looked successful. Sparse left/edge holes in rg02 and several
outer-foil settings motivated the 20→10 change (`2c046d7`). The resulting `min10`
selection was accepted for further fitting, without claiming optimal thresholds.
Its complete frozen records were later exposed (`ce10092`).

Count histograms did not adequately convey spatial acceptance. The review moved
to numbered circles at nominal hole positions, log colors, available/training/
proposed counts, and a 1.5×P25 cap for every rg01 delta slice. These were diagnostic
proposals, not changes to the original `min10` allocation. “Available” excludes
protected events, so a proposed count equal to availability does not consume
the protected reserve.

Blocked HMS labels `[2,3]` and `[5,5]` were confirmed by the user and saved in
`config/sieve_mask.json`; see [the schematic](../../docs/HMS-Sieve.pdf).
Absent labels, observed labels without eligible cores, and blocked positions
are different categories. An absent label is not proof of an illuminated hole
that the selector lost. The legacy selector did not apply the geometry mask;
the later strict allocator does.

The initial 80%/400-event/200,000 allocation was superseded for fitting by
[equal15 foil balancing](FOIL_BALANCE_RECAP.md). Density classifications and the
protected split stayed frozen. Density centroids select events; **nominal
physical hole coordinates remain the fit targets**.

## Records and checks

- Frozen definitions: `core/config.json`, `min10/config.json`, their manifests,
  `tsv/regions.tsv`, `tsv/counts.tsv`, `models/`, and `root/`.
- Code: `core_sample.py`, `plot_core_counts.py`, `plot_core_map.py` at repository root.
- Tests: `tests/test_core_sample.py` covers holdout independence, sparse/missing
  populations, repeatability, duplicate labels, 64-bit IDs and budget separation.
- No fresh final reserve can be created by renaming a tag. These held-out plots
  have informed development; later residual reviews are not an untouched final
  test. No additional data can be collected for this campaign.
