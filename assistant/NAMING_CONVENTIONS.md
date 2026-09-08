# Campaign naming conventions

Current HMS/SHMS implementation: [campaign switch and per-file changes](../docs/HMS_SHMS_CHANGES.md). The arm is selected solely by the `HMS_…` / `SHMS_…` campaign name. Target foil positions are independent metadata for either arm. Centered SHMS sieve is assumed; real-replay validation is pending. Executable repository files take precedence over historical copied source snippets below or in this folder.

Campaign folders use `<spectrometer>_<momentum>GeV`, with `p` replacing the decimal point: for example, `HMS_6p117GeV` and `SHMS_8p5695GeV`. Use the reported central momentum precision; do not put collaborator names or run numbers in the campaign folder name.

The corresponding campaign input table follows `config/rungroups_<momentum>_inputs.tsv`, for example `SHMS_8p5695GeV/config/rungroups_8p5695_inputs.tsv`. Individual run numbers and run groups belong in that table and the existing optics metadata.
