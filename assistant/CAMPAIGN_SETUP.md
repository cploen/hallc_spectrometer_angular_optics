# Campaign Setup

A campaign is a spectrometer optics-calibration dataset organized around one momentum setting or related group of runs.

The repository root is referred to as:

<repo>

The campaign directory is referred to as:

<campaign>

Example:

<repo>/HMS_6p117GeV

## Important distinction

Campaign setup currently involves three different kinds of information.

### 1. Physics and run metadata

The primary shared metadata file is:

<repo>/DATfiles/list_of_optics_run.dat

It contains run-level optics information such as:

- run number;
- optics identifier;
- central angle;
- number of foils;
- sieve flag;
- delta boundaries;
- foil positions;
- y-misalignment information.

Many ROOT macros read this file directly.

Before processing a new campaign, verify that all required runs and their foil and delta metadata are present and correct.

### 2. Campaign and run-group metadata

The current 6.117 GeV/c workflow also uses campaign-specific TSV files, including:

- rungroups_6p117.tsv
- rungroups_6p117_inputs.tsv
- run_manifest_6p117.tsv
- run_foil_manifest_6p117.tsv
- runlist_6p117_use.txt

These files describe:

- which runs are usable;
- which runs belong to each run group;
- spectrometer angle;
- foil set;
- expected grouped ROOT filename;
- full input ROOT path in some cases.

These campaign-specific files currently overlap in purpose.

Do not assume they are synchronized automatically.

### 3. User-specific ROOT-file location

The replayed ROOT files may live in any user-accessible location, including:

- /volatile
- /work
- another shared filesystem
- another local or project directory

The ROOT-file path is not universal and must not be hard-coded in generalized documentation or assistant answers.

The current 6.117 GeV/c workflow stores one example path in:

ML_angular_6p117_20260630/config/rootdir_6p117.txt

This path belongs to one user and one replay campaign.

## Two valid input modes

### Mode A: Direct ROOT-file mode

Use this when each calibration input is already represented by one ROOT file.

The user must provide:

- campaign name;
- spectrometer;
- momentum setting;
- run or logical dataset identifier;
- ROOT-file path;
- foil configuration;
- spectrometer angle;
- delta boundaries;
- corresponding metadata in list_of_optics_run.dat.

No hadd step is required.

### Mode B: Run-group and hadd mode

Use this when multiple replayed ROOT files should be combined into one logical calibration dataset.

The user must provide:

- source ROOT directory;
- run-group name;
- run list;
- spectrometer angle;
- foil configuration;
- output grouped ROOT filename;
- campaign name;
- run metadata.

The grouped file is typically created with hadd.

The current 6.117 GeV/c naming pattern is:

nps_hms_optics_6p117_<rungroup>.root

This pattern is current practice, not yet a universal spectrometer-independent contract.

## Run-group naming

The current 6.117 GeV/c workflow uses names such as:

rg01_theta12p385_foil0
rg04_theta12p395_foilpm3
rg06_theta12p385_foilpm8

The current components are:

- sequential run-group identifier;
- central angle encoded with p for decimal point;
- foil configuration.

Examples:

- theta12p385 means 12.385 degrees;
- foil0 means the central foil;
- foilpm3 means the ±3 cm foil pair;
- foilpm8 means the ±8 cm foil pair.

Future generalized runners should treat the run-group string as a campaign-defined tag rather than reconstructing it blindly.

## Recommended campaign directory structure

A campaign directory should contain:

<campaign>/
  00_inputs/
  01_ytar_cuts/
  02a_angle_scan_y/
  02b_angle_scan_x/
  03a_relabel_y/
  03b_relabel_x/
  04a_candidate_trees_y/
  04b_candidate_trees_x/
  05a_gmm_cleanup_y/
  05b_gmm_cleanup_x/
  06a_fit_ntuple/
  06b_svd_fit/
  issues/
  README.md

The current HMS_6p117GeV campaign does not yet have a 00_inputs directory.

For future campaigns, 00_inputs should contain campaign-local configuration or references to it, such as:

- rootdir.txt
- run_manifest.tsv
- rungroups.tsv
- rungroups_inputs.tsv
- runlist_use.txt
- README.md

Do not duplicate physics metadata unnecessarily if it already exists authoritatively in:

<repo>/DATfiles/list_of_optics_run.dat

## Minimum campaign README

Each campaign README should record:

- campaign name;
- spectrometer;
- central momentum;
- replay source;
- ROOT-file directory;
- whether direct-file or hadd mode is used;
- run groups;
- foil configurations;
- delta boundaries;
- source matrix;
- replay date or tag;
- known bad runs;
- known metadata uncertainties;
- current workflow status;
- nondefault parameters;
- important manual decisions.

## Hadd output directory

For hadd-based campaigns, a user may create a grouped-input directory outside the repository, for example:

<user_rootdir>/hadd_rungroups_<setting>/

That directory should contain:

- grouped ROOT files;
- README describing source runs and commands;
- optionally a copy or snapshot of the run-group manifest.

The repository should store the configuration and provenance, not necessarily the large ROOT files.

## Current technical debt

The current code contains several hard-coded 6.117 GeV/c and user-specific paths.

Examples include:

- 6.117-specific shell runners;
- Python runners with fixed ROOTDIR values;
- generated xrelabel_chain files containing absolute paths;
- macros with user-specific default rootDir values.

The assistant must warn users when a command depends on one of these hard-coded paths.

Generalized runners should accept the following explicitly:

- campaign;
- spectrometer;
- ROOT input directory;
- run-group manifest;
- optics metadata file;
- output campaign directory.

## New-campaign checklist

Before generating commands, determine:

1. spectrometer;
2. campaign name;
3. momentum setting;
4. replay ROOT location;
5. direct-file or hadd mode;
6. run list;
7. run-group definitions, if used;
8. foil configuration;
9. central angle;
10. delta boundaries;
11. metadata presence in list_of_optics_run.dat;
12. output campaign directory;
13. whether any current script still contains a hard-coded path.

Do not generate a final runnable command until all required values for that stage are known.
