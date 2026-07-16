#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 CAMPAIGN [--dry-run]

CAMPAIGN may be a campaign directory path or a campaign directory name
under the repository root.

Environment overrides:
  FILEID=-1
  MAX_PER_Y_BIN=1000
  MAX_PER_FOIL=15000
  MAX_PER_RUN_FOIL=30000
  NFIT_MAX=200000
  OPTICS_METADATA=/path/to/list_of_optics_run.dat
USAGE
}

[[ $# -ge 1 && $# -le 2 ]] || { usage; exit 2; }
CAMPAIGN_ARG=$1
MODE=${2:-}
[[ -z "$MODE" || "$MODE" == "--dry-run" ]] || { usage; exit 2; }
DRY_RUN=0
[[ "$MODE" == "--dry-run" ]] && DRY_RUN=1

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
REPO_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd -P)

if [[ -d "$CAMPAIGN_ARG" ]]; then
  CAMPAIGN_DIR=$(cd -- "$CAMPAIGN_ARG" && pwd -P)
elif [[ -d "$REPO_DIR/$CAMPAIGN_ARG" ]]; then
  CAMPAIGN_DIR=$(cd -- "$REPO_DIR/$CAMPAIGN_ARG" && pwd -P)
else
  echo "ERROR: campaign directory not found: $CAMPAIGN_ARG" >&2
  exit 1
fi

FILEID=${FILEID:--1}
MAX_PER_Y_BIN=${MAX_PER_Y_BIN:-1000}
MAX_PER_FOIL=${MAX_PER_FOIL:-15000}
MAX_PER_RUN_FOIL=${MAX_PER_RUN_FOIL:-30000}
NFIT_MAX=${NFIT_MAX:-200000}
OPTICS_METADATA=${OPTICS_METADATA:-$REPO_DIR/DATfiles/list_of_optics_run.dat}

FIT_TREE_DIR="$CAMPAIGN_DIR/06a_fit_ntuple/root"
OUTPUT_DIR="$CAMPAIGN_DIR/07_diagnostics/fit"
BALANCE_DIR="$OUTPUT_DIR/sample_balance"
RESIDUAL_DIR="$OUTPUT_DIR/residuals"
LOG_DIR="$OUTPUT_DIR/logs"

[[ -d "$FIT_TREE_DIR" ]] || { echo "ERROR: missing fit-tree directory: $FIT_TREE_DIR" >&2; exit 1; }
[[ -f "$OPTICS_METADATA" ]] || { echo "ERROR: missing optics metadata: $OPTICS_METADATA" >&2; exit 1; }

shopt -s nullglob
RUNGROUP_FILES=("$CAMPAIGN_DIR"/config/rungroups_*_inputs.tsv)
shopt -u nullglob
[[ ${#RUNGROUP_FILES[@]} -eq 1 ]] || {
  echo "ERROR: expected exactly one rungroups_*_inputs.tsv in $CAMPAIGN_DIR/config" >&2
  exit 1
}
RUNGROUPS_TSV=${RUNGROUP_FILES[0]}

mkdir -p "$BALANCE_DIR" "$RESIDUAL_DIR" "$LOG_DIR"
STAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$LOG_DIR/fit_diagnostics_${STAMP}.log"

run_cmd() {
  printf '%q ' "$@"
  printf '\n'
  if [[ $DRY_RUN -eq 0 ]]; then
    "$@"
  fi
}

{
  echo "Campaign:        $CAMPAIGN_DIR"
  echo "Rungroups TSV:   $RUNGROUPS_TSV"
  echo "Fit-tree dir:    $FIT_TREE_DIR"
  echo "Optics metadata: $OPTICS_METADATA"
  echo "Output dir:      $OUTPUT_DIR"
  echo

  balance_expr="$REPO_DIR/analyze_fit_sample_balance.C(\"$CAMPAIGN_DIR\",$FILEID,$MAX_PER_Y_BIN,$MAX_PER_FOIL,$MAX_PER_RUN_FOIL,$NFIT_MAX,\"$BALANCE_DIR\")"
  run_cmd root -l -b -q "$balance_expr"

  while IFS=$'\t' read -r rungroup optics_id angle foils nruns runs input_root; do
    [[ -z "$rungroup" || "$rungroup" == "rungroup" ]] && continue
    [[ "$optics_id" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid optics ID for $rungroup: $optics_id"; exit 1; }

    representative_run=${runs%%,*}
    representative_run=${representative_run//[[:space:]]/}
    [[ "$representative_run" =~ ^[0-9]+$ ]] || {
      echo "ERROR: cannot determine representative run for $rungroup from '$runs'"
      exit 1
    }

    fit_tree="$FIT_TREE_DIR/Optics_${optics_id}_${FILEID}_fit_tree_gmm.root"
    [[ -s "$fit_tree" ]] || { echo "ERROR: missing fit tree: $fit_tree"; exit 1; }

    tag=${rungroup//[^A-Za-z0-9_.-]/_}
    rungroup_out="$RESIDUAL_DIR/$tag"
    mkdir -p "$rungroup_out"

    for variable in xptar yptar ytar; do
      macro="$REPO_DIR/diagnostics/fit/plot_${variable}_residuals.C"
      [[ -f "$macro" ]] || { echo "ERROR: missing macro: $macro"; exit 1; }
      expr="$macro($representative_run,\"$fit_tree\",\"$OPTICS_METADATA\",\"$rungroup_out\",\"$tag\")"
      run_cmd root -l -b -q "$expr"
    done
  done < "$RUNGROUPS_TSV"

  echo
  echo "Fit diagnostics complete: $OUTPUT_DIR"
} 2>&1 | tee "$LOG_FILE"
