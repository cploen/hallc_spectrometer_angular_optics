#!/usr/bin/env bash
set -uo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 CAMPAIGN [--dry-run]"
  exit 2
fi

PROJECT_DIR=$(pwd -P)
CAMPAIGN="$1"
DRY_RUN="${2:-}"

if [[ -n "$DRY_RUN" && "$DRY_RUN" != "--dry-run" ]]; then
  echo "ERROR: unknown argument: $DRY_RUN"
  exit 2
fi

CAMPAIGN_DIR="${PROJECT_DIR}/${CAMPAIGN}"
HELPER="./run_relabel_xfp_from_plot_splits_v2.sh"

[[ -d "$CAMPAIGN_DIR" ]] || {
  echo "ERROR: missing campaign directory: $CAMPAIGN_DIR"
  exit 1
}

[[ -x "$HELPER" ]] || {
  echo "ERROR: missing executable helper: $HELPER"
  exit 1
}

shopt -s nullglob
RUNGROUP_FILES=("${CAMPAIGN_DIR}"/config/rungroups_*_inputs.tsv)
shopt -u nullglob

if [[ ${#RUNGROUP_FILES[@]} -ne 1 ]]; then
  echo "ERROR: expected exactly one rungroups_*_inputs.tsv in ${CAMPAIGN_DIR}/config"
  exit 1
fi

RUNGROUPS_TSV="${RUNGROUP_FILES[0]}"

RELABEL_DIR="${CAMPAIGN_DIR}/03b_relabel_x"
mkdir -p \
  "${RELABEL_DIR}/logs" \
  "${RELABEL_DIR}/cuts" \
  "${RELABEL_DIR}/plots" \
  "${RELABEL_DIR}/root"

n_ok=0
n_fail=0

while IFS=$'\t' read -r rungroup optics_id angle foils nruns runs input_root; do
  [[ -z "$rungroup" || "$rungroup" == "rungroup" ]] && continue

  if [[ -z "$optics_id" || ! "$optics_id" =~ ^[0-9]+$ ]]; then
    echo "FAIL [$rungroup]: invalid optics_id '$optics_id'"
    ((n_fail+=1))
    continue
  fi

  if [[ ! -f "$input_root" ]]; then
    echo "FAIL [$rungroup]: input ROOT missing: $input_root"
    ((n_fail+=1))
    continue
  fi

  logfile="${RELABEL_DIR}/logs/xrelabel_wrapper_${rungroup}.log"

  echo
  echo "============================================================"
  echo "X RELABEL:   $rungroup"
  echo "Optics ID:   $optics_id"
  echo "Input ROOT:  $input_root"
  echo "============================================================"

  if [[ "$DRY_RUN" == "--dry-run" ]]; then
    METADATA_RUN="$optics_id" \
    RUNGROUP="$rungroup" \
    CAMPAIGN="$CAMPAIGN" \
    INPUT_ROOT="$input_root" \
    DRYRUN=1 \
    "$HELPER"

    status=$?
  else
    METADATA_RUN="$optics_id" \
    RUNGROUP="$rungroup" \
    CAMPAIGN="$CAMPAIGN" \
    INPUT_ROOT="$input_root" \
    "$HELPER" 2>&1 | tee "$logfile"

    status=${PIPESTATUS[0]}
  fi

  if [[ $status -eq 0 ]]; then
    ((n_ok+=1))
  else
    echo "FAIL [$rungroup]"
    ((n_fail+=1))
  fi
done < "$RUNGROUPS_TSV"

echo
echo "Finished: ok=$n_ok failed=$n_fail"
[[ $n_fail -eq 0 ]]
