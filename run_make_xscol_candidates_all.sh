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
DATFILE="${PROJECT_DIR}/DATfiles/list_of_optics_run.dat"
MACRO="${PROJECT_DIR}/make_xscol_candidate_tree.C"

[[ -d "$CAMPAIGN_DIR" ]] || {
  echo "ERROR: missing campaign directory: $CAMPAIGN_DIR"
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

STEP_DIR="${CAMPAIGN_DIR}/04b_candidate_trees_x"
LOGDIR="${STEP_DIR}/logs"
ROOTDIR="${STEP_DIR}/root"
TSVDIR="${STEP_DIR}/tsv"

mkdir -p "$LOGDIR" "$ROOTDIR" "$TSVDIR"

[[ -f "$DATFILE" ]] || {
  echo "ERROR: missing $DATFILE"
  exit 1
}

[[ -f "$MACRO" ]] || {
  echo "ERROR: missing $MACRO"
  exit 1
}

n_ok=0
n_fail=0

while IFS=$'\t' read -r tag optics_id angle foils nruns runs input_root; do
  [[ -z "$tag" || "$tag" == "rungroup" ]] && continue

  if [[ -z "$optics_id" || ! "$optics_id" =~ ^[0-9]+$ ]]; then
    echo "FAIL [$tag]: invalid optics_id '$optics_id'"
    ((n_fail+=1))
    continue
  fi

  metadata_exists=$(awk -F',' -v wanted="$optics_id" '
    {
      id=$1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", id)

      if (id == wanted) {
        print 1
        exit
      }
    }
  ' "$DATFILE")

  if [[ "$metadata_exists" != "1" ]]; then
    echo "FAIL [$tag]: optics ID $optics_id not found in $DATFILE"
    ((n_fail+=1))
    continue
  fi

  if [[ ! -f "$input_root" ]]; then
    echo "FAIL [$tag]: missing input ROOT: $input_root"
    ((n_fail+=1))
    continue
  fi

  output_root="${ROOTDIR}/XscolCandidates_${tag}.root"
  output_tsv="${TSVDIR}/XscolCandidates_${tag}_summary.tsv"
  logfile="${LOGDIR}/make_xscol_candidates_${tag}.log"

  expression="make_xscol_candidate_tree.C(${optics_id},\"${tag}\",\"${CAMPAIGN}\",\"${input_root}\")"

  echo "RUN  [$tag] optics_id=$optics_id"

  if [[ "$DRY_RUN" == "--dry-run" ]]; then
    printf 'hcana -b -l -q %q\n' "$expression"
    ((n_ok+=1))
    continue
  fi

  rm -f "$output_root" "$output_tsv"

  hcana -b -l -q "$expression" >"$logfile" 2>&1
  status=$?

  if [[ $status -eq 0 && -s "$output_root" && -s "$output_tsv" ]] &&
     ! grep -qE '(^|[[:space:]])ERROR:' "$logfile"; then
    echo "OK   [$tag]"
    ((n_ok+=1))
  else
    echo "FAIL [$tag] — inspect $logfile"
    ((n_fail+=1))
  fi
done < "$RUNGROUPS_TSV"

echo
echo "Finished: ok=$n_ok failed=$n_fail"
[[ $n_fail -eq 0 ]]
