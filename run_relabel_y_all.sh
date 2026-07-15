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
MACRO="relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch.C"

MIN_EVENTS="${MIN_EVENTS:-30}"
OVERWRITE="${OVERWRITE:-0}"

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

STEP_DIR="${CAMPAIGN_DIR}/03a_relabel_y"
LOGDIR="${STEP_DIR}/logs"
CUTDIR="${STEP_DIR}/cuts"
PLOTDIR="${STEP_DIR}/plots"
ROOTDIR="${STEP_DIR}/root"

mkdir -p "$LOGDIR" "$CUTDIR" "$PLOTDIR" "$ROOTDIR"

[[ -f "$RUNGROUPS_TSV" ]] || {
  echo "ERROR: missing $RUNGROUPS_TSV"
  exit 1
}

[[ -f "$DATFILE" ]] || {
  echo "ERROR: missing $DATFILE"
  exit 1
}

[[ -f "$MACRO" ]] || {
  echo "ERROR: missing $MACRO"
  exit 1
}

n_ok=0
n_skip=0
n_fail=0

while IFS=$'\t' read -r tag optics_id angle foils nruns runs input_root; do
  [[ -z "$tag" || "$tag" == "rungroup" ]] && continue

  if [[ -z "$optics_id" || ! "$optics_id" =~ ^[0-9]+$ ]]; then
    echo "FAIL [$tag]: invalid optics_id '$optics_id'"
    ((n_fail+=1))
    continue
  fi

  num_foils=$(awk -F',' -v wanted="$optics_id" '
    {
      id=$1
      nf=$4
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", id)
      gsub(/[[:space:]]/, "", nf)

      if (id == wanted) {
        print nf
        exit
      }
    }
  ' "$DATFILE")

  if [[ -z "$num_foils" || ! "$num_foils" =~ ^[0-9]+$ ]]; then
    echo "FAIL [$tag]: no valid foil count for optics ID $optics_id"
    ((n_fail+=1))
    continue
  fi

  if [[ ! -f "$input_root" ]]; then
    echo "FAIL [$tag]: input ROOT file missing: $input_root"
    ((n_fail+=1))
    continue
  fi

  final_cut="${CUTDIR}/YpFpYFp_${tag}_auto_band_cut_yscolRelabeled.root"

  if [[ "$OVERWRITE" == "1" ]]; then
    rm -f "$final_cut"
    rm -f "${PLOTDIR}/yfp_ypfp_yscolRelabel_${tag}_"*.pdf
    rm -f "${ROOTDIR}/yfp_ypfp_yscolRelabel_${tag}_"*.root
  elif [[ -f "$final_cut" ]]; then
    echo "SKIP [$tag]: final relabeled cut exists"
    ((n_skip+=1))
    continue
  fi

  rungroup_failed=0

  for ((foil=0; foil<num_foils; foil++)); do
    logfile="${LOGDIR}/relabel_y_${tag}_foil${foil}.log"

    expression="${MACRO}(${optics_id},\"${tag}\",\"${CAMPAIGN}\",\"${input_root}\",${foil},0,\"\",true,-1,6.0,0.65,${MIN_EVENTS},true,true)"

    echo "RUN  [$tag] foil=$foil optics_id=$optics_id"

    if [[ "$DRY_RUN" == "--dry-run" ]]; then
      printf 'hcana -b -l -q %q\n' "$expression"
      continue
    fi

    hcana -b -l -q "$expression" >"$logfile" 2>&1
    status=$?

    if [[ $status -ne 0 ]] || grep -qE '(^|[[:space:]])ERROR:' "$logfile"; then
      echo "FAIL [$tag] foil=$foil — inspect $logfile"
      rungroup_failed=1
      break
    fi
  done

  if [[ "$DRY_RUN" == "--dry-run" ]]; then
    ((n_ok+=1))
  elif [[ $rungroup_failed -eq 0 && -s "$final_cut" ]]; then
    echo "OK   [$tag]"
    ((n_ok+=1))
  else
    echo "FAIL [$tag]: final cut missing or empty"
    ((n_fail+=1))
  fi
done < "$RUNGROUPS_TSV"

echo
echo "Finished: ok=$n_ok skipped=$n_skip failed=$n_fail"
[[ $n_fail -eq 0 ]]
