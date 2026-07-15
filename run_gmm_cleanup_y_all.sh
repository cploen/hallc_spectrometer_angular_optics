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
SCRIPT="${PROJECT_DIR}/gmm_cleanup_yscol_candidates_v2.py"
CANDIDATE_DIR="${CAMPAIGN_DIR}/04a_candidate_trees_y"
OUTPUT_DIR="${CAMPAIGN_DIR}/05a_gmm_cleanup_y"

LOGDIR="${OUTPUT_DIR}/logs"
PLOTDIR="${OUTPUT_DIR}/plots"
ROOTDIR="${OUTPUT_DIR}/root"
TSVDIR="${OUTPUT_DIR}/tsv"

MIN_EVENTS="${MIN_EVENTS:-30}"
MAX_COMPONENTS="${MAX_COMPONENTS:-9}"
KEEP_FRAC="${KEEP_FRAC:-0.90}"

[[ -d "$CAMPAIGN_DIR" ]] || {
  echo "ERROR: missing campaign directory: $CAMPAIGN_DIR"
  exit 1
}

[[ -f "$SCRIPT" ]] || {
  echo "ERROR: missing $SCRIPT"
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

mkdir -p "$LOGDIR" "$PLOTDIR" "$ROOTDIR" "$TSVDIR"

n_ok=0
n_fail=0
n_skip=0

while IFS=$'\t' read -r tag optics_id angle foils nruns runs input_root; do
  [[ -z "$tag" || "$tag" == "rungroup" ]] && continue

  candidate_root="${CANDIDATE_DIR}/root/YscolCandidates_${tag}.root"
  summary_tsv="${CANDIDATE_DIR}/tsv/YscolCandidates_${tag}_summary.tsv"

  if [[ ! -s "$candidate_root" || ! -s "$summary_tsv" ]]; then
    echo "FAIL [$tag]: missing candidate ROOT or TSV"
    ((n_fail+=1))
    continue
  fi

  while read -r foil ndel total; do
    [[ -z "$foil" ]] && continue

    if [[ "$total" -le 0 ]]; then
      echo "SKIP [$tag] foil=$foil ndel=$ndel: no candidates"
      ((n_skip+=1))
      continue
    fi

    logfile="${LOGDIR}/gmm_${tag}_foil${foil}_ndel${ndel}.log"

    echo "RUN  [$tag] foil=$foil ndel=$ndel N=$total"

    command=(
      python3 "$SCRIPT"
      --rungroup "$tag"
      --campaign "$CAMPAIGN"
      --foil "$foil"
      --ndel "$ndel"
      --min-events "$MIN_EVENTS"
      --max-components "$MAX_COMPONENTS"
      --keep-frac "$KEEP_FRAC"
    )

    if [[ "$DRY_RUN" == "--dry-run" ]]; then
      printf '%q ' "${command[@]}"
      printf '\n'
      ((n_ok+=1))
      continue
    fi

    "${command[@]}" >"$logfile" 2>&1
    status=$?

    if [[ $status -eq 0 ]] &&
       grep -q "GMM cleanup complete" "$logfile" &&
       ! grep -qE "Traceback|ERROR|RuntimeError" "$logfile"; then
      echo "OK   [$tag] foil=$foil ndel=$ndel"
      ((n_ok+=1))
    else
      echo "FAIL [$tag] foil=$foil ndel=$ndel — inspect $logfile"
      ((n_fail+=1))
    fi
  done < <(
    awk -F'\t' '
      NR > 1 {
        key=$3 FS $4
        sum[key]+=$6
      }
      END {
        for (key in sum) {
          split(key,a,FS)
          print a[1],a[2],sum[key]
        }
      }
    ' "$summary_tsv" | sort -n -k1,1 -k2,2
  )

done < "$RUNGROUPS_TSV"

echo
echo "Finished: ok=$n_ok skipped=$n_skip failed=$n_fail"
[[ $n_fail -eq 0 ]]
