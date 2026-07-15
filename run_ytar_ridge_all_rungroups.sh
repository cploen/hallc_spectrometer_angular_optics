#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 CAMPAIGN [SKIP_REGEX]"
  echo "Example: $0 HMS_5p878GeV"
  exit 2
fi

CAMPAIGN="$1"
SKIP_RE="${2:-^$}"

if [[ ! -d "$CAMPAIGN" ]]; then
  echo "ERROR: campaign directory not found: $CAMPAIGN" >&2
  exit 1
fi

shopt -s nullglob
configs=("${CAMPAIGN}"/config/rungroups_*_inputs.tsv)
shopt -u nullglob

if [[ ${#configs[@]} -ne 1 ]]; then
  echo "ERROR: expected exactly one rungroups_*_inputs.tsv in ${CAMPAIGN}/config" >&2
  printf 'Found: %s\n' "${configs[@]:-none}" >&2
  exit 1
fi

RGTSV="${configs[0]}"
LOGDIR="${CAMPAIGN}/01_ytar_cuts/notes/logs"
mkdir -p "$LOGDIR"

echo "Campaign: $CAMPAIGN"
echo "Config:   $RGTSV"
echo "Logs:     $LOGDIR"

while IFS=$'\t' read -r rungroup optics_id hms_angle_deg foils nruns runs rootfile; do
  [[ -z "$rungroup" ]] && continue

  if [[ "$rungroup" =~ $SKIP_RE ]]; then
    echo "Skipping: $rungroup"
    continue
  fi

  if [[ ! "$optics_id" =~ ^[0-9]+$ ]]; then
    echo "ERROR: invalid optics_id '$optics_id' for $rungroup" >&2
    exit 1
  fi

  if [[ ! -e "$rootfile" ]]; then
    echo "ERROR: ROOT file not found for $rungroup:" >&2
    echo "       $rootfile" >&2
    exit 1
  fi

  echo "Running ytar ridge cut: $rungroup (optics ID $optics_id)"
  log="${LOGDIR}/ytar_ridge_${rungroup}_$(date +%Y%m%d_%H%M%S).log"

  hcana -b -l -q \
    "ytar_ridge_cut.C(${optics_id},\"${rungroup}\",\"-1\",\"-1\",\"${rootfile}\",\"${CAMPAIGN}\")" \
    2>&1 | tee "$log"

done < <(tail -n +2 "$RGTSV")
