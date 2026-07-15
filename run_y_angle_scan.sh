#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 CAMPAIGN RUNGROUP FOIL_INDEX [--dry-run]"
  echo "Example: $0 HMS_5p878GeV rg01_theta16p500_foil0 0"
  echo "Example triple foil center: $0 HMS_EXPERIMENT rg_triplefoil 1"
  exit 2
fi

CAMPAIGN="$1"
RUNGROUP="$2"
FOIL_INDEX="$3"
DRY_RUN="${4:-}"

MACRO="assign_yfp_ypfp_angleScanBands.C"

[[ "$FOIL_INDEX" =~ ^[0-9]+$ ]] || {
  echo "ERROR: FOIL_INDEX must be a nonnegative integer" >&2
  exit 1
}

shopt -s nullglob
configs=("${CAMPAIGN}"/config/rungroups_*_inputs.tsv)
shopt -u nullglob

if [[ ${#configs[@]} -ne 1 ]]; then
  echo "ERROR: expected exactly one rungroups_*_inputs.tsv in ${CAMPAIGN}/config" >&2
  exit 1
fi

RGTSV="${configs[0]}"

row=$(awk -F'\t' -v rg="$RUNGROUP" 'NR > 1 && $1 == rg {print; exit}' "$RGTSV")

[[ -n "$row" ]] || {
  echo "ERROR: rungroup not found: $RUNGROUP" >&2
  exit 1
}

IFS=$'\t' read -r rungroup optics_id hms_angle_deg foils nruns runs rootfile <<< "$row"

[[ "$optics_id" =~ ^[0-9]+$ ]] || {
  echo "ERROR: invalid optics_id: $optics_id" >&2
  exit 1
}

[[ -e "$rootfile" ]] || {
  echo "ERROR: ROOT file not found: $rootfile" >&2
  exit 1
}

expr="${MACRO}(${optics_id},-10,10,\"${rungroup}\",9,1.0,0.18,0.08,0.25,1,0.005,0.08,1.0,true,${FOIL_INDEX},-1,-1,-999,true,-999.0,\"${CAMPAIGN}\",\"${rootfile}\")"

echo "Campaign:   $CAMPAIGN"
echo "Rungroup:   $rungroup"
echo "Optics ID:  $optics_id"
echo "Foil index: $FOIL_INDEX"
echo "ROOT file:  $rootfile"

if [[ "$DRY_RUN" == "--dry-run" ]]; then
  printf 'hcana -b -l -q %q\n' "$expr"
else
  hcana -b -l -q "$expr"
fi
