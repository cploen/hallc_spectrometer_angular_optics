#!/usr/bin/env bash
# Generate protected beam plots and legacy afterburner PDFs, campaign by campaign.
set -euo pipefail
PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
cd "$PROJECT_DIR"
if ! command -v "${ROOT:-root}" >/dev/null; then
  printf 'ROOT is required for the legacy PDFs; load ROOT or set ROOT=/path/to/root.\n' >&2
  exit 1
fi
CHECK_ONLY=0
if [[ ${1:-} == --check ]]; then
  CHECK_ONLY=1
  shift
fi
if [[ $# -eq 0 ]]; then
  set -- HMS_5p878GeV HMS_6p117GeV HMS_6p667GeV_joint_min5
fi
settings() {
  C=${1%/}
  case "$C" in
    HMS_5p878GeV) TAG=inner10 ;;
    HMS_6p117GeV|HMS_6p667GeV_joint_min5) TAG=inner20 ;;
    *) printf 'Unsupported campaign: %s\n' "$C" >&2; return 1 ;;
  esac
  SOURCE="beam20_${TAG}"
  MODERN="${SOURCE}_holdout"
  LEGACY="${SOURCE}_legacy_holdout"
  MATRIX="$C/06e_beam_search/$SOURCE/matrices/beam.dat"
  FIT_DIR="$C/06c_core_ntuple/$TAG/holdout/root"
  METADATA="$C/05c_core_sample/$TAG/metadata/optics.dat"
}
# Check every requested campaign before evaluating any protected residuals.
for campaign in "$@"; do
  settings "$campaign"
  ./run_beam_holdout.sh "$C" "$TAG" "$MODERN" --source "$SOURCE" --check
  ./run_legacy_diagnostics.sh "$C" "$LEGACY" \
    --matrix "$MATRIX" --fit-dir "$FIT_DIR" --metadata "$METADATA" \
    --offsets 0 0 0 --dry-run
done
if [[ $CHECK_ONLY -eq 1 ]]; then
  exit 0
fi
mkdir -p batch_logs
for campaign in "$@"; do
  settings "$campaign"
  ./run_beam_holdout.sh "$C" "$TAG" "$MODERN" --source "$SOURCE" \
    2>&1 | tee -a "batch_logs/${C}_${MODERN}.log"
  ./run_legacy_diagnostics.sh "$C" "$LEGACY" \
    --matrix "$MATRIX" --fit-dir "$FIT_DIR" --metadata "$METADATA" \
    --offsets 0 0 0 \
    2>&1 | tee -a "batch_logs/${C}_${LEGACY}.log"
done
