#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
cd "$PROJECT_DIR"
CHECK_ONLY=0
if [[ ${1:-} == --check ]]; then CHECK_ONLY=1; shift; fi
if [[ $# -eq 0 ]]; then set -- HMS_5p878GeV HMS_6p117GeV HMS_6p667GeV_joint_min5; fi
settings() {
  C=${1%/}
  case "$C" in
    HMS_5p878GeV) TAG=inner10 ;;
    HMS_6p117GeV|HMS_6p667GeV_joint_min5) TAG=inner20 ;;
    *) printf 'Unsupported campaign: %s\n' "$C" >&2; return 1 ;;
  esac
  SOURCE="beam20_${TAG}"
  NAME="${SOURCE}_standard"
}
for campaign in "$@"; do
  settings "$campaign"
  "${PYTHON:-python3}" matrix_diagnostic_suite.py "$C" "$TAG" "$NAME" --source "$SOURCE" --check
done
if [[ $CHECK_ONLY -eq 1 ]]; then exit 0; fi
mkdir -p batch_logs
for campaign in "$@"; do
  settings "$campaign"
  "${PYTHON:-python3}" matrix_diagnostic_suite.py "$C" "$TAG" "$NAME" --source "$SOURCE" \
    2>&1 | tee -a "batch_logs/${C}_${NAME}.log"
done
