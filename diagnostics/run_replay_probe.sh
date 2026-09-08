#!/usr/bin/env bash
# HMS/SHMS preparation: inspect the supplied replay paths without modifying them.
# The optional campaign name selects H/P branches; no spectrometer flag.
set -uo pipefail
if [[ $# -gt 1 ]]; then
  echo "Usage: bash diagnostics/run_replay_probe.sh [CAMPAIGN]" >&2
  exit 2
fi
REPO_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd) || exit 2
CAMPAIGN=${1:-SHMS_8p5695GeV}
REPLAY_DIR=${REPLAY_INPUT_DIR:-/volatile/hallc/c-deuteron/gvill/ROOTfiles/prod}
PROBE_MAX_EVENTS=${PROBE_MAX_EVENTS:-50000}
if [[ ! "$PROBE_MAX_EVENTS" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: PROBE_MAX_EVENTS must be a positive integer" >&2
  exit 2
fi
if command -v root >/dev/null 2>&1; then
  ROOT_BIN=root
elif command -v hcana >/dev/null 2>&1; then
  ROOT_BIN=hcana
else
  echo "ERROR: load a ROOT or HCANA environment first" >&2
  exit 2
fi
escape_cpp() {
  local value="$1"
  value=${value//\\/\\\\}
  value=${value//\"/\\\"}
  value=${value//$'\n'/\\n}
  value=${value//$'\r'/\\r}
  printf '%s' "$value"
}
macro=$(escape_cpp "$REPO_DIR/diagnostics/validation/replay/probe_replay_inputs.C")
campaign_cpp=$(escape_cpp "$CAMPAIGN")
failures=0
for run in 3283 3284 3285 3286; do
  input=$(escape_cpp "$REPLAY_DIR/deut_replay_prod_${run}_-1.root")
  echo "RUN_BEGIN $run"
  "$ROOT_BIN" -l -b -q -e "#include \"$macro\"" \
    -e "gSystem->Exit(probe_replay_inputs(\"$campaign_cpp\",\"$input\",$PROBE_MAX_EVENTS));"
  status=$?
  echo "RUN_END $run status=$status"
  if [[ $status -ne 0 ]]; then failures=$((failures+1)); fi
done
echo "PROBE_SUMMARY failed_or_needing_review=$failures total=4"
[[ $failures -eq 0 ]]
