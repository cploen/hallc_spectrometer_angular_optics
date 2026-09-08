#!/usr/bin/env bash
# Run on ifarm. Merge the four campaign inputs without modifying source files.
set -euo pipefail
cd "$(dirname -- "${BASH_SOURCE[0]}")/.."
command -v hadd >/dev/null || { echo "ERROR: hadd is required" >&2; exit 1; }
source_dir="${REPLAY_INPUT_DIR:-/volatile/hallc/c-deuteron/gvill/ROOTfiles/prod}"
output="SHMS_8p5695GeV/inputs/shms_optics_8p5695_rg01_theta8p915_foilpm10z0.root"
inputs=()
for run in 3283 3284 3285 3286; do
  input="$source_dir/deut_replay_prod_${run}_-1.root"
  [[ -r "$input" ]] || { echo "ERROR: unreadable input: $input" >&2; exit 1; }
  inputs+=("$input")
done
[[ ! -e "$output" ]] || { echo "ERROR: output already exists: $output" >&2; exit 1; }
mkdir -p "$(dirname -- "$output")"
# No -f: never overwrite an existing merged input.
hadd "$output" "${inputs[@]}"
echo "Prepared $output"
echo "Run campaign commands from the repository root."
