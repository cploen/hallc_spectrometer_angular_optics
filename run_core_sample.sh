#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: $0 CAMPAIGN [TAG] [--check]"
  exit 2
fi

PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec "${PYTHON:-python3}" "${PROJECT_DIR}/core_sample.py" "$@"
