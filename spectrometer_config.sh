# Shared campaign selection for runners. No new spectrometer argument.
HALLC_CONFIG_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
hallc_campaign() {
  local profile
  profile=$(python3 "$HALLC_CONFIG_DIR/spectrometer_config.py" "$1") || return 1
  IFS=$'\t' read -r HALLC_SPECTROMETER HALLC_ARM HALLC_NX HALLC_NY HALLC_DELTA_MIN HALLC_DELTA_MAX <<< "$profile"
  echo "Spectrometer: $HALLC_SPECTROMETER (centered sieve), delta=($HALLC_DELTA_MIN,$HALLC_DELTA_MAX)%"
}
