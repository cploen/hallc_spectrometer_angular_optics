#!/usr/bin/env bash
set -uo pipefail

usage()
{
  echo "Usage: $0 CAMPAIGN TAG [--dry-run]"
  echo
  echo "Environment overrides:"
  echo "  FILEID=-1"
  echo "  NFIT_MAX=200000"
  echo "  NSETTINGS=-1"
  echo "  RUN_SVD=1"
  echo "  OVERWRITE=0"
  echo "  VETO_FILE="
  echo "  CER_CUT=6.0"
  echo "  CAL_CUT=0.65"
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage
  exit 2
fi

PROJECT_DIR=$(pwd -P)
CAMPAIGN="$1"
TAG="$2"
MODE="${3:-}"

if [[ -n "$MODE" && "$MODE" != "--dry-run" ]]; then
  echo "ERROR: unknown argument: $MODE"
  usage
  exit 2
fi

DRYRUN=0
[[ "$MODE" == "--dry-run" ]] && DRYRUN=1

CAMPAIGN_DIR="${PROJECT_DIR}/${CAMPAIGN}"
FILEID="${FILEID:--1}"

YBASE="${CAMPAIGN_DIR}/05a_gmm_cleanup_y/root"
XBASE="${CAMPAIGN_DIR}/05b_gmm_cleanup_x/root"

NTUPLE_DIR="${CAMPAIGN_DIR}/06a_fit_ntuple"
SVD_DIR="${CAMPAIGN_DIR}/06b_svd_fit"

OLD_COEFFS="${OLD_COEFFS:-${CAMPAIGN_DIR}/config/oldfit.dat}"
VETO_FILE="${VETO_FILE:-}"

CER_CUT="${CER_CUT:-6.0}"
CAL_CUT="${CAL_CUT:-0.65}"

NFIT_MAX="${NFIT_MAX:-200000}"
NSETTINGS="${NSETTINGS:--1}"

RUN_SVD="${RUN_SVD:-1}"
OVERWRITE="${OVERWRITE:-0}"

NTUPLE_MACRO="${PROJECT_DIR}/make_fit_ntuple_from_gmm.C"
SVD_MACRO="${PROJECT_DIR}/fit_opt_matrix_gmm.C"
OPTICS_METADATA="${OPTICS_METADATA:-${PROJECT_DIR}/DATfiles/list_of_optics_run.dat}"

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

mkdir -p \
  "${NTUPLE_DIR}/root" \
  "${NTUPLE_DIR}/plots" \
  "${NTUPLE_DIR}/logs" \
  "${NTUPLE_DIR}/summaries" \
  "${SVD_DIR}/matrices" \
  "${SVD_DIR}/root" \
  "${SVD_DIR}/plots" \
  "${SVD_DIR}/logs" \
  "${SVD_DIR}/summaries"

STAMP=$(date +%Y%m%d_%H%M%S)
SUMMARY="${SVD_DIR}/summaries/build_summary_${TAG}_${STAMP}.txt"

log()
{
  printf '%s\n' "$*" | tee -a "$SUMMARY"
}

print_command()
{
  printf '  ' | tee -a "$SUMMARY"
  printf '%q ' "$@" | tee -a "$SUMMARY"
  printf '\n' | tee -a "$SUMMARY"
}

validate_tfit()
{
  local rootfile="$1"
  local output
  local status

  if [[ ! -s "$rootfile" ]]; then
    log "VALIDATION ERROR: missing or empty file: $rootfile"
    return 1
  fi

  output=$(
    root -l -b -q -e "
      TFile f(\"${rootfile}\");
      TTree *t = (TTree*)f.Get(\"TFit\");
      if (!t) {
        cout << \"No TFit tree\" << endl;
        gSystem->Exit(2);
      }
      cout << \"TFit entries \" << t->GetEntries() << endl;
      gSystem->Exit(t->GetEntries() > 0 ? 0 : 3);
    " 2>&1
  )
  status=$?

  printf '%s\n' "$output" | tee -a "$SUMMARY"
  return "$status"
}

log "============================================================"
log "GMM fit-tree and SVD build"
log "Started: $(date)"
log "============================================================"
log "CAMPAIGN      = $CAMPAIGN"
log "TAG           = $TAG"
log "RUNGROUPS_TSV = $RUNGROUPS_TSV"
log "OPTICS_METADATA = $OPTICS_METADATA"
log "FILEID        = $FILEID"
log "YBASE         = $YBASE"
log "XBASE         = $XBASE"
log "NTUPLE_DIR    = $NTUPLE_DIR"
log "SVD_DIR       = $SVD_DIR"
log "OLD_COEFFS    = $OLD_COEFFS"
log "VETO_FILE     = ${VETO_FILE:-<none>}"
log "CER_CUT       = $CER_CUT"
log "CAL_CUT       = $CAL_CUT"
log "NFIT_MAX      = $NFIT_MAX"
log "NSETTINGS     = $NSETTINGS"
log "RUN_SVD       = $RUN_SVD"
log "DRYRUN        = $DRYRUN"
log "OVERWRITE     = $OVERWRITE"
log

[[ -f "$NTUPLE_MACRO" ]] || {
  log "ERROR: missing $NTUPLE_MACRO"
  exit 1
}

[[ -f "$SVD_MACRO" ]] || {
  log "ERROR: missing $SVD_MACRO"
  exit 1
}

[[ -f "$OPTICS_METADATA" ]] || {
  log "ERROR: missing $OPTICS_METADATA"
  exit 1
}

[[ -f "$OLD_COEFFS" ]] || {
  log "ERROR: missing seed matrix $OLD_COEFFS"
  log "Place or symlink the input matrix at:"
  log "  ${CAMPAIGN_DIR}/config/oldfit.dat"
  exit 1
}

[[ -d "$YBASE" ]] || {
  log "ERROR: missing Y GMM directory: $YBASE"
  exit 1
}

[[ -d "$XBASE" ]] || {
  log "ERROR: missing X GMM directory: $XBASE"
  exit 1
}

nfailed=0
nrequested=0

while IFS=$'\t' read -r rungroup optics_id angle foils nruns runs input_root; do
  [[ -z "$rungroup" || "$rungroup" == "rungroup" ]] && continue

  ((nrequested+=1))

  if [[ -z "$optics_id" || ! "$optics_id" =~ ^[0-9]+$ ]]; then
    log "ERROR [$rungroup]: invalid optics_id '$optics_id'"
    ((nfailed+=1))
    continue
  fi

  if [[ ! -f "$input_root" ]]; then
    log "ERROR [$rungroup]: missing input ROOT file:"
    log "  $input_root"
    ((nfailed+=1))
    continue
  fi

  outfile="${NTUPLE_DIR}/root/Optics_${optics_id}_${FILEID}_fit_tree_gmm.root"
  plotfile="${NTUPLE_DIR}/plots/gmm_fit_ntuple_${optics_id}_${TAG}.pdf"
  logfile="${NTUPLE_DIR}/logs/make_fit_ntuple_${optics_id}_${TAG}_${STAMP}.log"

  log "============================================================"
  log "Building TFit tree"
  log "rungroup   = $rungroup"
  log "optics ID  = $optics_id"
  log "input ROOT = $input_root"
  log "output     = $outfile"
  log "============================================================"

  if [[ "$DRYRUN" != 1 && -s "$outfile" && "$OVERWRITE" != 1 ]]; then
    log "Existing output found; validating instead of overwriting."

    if validate_tfit "$outfile"; then
      log "VALID: existing TFit output retained."
      log
      continue
    fi

    log "Existing output is invalid; refusing to overwrite without OVERWRITE=1."
    ((nfailed+=1))
    log
    continue
  fi

  if [[ "$DRYRUN" != 1 && "$OVERWRITE" == 1 ]]; then
    rm -f -- "$outfile" "$plotfile"
  fi

  root_dir=$(dirname "$input_root")

  expr="${NTUPLE_MACRO}(${optics_id},${FILEID},\"${TAG}\",\"${root_dir}\",\"${YBASE}\",\"${XBASE}\",${CER_CUT},${CAL_CUT},true,\"${VETO_FILE}\",\"${NTUPLE_DIR}\",\"${input_root}\",\"${rungroup}\")"
  cmd=(hcana -b -l -q "$expr")

  print_command "${cmd[@]}"

  if [[ "$DRYRUN" == 1 ]]; then
    log "DRYRUN: command not executed."
    log
    continue
  fi

  "${cmd[@]}" >"$logfile" 2>&1
  hcana_status=$?

  log "hcana exit status = $hcana_status"
  log "log               = $logfile"

  if validate_tfit "$outfile"; then
    if [[ "$hcana_status" -ne 0 ]]; then
      log "WARNING: hcana returned nonzero, but the verified TFit output is valid."
    else
      log "SUCCESS: valid TFit output created."
    fi
  else
    log "FAILED: no valid TFit output."
    log "Last 40 log lines:"
    tail -40 "$logfile" | tee -a "$SUMMARY"
    ((nfailed+=1))
  fi

  log
done < "$RUNGROUPS_TSV"

if [[ "$nrequested" -eq 0 ]]; then
  log "ERROR: no rungroups were read from $RUNGROUPS_TSV"
  exit 1
fi

if [[ "$DRYRUN" == 1 ]]; then
  if [[ "$nfailed" -ne 0 ]]; then
    log "DRYRUN finished with $nfailed invalid rungroup input(s)."
    exit 1
  fi

  log "DRYRUN complete; SVD was not executed."
  exit 0
fi

if [[ "$nfailed" -ne 0 ]]; then
  log "ERROR: $nfailed requested TFit build(s) failed."
  log "SVD will not run."
  exit 1
fi

if [[ "$RUN_SVD" != 1 ]]; then
  log "RUN_SVD=$RUN_SVD; stopping after verified TFit trees."
  exit 0
fi

matrix_out="${SVD_DIR}/matrices/nps_hms_newfit_${TAG}.dat"
qa_out="${SVD_DIR}/root/fit_opt_matrix_${TAG}_qa.root"
old_pdf="${SVD_DIR}/plots/fit_opt_matrix_${TAG}_old_matrix_diff.pdf"
new_pdf="${SVD_DIR}/plots/fit_opt_matrix_${TAG}_new_matrix_diff.pdf"
svd_log="${SVD_DIR}/logs/fit_opt_matrix_${TAG}_${STAMP}.log"

if [[ "$OVERWRITE" != 1 && ( -e "$matrix_out" || -e "$qa_out" ) ]]; then
  log "ERROR: SVD output already exists."
  log "  $matrix_out"
  log "  $qa_out"
  log "Use a new TAG or set OVERWRITE=1."
  exit 1
fi

if [[ "$OVERWRITE" == 1 ]]; then
  rm -f -- "$matrix_out" "$qa_out" "$old_pdf" "$new_pdf"
fi

svd_expr="${SVD_MACRO}(\"${TAG}\",${FILEID},${NFIT_MAX},${NSETTINGS},\"${NTUPLE_DIR}/root\",\"${SVD_DIR}\",\"${OLD_COEFFS}\",\"${RUNGROUPS_TSV}\",\"${OPTICS_METADATA}\")"
svd_cmd=(hcana -b -l -q "$svd_expr")

log "============================================================"
log "Running SVD matrix fit"
log "============================================================"
print_command "${svd_cmd[@]}"

"${svd_cmd[@]}" >"$svd_log" 2>&1
svd_status=$?

log "hcana exit status = $svd_status"
log "SVD log           = $svd_log"

svd_ok=1

if [[ ! -s "$matrix_out" ]]; then
  log "ERROR: missing matrix output: $matrix_out"
  svd_ok=0
fi

if [[ ! -s "$qa_out" ]]; then
  log "ERROR: missing QA ROOT output: $qa_out"
  svd_ok=0
fi

if [[ "$svd_ok" == 1 ]]; then
  log "SUCCESS: verified SVD outputs:"
  log "  $matrix_out"
  log "  $qa_out"

  [[ -s "$old_pdf" ]] && log "  $old_pdf"
  [[ -s "$new_pdf" ]] && log "  $new_pdf"

  if [[ "$svd_status" -ne 0 ]]; then
    log "WARNING: hcana returned nonzero, but required SVD outputs exist."
  fi
else
  log "FAILED: incomplete SVD output."
  log "Last 60 SVD log lines:"
  tail -60 "$svd_log" | tee -a "$SUMMARY"
  exit 1
fi

log
log "Completed: $(date)"
log "Summary: $SUMMARY"
