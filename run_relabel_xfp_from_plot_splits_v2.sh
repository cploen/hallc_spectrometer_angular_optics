#!/usr/bin/env bash
set -euo pipefail

METADATA_RUN=${METADATA_RUN:?ERROR: METADATA_RUN is required}
RUNGROUP=${RUNGROUP:?ERROR: RUNGROUP is required}
CAMPAIGN=${CAMPAIGN:?ERROR: CAMPAIGN is required}
INPUT_ROOT=${INPUT_ROOT:?ERROR: INPUT_ROOT is required}

MACRO=${MACRO:-relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch.C}

PLOT_DIR="${CAMPAIGN}/02b_angle_scan_x/plots"
ORIGINAL_CUT="${CAMPAIGN}/02b_angle_scan_x/cuts/XpFpXFp_${RUNGROUP}_auto_band_cut.root"
RELABEL_DIR="${CAMPAIGN}/03b_relabel_x"
FINAL_CUT="${RELABEL_DIR}/cuts/XpFpXFp_${RUNGROUP}_auto_band_cut_xscolRelabeled.root"

USE_YTAR=${USE_YTAR:-true}
MAX_EVENTS=${MAX_EVENTS:--1}
CER_CUT=${CER_CUT:-6.0}
CAL_CUT=${CAL_CUT:-0.65}
MIN_EVENTS=${MIN_EVENTS:-25}
WRITE_RELABEL=${WRITE_RELABEL:-true}
ZONE_MODE=${ZONE_MODE:-all}
WRITE_SINGLETON_ALIASES=${WRITE_SINGLETON_ALIASES:-false}
OVERRIDE_FILE=${OVERRIDE_FILE:-}

DRYRUN=${DRYRUN:-0}
KEEP_CHAIN_MACRO=${KEEP_CHAIN_MACRO:-1}

mkdir -p "${RELABEL_DIR}/cuts" "${RELABEL_DIR}/plots" "${RELABEL_DIR}/root" "${RELABEL_DIR}/logs"

[[ -f "$MACRO" ]] || { echo "ERROR: missing macro $MACRO"; exit 1; }
[[ -f "$INPUT_ROOT" ]] || { echo "ERROR: missing replay ROOT $INPUT_ROOT"; exit 1; }
[[ -f "$ORIGINAL_CUT" ]] || { echo "ERROR: missing provisional cut file $ORIGINAL_CUT"; exit 1; }
[[ -d "$PLOT_DIR" ]] || { echo "ERROR: missing plot directory $PLOT_DIR"; exit 1; }

token_to_num() {
  local tok="$1"
  tok="${tok/#m/-}"
  tok="${tok//p/.}"
  printf '%s\n' "$tok"
}

delta_tag_to_ndel() {
  case "$1" in
    m10_to_m8) echo 0 ;;
    m8_to_m5)  echo 1 ;;
    m5_to_0)   echo 2 ;;
    0_to_5)    echo 3 ;;
    5_to_10)   echo 4 ;;
    *) echo UNKNOWN ;;
  esac
}

shopt -s nullglob

declare -A latest_file
declare -A latest_time
declare -A latest_split

pattern="${PLOT_DIR}/xfp_xpfp_angleScan_${RUNGROUP}_foil*_delta_*_xpfp_*_to_999.pdf"

for f in $pattern; do
  base=$(basename "$f")

  if [[ "$base" =~ ^xfp_xpfp_angleScan_${RUNGROUP}_foil([0-9]+)_delta_(.+)_xpfp_([^_]+)_to_999\.pdf$ ]]; then
    foil=${BASH_REMATCH[1]}
    dtag=${BASH_REMATCH[2]}
    split_tok=${BASH_REMATCH[3]}
    key="${foil}|${dtag}"
    mtime=$(stat -c %Y "$f")

    if [[ -z "${latest_time[$key]+x}" || "$mtime" -gt "${latest_time[$key]}" ]]; then
      latest_time[$key]=$mtime
      latest_file[$key]=$f
      latest_split[$key]=$split_tok
    fi
  fi
done

if [[ ${#latest_file[@]} -eq 0 ]]; then
  echo "ERROR: no high-zone split PDFs found:"
  echo "  $pattern"
  exit 1
fi

rows=$(mktemp)
sorted=$(mktemp)
trap 'rm -f "$rows" "$sorted"' EXIT

for key in "${!latest_file[@]}"; do
  foil=${key%%|*}
  dtag=${key#*|}
  ndel=$(delta_tag_to_ndel "$dtag")

  if [[ "$ndel" == UNKNOWN ]]; then
    echo "WARNING: skipping unknown delta tag $dtag"
    continue
  fi

  split_tok=${latest_split[$key]}
  split=$(token_to_num "$split_tok")
  printf '%s %s %s %s %s\n' "$foil" "$ndel" "$dtag" "$split_tok" "$split" >> "$rows"
done

sort -k1,1n -k2,2n "$rows" > "$sorted"

[[ -s "$sorted" ]] || { echo "ERROR: no valid split rows found"; exit 1; }

tmpA="${RELABEL_DIR}/cuts/.tmp_xrelabel_${RUNGROUP}_A.root"
tmpB="${RELABEL_DIR}/cuts/.tmp_xrelabel_${RUNGROUP}_B.root"
chain_macro="xrelabel_chain_${RUNGROUP}.C"
nsteps=$(wc -l < "$sorted" | tr -d ' ')

echo
echo "Rungroup       : $RUNGROUP"
echo "Metadata run   : $METADATA_RUN"
echo "Replay ROOT    : $INPUT_ROOT"
echo "Provisional cut: $ORIGINAL_CUT"
echo "Relabeled cut  : $FINAL_CUT"
echo "Split steps    : $nsteps"
echo

cat "$sorted"

cat > "$chain_macro" <<EOF
#include "${MACRO}"

void xrelabel_chain_${RUNGROUP}()
{
  gROOT->SetBatch(kTRUE);

  XRELABEL_RUNGROUP_TAG = "${RUNGROUP}";
  XRELABEL_CAMPAIGN_DIR = "${CAMPAIGN}";
  XRELABEL_INPUT_TREE_OVERRIDE = "${INPUT_ROOT}";

  TString originalInput = "${ORIGINAL_CUT}";
  TString finalOutput = "${FINAL_CUT}";
  TString tmpA = "${tmpA}";
  TString tmpB = "${tmpB}";

  gSystem->Unlink(finalOutput.Data());
  gSystem->Unlink(tmpA.Data());
  gSystem->Unlink(tmpB.Data());

  struct Step { int foil; int ndel; double split; };
  Step steps[] = {
EOF

i=0
while read -r foil ndel dtag split_tok split; do
  ((i+=1))
  comma=","
  [[ $i -eq $nsteps ]] && comma=""
  printf '    {%s, %s, %.12g}%s\n' "$foil" "$ndel" "$split" "$comma" >> "$chain_macro"
done < "$sorted"

cat >> "$chain_macro" <<EOF
  };

  const int nSteps = sizeof(steps) / sizeof(steps[0]);
  TString currentInput = originalInput;

  for (int i = 0; i < nSteps; ++i) {
    TString currentOutput =
      (i == nSteps - 1) ? finalOutput : ((i % 2 == 0) ? tmpA : tmpB);

    XRELABEL_INPUT_CUT_OVERRIDE = currentInput;
    XRELABEL_OUTPUT_CUT_OVERRIDE = currentOutput;

    cout << "\\nX relabel step " << i << "/" << nSteps - 1
         << " foil=" << steps[i].foil
         << " ndel=" << steps[i].ndel
         << " split=" << steps[i].split << endl;

    relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch_oneDelta(
      ${METADATA_RUN},
      steps[i].foil,
      steps[i].ndel,
      "",
      -1,
      "${OVERRIDE_FILE}",
      ${USE_YTAR},
      ${MAX_EVENTS},
      ${CER_CUT},
      ${CAL_CUT},
      ${MIN_EVENTS},
      ${WRITE_RELABEL},
      "${ZONE_MODE}",
      ${WRITE_SINGLETON_ALIASES},
      steps[i].split
    );

    TFile* fcheck = TFile::Open(currentOutput, "READ");
    bool ok = fcheck && !fcheck->IsZombie() &&
              fcheck->GetListOfKeys() &&
              fcheck->GetListOfKeys()->GetEntries() > 0;
    if (fcheck) fcheck->Close();

    if (!ok) {
      cerr << "ERROR: no usable relabel output for foil="
           << steps[i].foil << " ndel=" << steps[i].ndel << endl;
      return;
    }

    currentInput = currentOutput;
  }

  XRELABEL_INPUT_CUT_OVERRIDE = "";
  XRELABEL_OUTPUT_CUT_OVERRIDE = "";
  XRELABEL_RUNGROUP_TAG = "";
  XRELABEL_CAMPAIGN_DIR = "";
  XRELABEL_INPUT_TREE_OVERRIDE = "";

  gSystem->Unlink(tmpA.Data());
  gSystem->Unlink(tmpB.Data());

  cout << "\\nDONE: ${FINAL_CUT}" << endl;
}
EOF

echo
echo "Generated $chain_macro"

if [[ "$DRYRUN" == 1 ]]; then
  echo "DRYRUN complete."
  exit 0
fi

logfile="${RELABEL_DIR}/logs/xrelabel_${RUNGROUP}.log"
hcana -b -l -q "$chain_macro" 2>&1 | tee "$logfile"

[[ -s "$FINAL_CUT" ]] || { echo "ERROR: final relabeled cut missing"; exit 1; }

if [[ "$KEEP_CHAIN_MACRO" != 1 ]]; then
  rm -f "$chain_macro"
fi

echo
echo "Wrote: $FINAL_CUT"
echo "Log:   $logfile"
