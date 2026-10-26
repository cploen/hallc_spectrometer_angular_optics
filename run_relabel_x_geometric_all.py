#!/usr/bin/env python3
"""Relabel unsplit geometric X bands using the existing sieve-space relabeler.

Foils and delta intervals come directly from optics metadata. The stable split
runner is retained separately; no split or angle is inferred from PDF names.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

from spectrometer_config import from_campaign, run_metadata

REPO = Path(__file__).resolve().parent


def load_groups(campaign: Path, requested=None):
    tables = sorted((campaign / "config").glob("rungroups_*_inputs.tsv"))
    if len(tables) != 1:
        raise ValueError(f"Expected one rungroups_*_inputs.tsv in {campaign / 'config'}")
    with tables[0].open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"rungroup", "optics_id", "rootfile"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing required columns in {tables[0]}")
        rows = list(reader)
    tags = [row["rungroup"] for row in rows]
    if len(tags) != len(set(tags)) or any(not re.fullmatch(r"[A-Za-z0-9_.-]+", tag) for tag in tags):
        raise ValueError("Rungroup names must be unique and contain only letters, digits, _, ., or -")
    unknown = set(requested or []) - set(tags)
    if unknown:
        raise ValueError(f"Unknown rungroup(s): {', '.join(sorted(unknown))}")
    selected = [row for row in rows if not requested or row["rungroup"] in requested]
    if not selected:
        raise ValueError("No rungroups selected")
    for row in selected:
        if not row["optics_id"].isdigit() or not row["rootfile"]:
            raise ValueError(f"Invalid optics_id or rootfile for {row['rungroup']}")
    return selected


def require_provisional_campaign(campaign: Path):
    spec = from_campaign(campaign)
    marker_path = campaign / "config/geometric_selection.json"
    if not marker_path.is_file():
        raise ValueError(f"Missing provisional campaign marker: {marker_path}; use setup_geometric_campaign.py")
    marker = json.loads(marker_path.read_text())
    if (marker.get("selection_method") != "geometric_provisional" or
            marker.get("experimental") is not True or marker.get("spectrometer") != spec.name):
        raise ValueError(f"Invalid or incompatible provisional campaign marker: {marker_path}")
    return spec


def output_paths(campaign: Path, tag: str, meta: dict):
    stage = campaign / "03b_relabel_x"
    paths = [stage / "cuts" / f"XpFpXFp_{tag}_auto_band_cut_xscolRelabeled.root"]
    for foil in range(len(meta["foils"])):
        for ndel in range(len(meta["edges"])-1):
            base = f"xfp_xpfp_xscolRelabel_{tag}_foil{foil}_ndel{ndel}"
            paths.extend([stage / "plots" / (base + ".pdf"), stage / "root" / (base + ".root")])
    return paths


def chain_source(campaign: Path, row: dict, meta: dict, staging: Path, args) -> str:
    quote = json.dumps
    original = campaign / "02b_angle_scan_x/cuts" / f"XpFpXFp_{row['rungroup']}_auto_band_cut.root"
    replay = Path(row["rootfile"]).expanduser().resolve()
    steps = [(foil, ndel) for foil in range(len(meta["foils"]))
             for ndel in range(len(meta["edges"])-1)]
    source = [f'#include {quote(str(REPO / "relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch.C"))}',
              "void xrelabel_geometric_chain() {", "  gROOT->SetBatch(kTRUE);",
              f'  XRELABEL_RUNGROUP_TAG = {quote(row["rungroup"])};',
              f'  XRELABEL_CAMPAIGN_DIR = {quote(str(campaign))};',
              f'  XRELABEL_INPUT_TREE_OVERRIDE = {quote(str(replay))};',
              f'  TString currentInput = {quote(str(original))};']
    for index, (foil, ndel) in enumerate(steps):
        destination = staging / ("final.root" if index == len(steps)-1 else f"part{index % 2}.root")
        source.extend([
            f'  XRELABEL_INPUT_CUT_OVERRIDE = currentInput;',
            f'  XRELABEL_OUTPUT_CUT_OVERRIDE = {quote(str(destination))};',
            "  relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch_oneDelta(" +
            f'{int(row["optics_id"])}, {foil}, {ndel}, "", -1, {quote(str(args.override_file) if args.override_file else "")}, ' +
            f'true, {args.max_events}, {args.cer_cut:.17g}, {args.cal_cut:.17g}, {args.min_events}, true, "global", false, -999.0);',
            "  {",
            f'    TFile check({quote(str(destination))}, "READ");',
            "    int nRelabeled = 0;",
            "    if (!check.IsZombie()) {",
            "      TIter next(check.GetListOfKeys()); TKey* key = nullptr;",
            "      while ((key = (TKey*)next())) { TString name = key->GetName();",
            f'        if (name.BeginsWith("hXpFpXFp_cut_xscol_") && name.Contains("_nfoil_{foil}_ndel_{ndel}_part")) ++nRelabeled;',
            "      }",
            "    }",
            f'    if (nRelabeled == 0) {{ std::cerr << "ERROR: no relabeled global cuts for foil {foil}, ndel {ndel}" << std::endl; gSystem->Exit(1); return; }}',
            "  }",
            f'  currentInput = {quote(str(destination))};',
        ])
    source.extend(['  XRELABEL_INPUT_CUT_OVERRIDE = "";', '  XRELABEL_OUTPUT_CUT_OVERRIDE = "";',
                   '  XRELABEL_RUNGROUP_TAG = "";', '  XRELABEL_CAMPAIGN_DIR = "";',
                   '  XRELABEL_INPUT_TREE_OVERRIDE = "";', "}", ""])
    return "\n".join(source)


def run_group(campaign: Path, row: dict, meta: dict, args) -> None:
    tag = row["rungroup"]
    paths = output_paths(campaign, tag, meta)
    existing = [path for path in paths if path.exists() or path.is_symlink()]
    if existing and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {existing[0]}; use a fresh campaign or explicit --overwrite")
    original = campaign / "02b_angle_scan_x/cuts" / f"XpFpXFp_{tag}_auto_band_cut.root"
    replay = Path(row["rootfile"]).expanduser().resolve()
    if not args.dry_run:
        for path in (original, replay):
            if not path.is_file():
                raise FileNotFoundError(f"Missing input: {path}")
        foil_cuts = campaign / "01_ytar_cuts/cuts" / f"ytar_ridge_cut_{tag}.root"
        if not foil_cuts.is_file():
            raise FileNotFoundError(f"Missing foil cuts: {foil_cuts}")
    print(f"{'PLAN' if args.dry_run else 'RUN'} [{tag}] {len(meta['foils'])} foils, {len(meta['edges'])-1} delta intervals; global X bands", flush=True)
    for foil in range(len(meta["foils"])):
        for ndel, (low, high) in enumerate(zip(meta["edges"], meta["edges"][1:])):
            print(f"  foil={foil} z={meta['foils'][foil]:g} cm ndel={ndel} delta=[{low:g},{high:g})", flush=True)
    stage = campaign / "03b_relabel_x"
    if args.dry_run:
        print("  " + shlex.join(shlex.split(args.root_command) + ["-b", "-l", "-q", str(stage / "logs/xrelabel_geometric_chain.C")]))
        print(f"  output: {paths[0]}")
        return
    for folder in ("cuts", "plots", "root", "logs"):
        (stage / folder).mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    logfile = stage / "logs" / f"xrelabel_geometric_{tag}_{stamp}.log"
    saved_macro = stage / "logs" / f"xrelabel_geometric_{tag}_{stamp}.C"
    with tempfile.TemporaryDirectory(prefix=".geometric-xrelabel-", dir=stage / "cuts") as folder:
        staging = Path(folder)
        # The legacy macro reads this fixed metadata path; isolate it per job.
        (staging / "DATfiles").mkdir()
        (staging / "DATfiles/list_of_optics_run.dat").symlink_to(args.metadata)
        macro = staging / "xrelabel_geometric_chain.C"
        source = chain_source(campaign, row, meta, staging, args)
        macro.write_text(source)
        saved_macro.write_text(source)
        command = shlex.split(args.root_command) + ["-b", "-l", "-q", str(macro)]
        print("  " + shlex.join(command), flush=True)
        with logfile.open("x") as stream:
            result = subprocess.run(command, cwd=staging, stdout=stream, stderr=subprocess.STDOUT)
        final = staging / "final.root"
        log = logfile.read_text(errors="replace")
        if result.returncode or "ERROR:" in log or not final.is_file() or final.stat().st_size == 0:
            raise RuntimeError(f"Relabeling failed for {tag}; inspect {logfile}")
        if args.overwrite:
            os.replace(final, paths[0])
        else:
            # Atomic exclusive publication prevents accidental replacement.
            os.link(final, paths[0])
        print(f"OK [{tag}] {paths[0]}", flush=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--rungroup", action="append", help="Process only this rungroup; repeatable")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite", action="store_true", help="Explicitly replace this alternative campaign's relabel outputs")
    parser.add_argument("--metadata", type=Path, default=REPO / "DATfiles/list_of_optics_run.dat")
    parser.add_argument("--root-command", default="hcana", help="ROOT-compatible executable, default hcana")
    parser.add_argument("--override-file", type=Path, default=os.environ.get("OVERRIDE_FILE") or None)
    parser.add_argument("--min-events", type=int, default=int(os.environ.get("MIN_EVENTS", "25")))
    parser.add_argument("--max-events", type=int, default=int(os.environ.get("MAX_EVENTS", "-1")))
    parser.add_argument("--cer-cut", type=float, default=float(os.environ.get("CER_CUT", "6.0")))
    parser.add_argument("--cal-cut", type=float, default=float(os.environ.get("CAL_CUT", "0.65")))
    args = parser.parse_args(argv)
    args.campaign = args.campaign.expanduser().resolve()
    args.metadata = args.metadata.expanduser().resolve()
    if args.override_file:
        args.override_file = Path(args.override_file).expanduser().resolve()
    try:
        spec = require_provisional_campaign(args.campaign)
        if not args.metadata.is_file():
            raise FileNotFoundError(f"Missing optics metadata: {args.metadata}")
        if args.override_file and not args.override_file.is_file():
            raise FileNotFoundError(f"Missing label override file: {args.override_file}")
        if args.min_events < 1 or args.max_events == 0 or args.max_events < -1:
            raise ValueError("min-events must be positive; max-events must be -1 or positive")
        groups = load_groups(args.campaign, args.rungroup)
        # Preflight every requested group before writing any output.
        plan = []
        for row in groups:
            meta = run_metadata(int(row["optics_id"]), args.metadata)
            if spec.name == "SHMS" and meta["sieve_flag"] != 1:
                raise ValueError("Only the centered SHMS sieve (SieveFlag=1) is supported")
            existing = [path for path in output_paths(args.campaign, row["rungroup"], meta)
                        if path.exists() or path.is_symlink()]
            if existing and not args.overwrite:
                raise FileExistsError(f"Refusing to overwrite {existing[0]}; use a fresh campaign or explicit --overwrite")
            plan.append((row, meta))
        for row, meta in plan:
            run_group(args.campaign, row, meta, args)
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
