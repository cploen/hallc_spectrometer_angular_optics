#!/usr/bin/env python3
"""Run one stage of a provisional geometric-selection campaign.

The driver retains the stable ridge, candidate, GMM, fit-ntuple and SVD runners.
Run stages explicitly, inspecting the foil cuts, selection and physical labels
between stages. --plan and --dry-run only print commands and write no files.

A normal campaign runtime needs HCANA/ROOT and a Python installation with PyROOT,
NumPy, SciPy, Matplotlib and scikit-learn. In particular, the stable GMM scripts
need PyROOT and scikit-learn in the SAME Python environment. --root-python may
supply separate ROOT I/O for the new selector; it does not replace that GMM
runtime requirement.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

from run_relabel_x_geometric_all import load_groups, require_provisional_campaign
from spectrometer_config import run_metadata

REPO = Path(__file__).resolve().parent
STAGES = (
    "ridge", "selection", "relabel-y", "relabel-x", "candidates-y",
    "candidates-x", "gmm-y", "gmm-x", "fit",
)
STAGE_DIRS = {
    "ridge": ("01_ytar_cuts",),
    "relabel-y": ("03a_relabel_y",),
    "relabel-x": ("03b_relabel_x",),
    "candidates-y": ("04a_candidate_trees_y",),
    "candidates-x": ("04b_candidate_trees_x",),
    "gmm-y": ("05a_gmm_cleanup_y",),
    "gmm-x": ("05b_gmm_cleanup_x",),
    "fit": ("06a_fit_ntuple", "06b_svd_fit"),
}
WRAPPERS = {
    "ridge": "run_ytar_ridge_all_rungroups.sh",
    "candidates-y": "run_make_yscol_candidates_all.sh",
    "candidates-x": "run_make_xscol_candidates_all.sh",
    "gmm-y": "run_gmm_cleanup_y_all.sh",
    "gmm-x": "run_gmm_cleanup_x_all.sh",
    "fit": "run_build_gmm_fit_and_matrix.sh",
}


def resolve_python(value: str) -> Path:
    executable = shutil.which(value)
    if executable is None:
        candidate = Path(value).expanduser()
        if not candidate.is_file() or not os.access(candidate, os.X_OK):
            raise ValueError(f"Python interpreter is not executable: {value}")
        executable = str(candidate)
    # Keep the selected environment's bin directory, including virtualenv links.
    return Path(executable).absolute()


def stage_directories(stage: str, axis: str) -> tuple[str, ...]:
    if stage == "selection":
        return tuple(name for selected, name in (("y", "02a_angle_scan_y"), ("x", "02b_angle_scan_x"))
                     if axis in (selected, "both"))
    return STAGE_DIRS[stage]


def existing_stage_output(campaign: Path, stage: str, axis: str) -> Path | None:
    for folder in stage_directories(stage, axis):
        directory = campaign / folder
        if directory.is_symlink() or (directory.exists() and not directory.is_dir()):
            return directory
        if directory.is_dir():
            for path in directory.rglob("*"):
                if path.is_file() or path.is_symlink():
                    return path
    return None


def command_for(stage: str, campaign: Path, args, python: Path) -> list[str]:
    # Stable shell runners prepend their working directory to CAMPAIGN.
    relative = os.path.relpath(campaign, REPO)
    if stage == "selection":
        command = [str(python), str(REPO / "run_geometric_selection.py"), relative, "--axis", args.axis]
        if args.root_python:
            command += ["--root-python", args.root_python]
        return command
    if stage in ("relabel-x", "relabel-y"):
        command = [str(python), str(REPO / f"run_relabel_{stage[-1]}_geometric_all.py"), relative]
        if args.label_overrides:
            command += ["--override-file", str(args.label_overrides)]
        return command
    command = ["./" + WRAPPERS[stage], relative]
    if stage == "fit":
        command.append(args.tag or "geometric_trial")
    return command


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing input: {path}")


def check_stage_inputs(stage: str, campaign: Path, groups: list[dict]) -> None:
    """Check immediate inputs only when executing; planning needs no farm mount."""
    for row in groups:
        tag = row["rungroup"]
        if stage not in ("gmm-x", "gmm-y"):
            replay = Path(row["rootfile"]).expanduser()
            require_file(replay if replay.is_absolute() else REPO / replay)
        if stage in ("selection", "relabel-y", "relabel-x", "candidates-y", "candidates-x"):
            require_file(campaign / "01_ytar_cuts/cuts" / f"ytar_ridge_cut_{tag}.root")
        if stage == "relabel-y":
            require_file(campaign / "02a_angle_scan_y/cuts" / f"YpFpYFp_{tag}_auto_band_cut.root")
        elif stage == "relabel-x":
            require_file(campaign / "02b_angle_scan_x/cuts" / f"XpFpXFp_{tag}_auto_band_cut.root")
        elif stage == "candidates-y":
            require_file(campaign / "03a_relabel_y/cuts" / f"YpFpYFp_{tag}_auto_band_cut_yscolRelabeled.root")
        elif stage == "candidates-x":
            require_file(campaign / "03b_relabel_x/cuts" / f"XpFpXFp_{tag}_auto_band_cut_xscolRelabeled.root")
        elif stage in ("gmm-y", "gmm-x"):
            axis = stage[-1]
            candidate = "YscolCandidates" if axis == "y" else "XscolCandidates"
            directory = campaign / ("04a_candidate_trees_y" if axis == "y" else "04b_candidate_trees_x")
            require_file(directory / "root" / f"{candidate}_{tag}.root")
            require_file(directory / "tsv" / f"{candidate}_{tag}_summary.tsv")
    if stage == "fit":
        require_file(campaign / "config/oldfit.dat")
        for folder in ("05a_gmm_cleanup_y/root", "05b_gmm_cleanup_x/root"):
            directory = campaign / folder
            if not directory.is_dir() or not any(directory.glob("gmm_clean_*.root")):
                raise FileNotFoundError(f"Missing GMM mask inputs: {directory}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("stage", nargs="?", choices=(*STAGES, "plan"), help="Run exactly one stage, or print the full plan")
    parser.add_argument("--plan", action="store_true", help="Print the full ordered campaign commands without executing")
    parser.add_argument("--dry-run", action="store_true", help="Print the selected stage command without executing")
    parser.add_argument("--axis", choices=("x", "y", "both"), default="both", help="Selection axes (default both)")
    parser.add_argument("--python", default=sys.executable, help="Python interpreter for selector and stable runner Python calls")
    parser.add_argument("--root-python", default=os.environ.get("HALLC_ROOT_PYTHON"), help="Optional separate PyROOT interpreter for selector I/O")
    parser.add_argument("--tag", help="Required for fit; also sets the fit tag displayed by --plan")
    parser.add_argument("--label-overrides", type=Path,
                        help="Physical-label override TSV for the selected relabel-x or relabel-y stage")
    args = parser.parse_args(argv)
    if args.stage is None and not args.plan:
        parser.error("choose one explicit stage or --plan")
    if args.plan and args.stage not in (None, "plan"):
        parser.error("--plan cannot be combined with an execution stage")
    planning = args.plan or args.stage == "plan"
    if args.stage == "fit" and not args.tag:
        parser.error("fit requires a descriptive --tag")
    if args.tag and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.tag):
        parser.error("tag must begin with a letter or digit and contain only letters, digits, _, ., or -")
    if args.label_overrides and (planning or args.stage not in ("relabel-x", "relabel-y")):
        parser.error("--label-overrides applies to one explicit relabel-x or relabel-y stage")
    try:
        campaign = args.campaign.expanduser().resolve()
        spec = require_provisional_campaign(campaign)
        groups = load_groups(campaign)
        metadata = REPO / "DATfiles/list_of_optics_run.dat"
        require_file(metadata)
        for row in groups:
            meta = run_metadata(int(row["optics_id"]), metadata)
            if spec.name == "SHMS" and meta["sieve_flag"] != 1:
                raise ValueError("Only the centered SHMS sieve (SieveFlag=1) is supported")
        python = resolve_python(args.python)
        if args.root_python:
            args.root_python = str(resolve_python(args.root_python))
        if args.label_overrides:
            args.label_overrides = args.label_overrides.expanduser().resolve()
            require_file(args.label_overrides)
        stages = STAGES if planning else (args.stage,)
        commands = [(stage, command_for(stage, campaign, args, python)) for stage in stages]
        for stage, command in commands:
            require_file(REPO / (WRAPPERS[stage] if stage in WRAPPERS else
                                "run_geometric_selection.py" if stage == "selection" else
                                f"run_relabel_{stage[-1]}_geometric_all.py"))
            if not planning and not args.dry_run:
                existing = existing_stage_output(campaign, stage, args.axis)
                if existing is not None:
                    raise FileExistsError(f"Refusing to overwrite existing {stage} output: {existing}; use a fresh campaign")
                check_stage_inputs(stage, campaign, groups)
        env = os.environ.copy()
        env["PATH"] = str(python.parent) + os.pathsep + env.get("PATH", "")
        env["OVERWRITE"] = "0"
        if args.root_python:
            env["HALLC_ROOT_PYTHON"] = args.root_python
        if "fit" in stages:
            env["RUN_SVD"] = "1"
        print(f"Campaign: {campaign} ({spec.name}, provisional geometric selection)", flush=True)
        print(f"Working directory: {REPO}", flush=True)
        print(f"Python: {python}; its bin directory is first on PATH", flush=True)
        print("export PATH=" + shlex.quote(str(python.parent)) + ':"$PATH"', flush=True)
        print("export OVERWRITE=0", flush=True)
        if planning and not args.tag:
            print("Plan uses example fit tag geometric_trial; choose --tag when executing fit.", flush=True)
        for stage, command in commands:
            shown = ["env", "RUN_SVD=1"] if stage == "fit" else []
            print(f"[{stage}] " + shlex.join(shown + command), flush=True)
            if not planning and not args.dry_run:
                subprocess.run(command, cwd=REPO, env=env, check=True)
        return 0
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
