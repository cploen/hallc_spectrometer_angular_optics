#!/usr/bin/env python3
"""Build HMS TFit trees from a saved core sample without running the SVD."""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import uproot

from core_sample import PROJECT, digest, read_tsv


def build(campaign, tag, sample):
    source = campaign / "05c_core_sample" / tag
    manifest = json.loads((source / "manifest.json").read_text())
    # Keep the build tied to the saved campaign and masks, even if config changes.
    tables = list(source.glob("rungroups_*_inputs.tsv"))
    if len(tables) != 1:
        raise ValueError("Missing saved campaign table")
    if digest(tables[0]) != manifest["outputs"][tables[0].name]:
        raise ValueError("Saved campaign table changed after selection")
    settings = read_tsv(tables[0])
    executable = shutil.which(os.environ.get("HCANA", "hcana"))
    if not executable:
        raise ValueError("hcana not found; load the campaign ROOT/HCANA environment first")
    code = {"fit": 1, "holdout": 2, "surplus": 3}[sample]
    target = campaign / "06c_core_ntuple" / tag / sample
    if target.exists():
        raise FileExistsError(f"Output exists: {target}")
    for row in settings:
        if not Path(row["rootfile"]).is_file():
            raise FileNotFoundError(row["rootfile"])
    target.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".build_", dir=target.parent))
    try:
        for sub in ("root", "plots", "logs"):
            (stage / sub).mkdir()
        report = []
        for row in settings:
            name, run = row["rungroup"], int(row["optics_id"])
            relative = f"root/CoreSample_{name}.root"
            mask = source / relative
            if digest(mask) != manifest["outputs"][relative]:
                raise ValueError(f"Core sample changed after selection: {mask}")
            with uproot.open(mask) as root:
                a = root["CoreSample"].arrays(["entry", "sample", "xscol", "yscol"], library="np")
            entries = a["entry"][a["sample"] == code]
            if np.any((a["xscol"] > 8) | (a["yscol"] > 8)):
                raise ValueError("This TFit adapter uses existing HMS 9x9 geometry; selection itself supports other hole counts")
            if not len(entries):
                report.append(dict(rungroup=name, entries=0, status="empty"))
                print(f"SKIP {name}: no {sample} events", flush=True)
                continue
            # JSON quoting here creates C++ string literals, never shell commands.
            args = [str(run), "-1", json.dumps(tag), json.dumps(str(Path(row["rootfile"]).parent)),
                    '""', '""', "6.0", "0.65", "false", '""', json.dumps(str(stage)),
                    json.dumps(row["rootfile"]), json.dumps(name), json.dumps(str(mask)), str(code)]
            expression = str(PROJECT / "make_fit_ntuple_from_gmm.C") + "(" + ",".join(args) + ")"
            logfile = stage / "logs" / f"core_{name}.log"
            print(f"Build {name}: {len(entries)} {sample} events", flush=True)
            with open(logfile, "w") as stream:
                result = subprocess.run([executable, "-b", "-l", "-q", expression], cwd=PROJECT,
                                        stdout=stream, stderr=subprocess.STDOUT)
            output = stage / "root" / f"Optics_{run}_-1_fit_tree_gmm.root"
            if result.returncode or not output.exists():
                raise RuntimeError(f"TFit build failed for {name}:\n{logfile.read_text()[-5000:]}")
            with uproot.open(output) as root:
                actual = root["TFit"]["entry"].array(library="np")
            if not np.array_equal(np.sort(entries), np.sort(actual)):
                raise RuntimeError(f"TFit membership differs from saved {sample} mask for {name}")
            report.append(dict(rungroup=name, entries=len(entries), status="verified"))
        metadata = PROJECT / "DATfiles/list_of_optics_run.dat"
        (stage / "build.json").write_text(json.dumps(dict(sample=sample, groups=report,
            sample_manifest=digest(source / "manifest.json"), metadata=digest(metadata),
            macro=digest(PROJECT / "make_fit_ntuple_from_gmm.C")), indent=2) + "\n")
        stage.rename(target)
        print(f"Core TFit build complete: {target}. SVD was not run.")
    except BaseException:
        shutil.rmtree(stage)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("tag", nargs="?", default="core")
    parser.add_argument("sample", nargs="?", choices=["fit", "holdout", "surplus"], default="fit")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.tag):
        parser.error("Invalid TAG")
    campaign = args.campaign if args.campaign.is_absolute() else PROJECT / args.campaign
    if not campaign.name.startswith("HMS_"):
        parser.error("The TFit adapter currently supports HMS campaigns only")
    try:
        build(campaign.resolve(), args.tag, args.sample)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
