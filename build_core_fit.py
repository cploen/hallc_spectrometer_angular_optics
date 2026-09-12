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
from preallocated_svd import membership_manifest


def build(campaign, tag, sample):
    source = campaign / "05c_core_sample" / tag
    manifest = json.loads((source / "manifest.json").read_text())
    balanced = manifest.get('schema') == 'core_balance_v1'
    if balanced and (manifest.get('mode') != 'event_allocation' or manifest.get('status') != 'quotas_ready' or manifest['totals']['fit'] <= 0):
        raise ValueError('Export requires a positive event-level allocation; preview/zero budget cannot be exported')
    if balanced:
        for rel in ('tsv/selected_ids.tsv','metadata/optics.dat','metadata/sieve_mask.json','allocation.json'):
            if digest(source/rel) != manifest['outputs'][rel]:
                raise ValueError(f'Allocation input changed: {rel}')
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
                a = root["CoreSample"].arrays(library="np")
            entries = a["entry"][a["sample"] == code]
            if balanced and sample == 'fit':
                wanted = [r for r in read_tsv(source/'tsv/selected_ids.tsv') if r['rungroup']==name]
                expected = np.array([int(r['entry']) for r in wanted], dtype=np.int64)
                if not np.array_equal(np.sort(entries), np.sort(expected)) or len(np.unique(entries)) != len(entries):
                    raise ValueError(f'Training mask/selected IDs disagree: {name}')
                mask_selected=a['sample']==1
                if np.any(a['quality'][mask_selected]!=2) or np.any(a['balance_excluded'][mask_selected]!=0):
                    raise ValueError('Noncore or excluded event in training mask')
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
        metadata = source / "metadata/optics.dat"
        build_manifest=dict(sample=sample, groups=report,
            sample_manifest=digest(source / "manifest.json"), metadata=digest(metadata),
            macro=digest(PROJECT / "make_fit_ntuple_from_gmm.C"),
            geometry={name: digest(PROJECT / name) for name in
                      ("spectrometer_config.h", "spectrometer_root.h")})
        if balanced and sample=='fit':
            build_manifest['schema']='core_preallocated_v1'
            shutil.copy2(source/'manifest.json',stage/'allocation_manifest.json')
            shutil.copy2(source/'tsv/selected_ids.tsv',stage/'selected_ids.tsv')
            shutil.copy2(source/'metadata/optics.dat',stage/'optics.dat')
            shutil.copy2(source/'metadata/sieve_mask.json',stage/'sieve_mask.json')
            shutil.copy2(tables[0],stage/tables[0].name)
            shutil.copy2(campaign/'config/oldfit.dat',stage/'oldfit.dat')
            membership_manifest(stage/'solver_input.tsv',read_tsv(stage/'selected_ids.tsv'),settings,build_manifest['sample_manifest'])
        # Holdout exports need the same immutable file checks as training exports.
        build_manifest['outputs']={str(p.relative_to(stage)):digest(p) for p in stage.rglob('*') if p.is_file()}
        (stage/'build.json').write_text(json.dumps(build_manifest,indent=2)+'\n')
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
