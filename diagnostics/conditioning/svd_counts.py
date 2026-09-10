#!/usr/bin/env python3
"""Recover historical fit counts in tree order and verify the saved fit log."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import ROOT
from preliminary_angular_conditioning import read_campaign_rows, read_optics_metadata, selected_indices

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from spectrometer_config import from_campaign


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("campaign", type=Path)
    p.add_argument("log", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--metadata", type=Path, default=PROJECT / "DATfiles/list_of_optics_run.dat")
    args = p.parse_args()
    campaign = args.campaign
    table, = campaign.glob("config/rungroups_*_inputs.tsv")
    settings = read_campaign_rows(table)
    metadata_path = args.metadata
    metadata = read_optics_metadata(metadata_path, {i for _, i in settings})
    spec = from_campaign(campaign)
    log = args.log.read_text()
    expected_total = int(re.search(r"number to fit =\s*(\d+)", log)[1])
    expected_y = {}
    run = foil = delta = None
    for line in log.splitlines():
        if m := re.search(r"INfile = .*Optics_(\d+)_-1_fit_tree", line):
            run = int(m[1])
        elif m := re.match(r"\s*ztar =\s*([-\d.]+)", line):
            foil = metadata[run][0].index(float(m[1]))
            delta = -1
        elif "Ndelta =" in line:
            delta += 1
        elif re.fullmatch(r"\s*\d+(?:\s+\d+){8}\s*", line) and run is not None:
            for iy, value in enumerate(map(int, line.split())):
                expected_y[(run, foil, delta, iy)] = value
    records, hashes, selected_ids = [], {}, {}
    total = 0
    # Check original numerical inputs against the retained reproduction record.
    reproduction = (campaign / "08_preliminary_conditioning/REPRODUCTION.md").read_text()
    recorded_hashes = {path: digest for digest, path in re.findall(r"^([0-9a-f]{64})  (.+)$", reproduction, re.M)}
    for path in [table, metadata_path]:
        digest = sha(path)
        rel = "DATfiles/list_of_optics_run.dat" if path == metadata_path else (str(path.relative_to(PROJECT)) if path.is_absolute() else str(path))
        assert recorded_hashes[rel] == digest, f"Historical input changed: {path}"
        hashes[rel] = digest
    actual_y = {}
    for name, run in settings:
        path = campaign / f"06a_fit_ntuple/root/Optics_{run}_-1_fit_tree_gmm.root"
        digest = sha(path)
        rel = str(path.relative_to(PROJECT)) if path.is_absolute() else str(path)
        assert recorded_hashes[rel] == digest, f"Historical input changed: {path}"
        hashes[rel] = digest
        a = ROOT.RDataFrame("TFit", str(path)).AsNumpy(
            ["entry", "foil", "ndel", "xscol", "yscol", "delta", "ztarT", "ysT"])
        indices, per_foil = selected_indices(a, *metadata[run], total, 200000, spec)
        total += len(indices)
        selected_ids[f"{run}_tree_index"] = indices
        selected_ids[f"{run}_entry"] = a["entry"][indices]
        available = Counter(zip(a["foil"], a["ndel"], a["xscol"], a["yscol"]))
        retained = Counter(zip(a["foil"][indices], a["ndel"][indices], a["xscol"][indices], a["yscol"][indices]))
        for (nf, nd, ix, iy), n in sorted(available.items()):
            fit = retained[(nf, nd, ix, iy)]
            key = (run, int(nf), int(nd), int(iy))
            actual_y[key] = actual_y.get(key, 0) + fit
            records.append(dict(rungroup=name, foil=int(nf), ndel=int(nd), xscol=int(ix),
                                yscol=int(iy), available=n, training=fit))
        print(f"{name}: {len(a['entry']):,} available; {len(indices):,} selected; foil counts {per_foil}")
    assert total == expected_total, (total, expected_total)
    assert expected_y, "No detailed counts found in log"
    assert all(actual_y.get(k, 0) == n for k, n in expected_y.items()), "Log cell counts differ"
    assert all(k in expected_y or n == 0 for k, n in actual_y.items()), "Unlogged selected cell"
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "optics.dat").write_bytes(metadata_path.read_bytes())
    with (args.out / "counts.tsv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(records)
    np.savez_compressed(args.out / "selected.npz", **selected_ids)
    (args.out / "provenance.json").write_text(json.dumps(dict(
        fit_log=str(args.log), log_sha256=sha(args.log), input_sha256=hashes,
        script_sha256=sha(Path(__file__)), selector_sha256=sha(Path(__file__).with_name("preliminary_angular_conditioning.py")),
        selected=total, available=sum(r["available"] for r in records),
        matched_log_cells=len(expected_y),
        method="Sequential TFit order; 1000 per foil/delta/Y column, 15000 per rungroup/foil, 200000 global",
        limitation="Reconstructed membership from historical inputs and selection code; no original event-ID list was saved."
    ), indent=2) + "\n")
    print(f"Verified {total:,} selected events and all {len(expected_y)} logged Y-column counts.")


if __name__ == "__main__":
    main()
