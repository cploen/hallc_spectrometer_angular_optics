#!/usr/bin/env python3
"""Allocate unit-weight cores plus 10% or 20% inner shoulders from a frozen tag.

Copy this file to the existing joint-surface checkout. --check writes nothing.
The default remains 20%; select 90/10 explicitly with --shoulder-percent 10.
The output is a mixed allocation mask tag, not a legacy core-only tag.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile

import numpy as np
import uproot
from core_balance import LEAF, allocate_counts
from core_sample import PROJECT, digest, groups, ranks, write_tsv
from reallocate_core import checked, event_allocate, load_inputs

SCHEMA = "inner20_allocation_v1"


def capped_proportional(capacity, target, budget):
    """Study allocator: minimum-one coverage, capped proportions, integer total."""
    raw = np.asarray(capacity)
    target = np.asarray(target, dtype=float)
    if (raw.ndim != 1 or target.shape != raw.shape or not np.isfinite(raw).all()
            or np.any(raw < 0) or np.any(raw != np.floor(raw))
            or not np.isfinite(target).all() or np.any(target < 0)):
        raise ValueError("Invalid capacities or target proportions")
    capacity = raw.astype(np.int64)
    if type(budget) is not int or not 0 <= budget <= int(capacity.sum()):
        raise ValueError("Requested population budget exceeds capacity")
    if not budget:
        return np.zeros_like(capacity)
    minimum = ((capacity > 0) & (budget >= np.count_nonzero(capacity))).astype(np.int64)
    quota = minimum.astype(float)
    room = (capacity - minimum).astype(float)
    remaining = float(budget - int(minimum.sum()))
    while remaining > 1e-10:
        live = room > 1e-12
        if not live.any():
            raise ArithmeticError("Capacity exhausted before reaching budget")
        proportion = np.where(live, target, 0.)
        if proportion.sum() == 0:
            proportion = room.copy()
        proposed = remaining * proportion / proportion.sum()
        capped = live & (proposed >= room - 1e-12)
        if not capped.any():
            quota += proposed
            break
        quota[capped] += room[capped]
        remaining -= float(room[capped].sum())
        room[capped] = 0.
    result = np.floor(quota + 1e-12).astype(np.int64)
    missing = budget - int(result.sum())
    if missing < 0:
        raise ArithmeticError("Floating allocation exceeded budget")
    if missing:
        eligible = np.flatnonzero(result < capacity)
        order = np.lexsort((eligible, -(quota[eligible] - result[eligible])))
        if missing > len(eligible):
            raise ArithmeticError("Integer allocation lost too much mass")
        result[eligible[order[:missing]]] += 1
    if int(result.sum()) != budget or np.any(result > capacity) or np.any(result < minimum):
        raise ArithmeticError("Invalid integer allocation")
    return result


def select_mixture(data, baseline, cfg, foils, per_foil=None, shoulder_percent=20):
    """Use reconstructed baseline core IDs and frozen quality-1 density scores.

    The caller already ran event_allocate on fresh arrays. No arrays are changed
    until all quotas and combined baseline-cell coverage have passed checks.
    """
    if type(shoulder_percent) is not int or shoulder_percent not in (10, 20):
        raise ValueError("Shoulder percentage must be 10 or 20")
    inner_parts = shoulder_percent // 10
    core_parts = 10 - inner_parts
    leaves = baseline["leaves"]
    index = {tuple(r[k] for k in LEAF): i for i, r in enumerate(leaves)}
    target = np.asarray([r["quota"] for r in leaves], dtype=np.int64)
    capacity = dict(core=np.zeros(len(leaves), dtype=np.int64),
                    inner=np.zeros(len(leaves), dtype=np.int64))
    pools, cell_events, protected = {}, {}, {}
    for name, a in data.items():
        protected[name] = a["sample"] == 2
        development = ~protected[name] & (a["balance_excluded"] == 0)
        core = a["sample"] == 1
        if np.any(core & (~development | (a["quality"] != 2))):
            raise ValueError("Invalid baseline core population")
        inner = (development & (a["quality"] == 1) & np.isfinite(a["core_score"])
                 & (a["core_score"] >= .4) & (a["core_score"] < .5))
        pools[name] = dict(core=core, inner=inner)
        for key, ids in groups(a, list(LEAF[1:])).items():
            cell = index[(name, *key)]
            cell_events[cell] = (name, ids)
            for population in capacity:
                capacity[population][cell] = int(pools[name][population][ids].sum())
    if len(cell_events) != len(leaves) or not np.array_equal(capacity["core"], target):
        raise ValueError("Baseline event counts do not match the saved quotas")
    z = np.asarray([r["zfoil"] for r in leaves])
    ceilings = []
    for foil in foils:
        h = z == foil
        ceilings.append(10 * min(int(capacity["core"][h].sum()) // core_parts,
                                 int(capacity["inner"][h].sum()) // inner_parts))
    ceiling = min([10 * (cfg["fit_max"] // len(foils) // 10)] + ceilings)
    budget = ceiling if per_foil is None else per_foil
    if type(budget) is not int or budget <= 0 or budget % 10 or budget > ceiling:
        raise ValueError(f"Per-foil budget must be a positive multiple of ten <= {ceiling}")
    quotas = {p: np.zeros(len(leaves), dtype=np.int64) for p in capacity}
    for foil in foils:
        h = z == foil
        quotas["core"][h] = capped_proportional(capacity["core"][h], target[h], core_parts * budget // 10)
        quotas["inner"][h] = capped_proportional(capacity["inner"][h], target[h], inner_parts * budget // 10)
    uncovered = (target > 0) & ((quotas["core"] + quotas["inner"]) == 0)
    if uncovered.any():
        examples = [tuple(leaves[i][k] for k in LEAF) for i in np.flatnonzero(uncovered)[:8]]
        raise ValueError(f"Budget loses {int(uncovered.sum())} baseline cells: {examples}")
    selected = {name: np.zeros(len(a["entry"]), dtype=bool) for name, a in data.items()}
    for cell, (name, ids) in cell_events.items():
        a = data[name]
        for population in capacity:
            pool = ids[pools[name][population][ids]]
            order = np.lexsort((a["entry"][pool], ranks(cfg["seed"], name, a["entry"][pool],
                                                       "literal_counts_" + ("core" if population == "core" else "shoulder"))))
            selected[name][pool[order[:quotas[population][cell]]]] = True
    chosen, by_foil = [], []
    for name, a in data.items():
        keep = selected[name]
        if np.any(keep & ~(pools[name]["core"] | pools[name]["inner"])) or np.any(keep & protected[name]):
            raise AssertionError("Reserved or ineligible event selected")
        a["sample"][pools[name]["core"]] = 3
        a["sample"][keep] = 1
        a[f"inner{shoulder_percent}_keep"] = keep.astype(np.int32)
        a["training_population"] = np.where(keep, a["quality"], 0).astype(np.int32)
        if not np.array_equal(a["sample"] == 2, protected[name]):
            raise AssertionError("Reserved identities changed")
        for i in np.flatnonzero(keep):
            row = dict(rungroup=name, entry=int(a["entry"][i]), foil=int(a["foil"][i]), ndel=int(a["ndel"][i]),
                       zfoil=float(a["zfoil"][i]), delta_low=float(a["delta_low"][i]), delta_high=float(a["delta_high"][i]),
                       xscol=int(a["xscol"][i]), yscol=int(a["yscol"][i]), quality=int(a["quality"][i]),
                       population="core" if a["quality"][i] == 2 else "inner_shoulder", core_score=float(a["core_score"][i]))
            chosen.append(row)
    for foil in foils:
        ncore = sum(int(np.sum((a["zfoil"] == foil) & (a["sample"] == 1) & (a["quality"] == 2))) for a in data.values())
        ninner = sum(int(np.sum((a["zfoil"] == foil) & (a["sample"] == 1) & (a["quality"] == 1))) for a in data.values())
        if (ncore, ninner) != (core_parts * budget // 10, inner_parts * budget // 10):
            raise AssertionError("Per-foil mixture counts disagree")
        h = z == foil
        by_foil.append(dict(zfoil=float(foil), total=budget, core=ncore, inner_shoulder=ninner,
                            baseline_cells=int(np.sum(h & (target > 0))),
                            baseline_cells_without_inner=int(np.sum(h & (target > 0) & (capacity["inner"] == 0)))))
    chosen.sort(key=lambda row: (row["rungroup"], row["entry"]))
    for i, row in enumerate(leaves):
        row.update(baseline_quota=int(target[i]), core_available=int(capacity["core"][i]),
                   inner_available=int(capacity["inner"][i]), core_quota=int(quotas["core"][i]),
                   inner_quota=int(quotas["inner"][i]), mixed_quota=int(quotas["core"][i] + quotas["inner"][i]))
    return chosen, dict(foil_budget=budget, count_ceiling=ceiling, training_total=len(chosen), shoulder_percent=shoulder_percent,
                       foils=by_foil, leaves=leaves, baseline_cells_preserved=True,
                       reserved_identities_preserved=True)


def prepare(campaign, parent, per_foil=None, shoulder_percent=20):
    inp = load_inputs(campaign, parent)
    if inp["manifest"].get("schema") == SCHEMA:
        raise ValueError("Use the original frozen density parent, not a mixed tag")
    if inp["cfg"]["min_events"] != 10 or inp["cfg"]["core"] != .5:
        raise ValueError("Expected frozen min_events=10 and core=0.50")
    active, = (campaign / "config").glob("rungroups_*_inputs.tsv")
    if digest(active) != digest(inp["table"]):
        raise ValueError("Density parent table differs from active campaign; use min10_7008 for the expanded 5.878 pass")
    baseline = allocate_counts(inp["rows"], inp["cfg"], inp["required"], zip(inp["edges"], inp["edges"][1:]))
    data, hashes = {}, {str(inp["source"] / "manifest.json"): digest(inp["source"] / "manifest.json")}
    for path in (inp["table"], inp["metadata"], inp["maskpath"],
                 inp["source"] / "config.json", inp["source"] / "tsv/allocation.tsv", inp["source"] / "tsv/counts.tsv"):
        hashes[str(path)] = digest(path)
    policy = campaign / "config/core_balance.json"
    if policy.exists():
        hashes[str(policy)] = digest(policy)
    for row in inp["settings"]:
        relative = f"root/CoreSample_{row['rungroup']}.root"
        path = checked(inp["source"] / relative, inp["manifest"]["outputs"][relative])
        hashes[str(path)] = digest(path)
        with uproot.open(path) as root:
            data[row["rungroup"]] = root["CoreSample"].arrays(library="np")
        a = data[row["rungroup"]]
        if np.any(a["run"] != int(row["optics_id"])) or np.any(a["entry"] < 0):
            raise ValueError(f"Invalid frozen run/entry identity: {row['rungroup']}")
        if not np.array_equal(a["core_keep"], a["quality"] == 2):
            raise ValueError(f"Frozen core_keep and quality disagree: {row['rungroup']}")
        if "core_score" not in a:
            raise ValueError("Frozen tag lacks density scores; do not substitute a radius cut")
    baseline_ids, geometry_audit = event_allocate(data, inp, baseline)
    # Keep the original baseline report separate from the mixed quota annotations.
    baseline_report = json.loads(json.dumps(baseline))
    chosen, result = select_mixture(data, baseline, inp["cfg"], inp["required"], per_foil, shoulder_percent)
    inp.update(data=data, hashes=hashes, baseline_ids=baseline_ids, baseline_report=baseline_report,
               geometry_audit=geometry_audit)
    return inp, chosen, result


def run(campaign, parent, tag, per_foil=None, check=False, shoulder_percent=20):
    if type(shoulder_percent) is not int or shoulder_percent not in (10, 20):
        raise ValueError("Shoulder percentage must be 10 or 20")
    mixture_name = f"inner{shoulder_percent}"
    schema = f"{mixture_name}_allocation_v1"
    if tag in ("inner10", "inner20") and tag != mixture_name:
        raise ValueError(f"Tag {tag} conflicts with --shoulder-percent {shoulder_percent}")
    campaign = campaign.expanduser().resolve()
    target = campaign / "05c_core_sample" / tag
    if not check and target.exists():
        raise FileExistsError(f"Output exists: {target}; choose a new tag")
    inp, chosen, result = prepare(campaign, parent, per_foil, shoulder_percent)
    print(f"{campaign.name}: {len(inp['required'])} physical foils; {result['foil_budget']} events per foil", flush=True)
    for row in result["foils"]:
        print(f"  foil {row['zfoil']:+g}: {row['core']} core + {row['inner_shoulder']} inner = {row['total']}", flush=True)
    print(f"Total {len(chosen)}; unit event weights; baseline cells covered; reserved IDs unchanged.", flush=True)
    if check:
        print("Check passed. No files written.", flush=True)
        return result
    target.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{mixture_name}_", dir=target.parent))
    try:
        for folder in ("root", "tsv", "metadata", "parent", "code"):
            (stage / folder).mkdir()
        shutil.copy2(inp["table"], stage / inp["table"].name)
        shutil.copy2(inp["metadata"], stage / "metadata/optics.dat")
        shutil.copy2(inp["maskpath"], stage / "metadata/sieve_mask.json")
        shutil.copy2(inp["source"] / "manifest.json", stage / "parent/manifest.json")
        shutil.copy2(inp["source"] / "config.json", stage / "parent/config.json")
        for directory in ("models", "code"):
            for path in sorted((inp["source"] / directory).rglob("*")):
                if path.is_file():
                    relative = str(path.relative_to(inp["source"]))
                    checked(path, inp["manifest"]["outputs"][relative])
            shutil.copytree(inp["source"] / directory, stage / ("models" if directory == "models" else "parent/code"))
        shutil.copy2(Path(__file__), stage / "code/allocate_inner20.py")
        for name in ("core_balance.py", "core_sample.py", "reallocate_core.py", "spectrometer_config.py", "spectrometer_profiles.def"):
            shutil.copy2(PROJECT / name, stage / "code" / name)
        write_tsv(stage / "tsv/selected_ids.tsv", chosen)
        write_tsv(stage / "tsv/baseline_ids.tsv", inp["baseline_ids"])
        write_tsv(stage / "tsv/geometry_exclusions.tsv", inp["geometry_audit"])
        write_tsv(stage / "tsv/foils.tsv", result["foils"])
        write_tsv(stage / "tsv/leaves.tsv", result["leaves"])
        for name, data in inp["data"].items():
            path = stage / "root" / f"CoreSample_{name}.root"
            with uproot.recreate(path) as root:
                root.mktree("CoreSample", data)
            with uproot.open(path) as root:
                back = root["CoreSample"].arrays(library="np")
            for field in data:
                if not np.array_equal(data[field], back[field], equal_nan=True):
                    raise ValueError(f"ROOT round-trip mismatch: {name}/{field}")
        inner_parts = shoulder_percent // 10
        policy = dict(name=mixture_name, core_fraction=(10-inner_parts)/10,
                      inner_shoulder_fraction=inner_parts/10, event_weight=1,
                      inner_quality=1, q_low_inclusive=.4, q_high_exclusive=.5,
                      core_pool="Reconstructed equal15 baseline core IDs",
                      allocation="Capped baseline-proportional counts with minimum-one coverage when affordable",
                      budget="Campaign-specific maximum count-feasible equal-foil budget" if per_foil is None else "Explicit per-foil budget")
        (stage / "config.json").write_text(json.dumps(inp["cfg"], indent=2) + "\n")
        (stage / "allocation.json").write_text(json.dumps(result, indent=2) + "\n")
        (stage / "baseline_equal15.json").write_text(json.dumps(inp["baseline_report"], indent=2) + "\n")
        (stage / "README.md").write_text(f"# {mixture_name} training allocation\n\n"
            f"Unit-weight {100-shoulder_percent}% cores and {shoulder_percent}% frozen quality-1 inner shoulders, 0.40 <= q < 0.50.\n"
            f"The quality/core_keep labels are unchanged; {mixture_name}_keep identifies training membership.\n"
            "All reserved event identities are preserved. See tsv/foils.tsv and selected_ids.tsv.\n"
            f"This {schema} tag requires the mixed TFit export adapter; do not use the legacy core-only builder.\n")
        for path, expected in inp["hashes"].items():
            checked(Path(path), expected)
        manifest = dict(schema=schema, mode="event_allocation", status="quotas_ready", parent=parent,
                        parent_manifest=digest(inp["source"] / "manifest.json"), config=inp["cfg"],
                        policy=policy, totals=dict(fit=len(chosen), core=(10-inner_parts) * len(chosen) // 10,
                                                 inner_shoulder=inner_parts * len(chosen) // 10),
                        foil_budget=result["foil_budget"], checks=dict(baseline_cells_preserved=True,
                        reserved_identities_preserved=True, quality_labels_unchanged=True, root_round_trip=True),
                        inputs=inp["hashes"], versions={name: importlib.metadata.version(name) for name in ("numpy", "uproot")})
        manifest["outputs"] = {str(path.relative_to(stage)): digest(path) for path in sorted(stage.rglob("*")) if path.is_file()}
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        stage.rename(target)
        print(f"Mixed allocation complete: {target}", flush=True)
        return result
    except BaseException:
        shutil.rmtree(stage)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("parent")
    parser.add_argument("tag", nargs="?", default="inner20")
    parser.add_argument("--per-foil", type=int, help="Override the campaign-specific count ceiling; multiple of ten")
    parser.add_argument("--shoulder-percent", type=int, choices=(10, 20), default=20,
                        help="Inner-shoulder event percentage per foil (default: 20)")
    parser.add_argument("--check", action="store_true", help="Compute and validate exact quotas/IDs without writing")
    args = parser.parse_args()
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", value) for value in (args.parent, args.tag)):
        parser.error("Invalid tag")
    try:
        run(args.campaign, args.parent, args.tag, args.per_foil, args.check, args.shoulder_percent)
    except (ValueError, OSError, KeyError, AssertionError) as error:
        parser.exit(1, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
