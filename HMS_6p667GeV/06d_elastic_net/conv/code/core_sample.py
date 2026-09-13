#!/usr/bin/env python3
"""Label-seeded sieve cores, protected holdout, and campaign fit allocation.

Run through run_core_sample.sh. See CORE_SAMPLE.md for the sample definitions.
The selector reads pre-GMM candidate trees; it never updates cuts or matrices.
"""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import re
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import uproot
from scipy import ndimage as ndi
from skimage.segmentation import watershed


PROJECT = Path(__file__).resolve().parent
SAMPLES = {0: "noncore", 1: "fit", 2: "holdout", 3: "surplus"}


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0]) if rows else ["entry", "reason"]
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def ranks(seed, tag, entries, purpose):
    # Separate random streams: changing a fit budget cannot change the holdout.
    prefix = f"{seed}:{purpose}:{tag}:"
    return np.array([int.from_bytes(hashlib.blake2b(
        (prefix + str(int(i))).encode(), digest_size=8).digest(), "big")
        for i in entries], dtype=np.uint64)


def groups(data, columns, indices=None):
    out = defaultdict(list)
    if indices is None:
        indices = range(len(data["entry"]))
    for i in indices:
        out[tuple(data[k][i].item() for k in columns)].append(i)
    return {k: np.asarray(v, dtype=int) for k, v in sorted(out.items())}


def load_config(campaign):
    cfg = json.loads((PROJECT / "core_sample.json").read_text())
    override = campaign / "config/core_sample.json"
    if override.exists():
        changes = json.loads(override.read_text())
        unknown = changes.keys() - cfg.keys()
        if unknown:
            raise ValueError(f"Unknown core settings: {sorted(unknown)}")
        cfg.update(changes)
    for name in ("holdout", "core", "fit_fraction"):
        if not 0 < cfg[name] < 1:
            raise ValueError(f"{name} must be between 0 and 1")
    if not .03 <= cfg["smooth"] <= .25:
        raise ValueError("smooth must be between 0.03 and 0.25 hole spacings")
    for name in ("min_events", "fit_cap", "fit_max"):
        if type(cfg[name]) is not int or cfg[name] < 2:
            raise ValueError(f"{name} must be an integer >= 2")
    if type(cfg["seed"]) is not int:
        raise ValueError("seed must be an integer")
    if cfg["peak_ratio"] <= 1 or not 0 < cfg["max_shift"] < .5:
        raise ValueError("Require peak_ratio > 1 and 0 < max_shift < 0.5")
    if cfg["spacing"] is not None:
        if len(cfg["spacing"]) != 2 or not all(
                np.isfinite(x) and x > 0 for x in cfg["spacing"]):
            raise ValueError("spacing must be null or [positive X cm, positive Y cm]")
    if not isinstance(cfg["pools"], dict) or not all(
            isinstance(v, str) and v for v in cfg["pools"].values()):
        raise ValueError("pools must map rungroup names to nonempty pool names")
    return cfg


def read_joint(xpath, ypath, optics_id):
    common = ["run", "entry", "foil", "ndel", "zfoil", "delta_low", "delta_high",
              "xsieve", "ysieve", "xfp", "xpfp", "yfp", "ypfp", "delta"]
    arrays = []
    for path, tree, col in ((xpath, "TXCand", "xscol"), (ypath, "TYCand", "yscol")):
        with uproot.open(path) as root:
            if tree not in root:
                raise ValueError(f"Missing {tree}: {path}")
            missing = set(common + [col]) - set(root[tree].keys())
            if missing:
                raise ValueError(f"{path}: missing branches {sorted(missing)}")
            arrays.append(root[tree].arrays(library="np"))
    x, y = arrays
    # The current producer stores the optics ID and the entry in the merged T.
    for a in arrays:
        if np.any(a["run"] != int(optics_id)) or np.any(a["entry"] < 0):
            raise ValueError("Candidate run/entry does not match campaign optics ID")
    ux, ix, nx = np.unique(x["entry"], return_index=True, return_counts=True)
    uy, iy, ny = np.unique(y["entry"], return_index=True, return_counts=True)
    joint, jx, jy = np.intersect1d(ux, uy, return_indices=True)
    valid = (nx[jx] == 1) & (ny[jy] == 1)
    audit = [{"entry": int(e), "reason": "missing_axis"}
             for e in np.setxor1d(ux, uy)]
    audit += [{"entry": int(e), "reason": "multiple_labels"} for e in joint[~valid]]
    xi, yi = ix[jx[valid]], iy[jy[valid]]
    # Refuse to join candidate trees made from different replays or slice cuts.
    for key in common:
        if not np.allclose(x[key][xi], y[key][yi], rtol=0, atol=1e-9, equal_nan=True):
            raise ValueError(f"X/Y candidates disagree in {key}; rebuild both stages")
    data = {k: v[yi].copy() for k, v in y.items() if v.dtype.kind in "iufb"}
    data["xscol"] = x["xscol"][xi].astype(np.int32)
    usable = np.ones(len(xi), dtype=bool)
    for key in common:
        usable &= np.isfinite(data[key])
    usable &= (data["xscol"] >= 0) & (data["yscol"] >= 0)
    usable &= (data["foil"] >= 0) & (data["ndel"] >= 0)
    usable &= data["delta_low"] < data["delta_high"]
    usable &= (data["delta"] >= data["delta_low"]) & (data["delta"] < data["delta_high"])
    audit += [{"entry": int(e), "reason": "invalid_values"} for e in data["entry"][~usable]]
    data = {k: v[usable] for k, v in data.items()}
    if not len(data["entry"]):
        raise ValueError("No unambiguous joint candidates")
    return data, sorted(audit, key=lambda row: row["entry"])


def reserve(data, cfg, tag):
    reserved = np.zeros(len(data["entry"]), dtype=bool)
    order = ranks(cfg["seed"], tag, data["entry"], "holdout")
    for idx in groups(data, ["run", "foil", "ndel", "xscol", "yscol"]).values():
        # A singleton goes to holdout and cannot establish a core by itself.
        n = min(len(idx), max(1, math.ceil(cfg["holdout"] * len(idx))))
        reserved[idx[np.argsort(order[idx], kind="stable")[:n]]] = True
    return reserved


def spacing_from_labels(data, dev):
    spacing = []
    for label, coordinate in (("xscol", "xsieve"), ("yscol", "ysieve")):
        medians = [(k[0], np.median(data[coordinate][idx]))
                   for k, idx in groups(data, [label], dev).items() if len(idx) >= 5]
        slopes = [abs((b[1] - a[1]) / (b[0] - a[0]))
                  for a, b in zip(medians, medians[1:]) if b[0] != a[0]]
        slopes = [v for v in slopes if np.isfinite(v) and v > .01]
        if not slopes:
            raise ValueError("Cannot infer both hole spacings; set spacing: [X cm, Y cm] in config/core_sample.json")
        spacing.append(float(np.median(slopes)))
    return np.array(spacing)


def connected(mask, peak):
    labels, _ = ndi.label(mask)
    number = labels[peak]
    return labels == number if number else np.zeros_like(mask)


def select_slice(data, idx, cfg, reserved, spacing):
    """Fit a frozen map on development events and score both partitions."""
    points = np.column_stack((data["xsieve"][idx], data["ysieve"][idx])) / spacing
    train = ~reserved[idx]
    holes = groups(data, ["xscol", "yscol"], idx)
    # Axis bounds and seeds use development data only.
    if not np.any(train):
        raise ValueError("Slice has no development events")
    lo = np.quantile(points[train], .001, axis=0) - .75
    hi = np.quantile(points[train], .999, axis=0) + .75
    lookup = {int(i): j for j, i in enumerate(idx)}
    # Global quantiles alone can clip an outer hole with a tiny occupancy share.
    # Include every sufficiently supported label's center in the map bounds.
    for members in holes.values():
        local = np.array([lookup[int(i)] for i in members])
        own = local[train[local]]
        if len(own) >= cfg["min_events"]:
            center = np.median(points[own], axis=0)
            lo = np.minimum(lo, center - .75)
            hi = np.maximum(hi, center + .75)
    bins = np.maximum(40, np.ceil((hi - lo) * 30).astype(int))
    if np.prod(bins) > 4_000_000:
        raise ValueError("Unreasonably large sieve grid; inspect coordinates/spacing")
    edges = [np.linspace(lo[d], hi[d], bins[d] + 1) for d in range(2)]
    axes = [(e[1:] + e[:-1]) / 2 for e in edges]
    gx, gy = np.meshgrid(*axes, indexing="ij")
    step = (hi - lo) / bins
    hist = np.histogram2d(*points[train].T, bins=edges)[0]
    density = ndi.gaussian_filter(hist, cfg["smooth"] / step, mode="constant")
    markers = np.zeros_like(hist, dtype=np.int32)
    pixels = np.floor((points - lo) / step).astype(int)
    inside = np.all((pixels >= 0) & (pixels < bins), axis=1)
    pixels = np.clip(pixels, 0, bins - 1)
    info = []
    for number, (hole, members) in enumerate(holes.items(), start=1):
        local = np.array([lookup[int(i)] for i in members])
        own = local[train[local]]
        row = dict(xscol=hole[0], yscol=hole[1], n_dev=len(own), status="low_stats")
        item = dict(number=number, local=local, own=own, row=row, peak=None)
        info.append(item)
        if len(own) < cfg["min_events"]:
            continue
        center = np.median(points[own], axis=0)
        search = ((gx - center[0]) ** 2 + (gy - center[1]) ** 2) < .4 ** 2
        # Use this label's map to seed the global map, preventing a busy neighbor
        # from supplying the seed of a sparse label.
        own_hist = np.histogram2d(*points[own].T, bins=edges)[0]
        own_density = ndi.gaussian_filter(own_hist, cfg["smooth"] / step, mode="constant")
        peak = np.unravel_index(np.argmax(np.where(search, own_density, -1)), bins)
        if markers[peak]:
            row["status"] = "seed_collision"
            continue
        markers[peak] = number
        item["peak"] = peak
    # A watershed assigns the intervening valleys. Core masks are additionally
    # bounded to 0.65 spacing from their marker and to their own label's events.
    basin = watershed(-density, markers) if markers.any() else markers.copy()
    coremap = np.zeros_like(markers)
    quality = np.zeros(len(idx), dtype=np.int32)  # 0 unsupported, 1 noncore, 2 core
    score = np.full(len(idx), np.nan)
    rows = []
    for item in info:
        number, local, own, row, peak = (item[k] for k in ("number", "local", "own", "row", "peak"))
        row.update(peak_x=float("nan"), peak_y=float("nan"), background=float("nan"),
                   peak_density=float("nan"), shift=float("nan"), center_x=float("nan"),
                   center_y=float("nan"), n_core_dev=0)
        if peak is not None:
            px, py = axes[0][peak[0]], axes[1][peak[1]]
            radius = np.hypot(gx - px, gy - py)
            annulus = (radius >= .45) & (radius <= .7)
            bg = float(np.median(density[annulus])) if annulus.any() else 0.
            height = float(density[peak])
            row.update(peak_x=px * spacing[0], peak_y=py * spacing[1],
                       background=bg, peak_density=height)
            row["status"] = "weak_peak"
            if height > max(1e-12, cfg["peak_ratio"] * bg):
                q = (density - bg) / (height - bg)
                region = (basin == number) & (radius <= .65)
                mask = connected(region & (q >= cfg["core"]), peak)
                tiers = [connected(region & (q >= t), peak)
                         for t in (max(.05, cfg["core"] - .15), min(.95, cfg["core"] + .15))]
                centers = []
                for tier in tiers:
                    selected = own[inside[own] & tier[tuple(pixels[own].T)]]
                    if len(selected) >= 5:
                        centers.append(np.mean(points[selected], axis=0))
                selected = own[inside[own] & mask[tuple(pixels[own].T)]]
                row["n_core_dev"] = len(selected)
                row["status"] = "unstable"
                if len(centers) == 2:
                    row["shift"] = float(np.linalg.norm(centers[1] - centers[0]))
                if len(selected) >= 5 and len(centers) == 2 and row["shift"] <= cfg["max_shift"]:
                    center = np.mean(points[selected], axis=0) * spacing
                    row.update(status="accepted", center_x=center[0], center_y=center[1])
                    quality[local] = 1
                    chosen = local[inside[local] & mask[tuple(pixels[local].T)]]
                    quality[chosen] = 2
                    score[local[inside[local]]] = q[tuple(pixels[local[inside[local]]].T)]
                    coremap[mask] = number
        rows.append(row)
    model = dict(xedges=edges[0] * spacing[0], yedges=edges[1] * spacing[1],
                 density=density, basin=basin, coremap=coremap, spacing=spacing,
                 hole_ids=np.array(list(holes)), quality=quality, score=score)
    return model, rows


def fair_counts(capacities, budget):
    """Equal shares, redistributing an undersupplied group's unused share."""
    out = np.zeros(len(capacities), dtype=int)
    remaining = min(budget, sum(capacities))
    while remaining:
        active = np.flatnonzero(out < capacities)
        share = max(1, remaining // len(active))
        for i in active:
            take = min(share, capacities[i] - out[i], remaining)
            out[i] += take
            remaining -= take
            if not remaining:
                break
    return out


def allocate(all_data, cfg):
    pools = defaultdict(list)
    for tag, data in sorted(all_data.items()):
        eligible = np.flatnonzero((data["quality"] == 2) & (data["sample"] != 2))
        data["sample"][eligible] = 3
        for key, idx in groups(data, ["zfoil", "delta_low", "delta_high", "xscol", "yscol"], eligible).items():
            # Without an explicit mapping, different settings never share a cap.
            pools[(cfg["pools"].get(tag, tag), *key)].append((tag, idx))
    allocation = []
    ordered = sorted(pools.items())
    pool_caps = [min(cfg["fit_cap"], sum(int(math.floor(cfg["fit_fraction"] * len(idx)))
                                       for _, idx in parts)) for _, parts in ordered]
    if cfg["fit_max"] < sum(cap > 0 for cap in pool_caps):
        raise ValueError("fit_max is too small to represent every populated allocation cell")
    budgets = fair_counts(np.array(pool_caps), cfg["fit_max"])
    for (key, parts), budget in zip(ordered, budgets):
        capacities = np.array([int(math.floor(cfg["fit_fraction"] * len(idx))) for _, idx in parts])
        counts = fair_counts(capacities, budget)
        for (tag, idx), nfit in zip(parts, counts):
            data = all_data[tag]
            order = np.argsort(ranks(cfg["seed"], tag, data["entry"][idx], "fit"), kind="stable")
            data["sample"][idx[order[:nfit]]] = 1
            allocation.append(dict(pool=key[0], rungroup=tag, zfoil=key[1], delta_low=key[2],
                                   delta_high=key[3], xscol=key[4], yscol=key[5],
                                   available=len(idx), fit=int(nfit), surplus=len(idx) - int(nfit)))
    return allocation


def plot_slice(path, data, idx, model, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    x, y = data["xsieve"][idx], data["ysieve"][idx]
    core, sample = data["quality"][idx] == 2, data["sample"][idx]
    masks = [np.ones(len(idx), bool), sample == 1, (sample == 2) & core,
             sample == 3, (sample != 2) & ~core, (sample == 2) & ~core]
    titles = ["Joint labels + core boundaries", "Fit cores", "Protected holdout: cores",
              "Surplus cores", "Development: noncore / unsupported", "Protected holdout: noncore / unsupported"]
    hist = np.histogram2d(x, y, bins=[model["xedges"], model["yedges"]])[0]
    norm = LogNorm(vmin=1, vmax=max(2, hist.max()))
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)
    xc = (model["xedges"][1:] + model["xedges"][:-1]) / 2
    yc = (model["yedges"][1:] + model["yedges"][:-1]) / 2
    for ax, mask, name in zip(axes.flat, masks, titles):
        # Match Hall C convention: horizontal Y sieve, vertical X sieve.
        h = np.histogram2d(x[mask], y[mask], bins=[model["xedges"], model["yedges"]])[0]
        artist = ax.pcolormesh(model["yedges"], model["xedges"], np.ma.masked_less(h, 1),
                               norm=norm, cmap="viridis", shading="auto")
        clipped = int(mask.sum() - h.sum())
        ax.set_title(f"{name}\nN={mask.sum():,}" + (f"; outside frame={clipped}" if clipped else ""), fontsize=10)
        ax.set_xlabel("Y sieve (cm)")
        ax.set_ylabel("X sieve (cm)")
        ax.set_aspect("equal")
    for number, hole in enumerate(model["hole_ids"], 1):
        region = model["coremap"] == number
        if region.any():
            axes.flat[0].contour(yc, xc, region, levels=[.5], colors=["crimson"], linewidths=.7)
            pos = np.argwhere(region).mean(axis=0).astype(int)
            axes.flat[0].text(yc[pos[1]], xc[pos[0]], f"{hole[0]},{hole[1]}", fontsize=6, color="black")
    fig.suptitle(title + " | labels shown as X,Y", fontsize=13)
    fig.colorbar(artist, ax=list(axes.flat), label="Events / bin (shared log scale)", fraction=.02)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(campaign, tag, check=False):
    cfg = load_config(campaign)
    tables = sorted((campaign / "config").glob("rungroups_*_inputs.tsv"))
    if len(tables) != 1:
        raise ValueError(f"Expected one rungroup table in {campaign / 'config'}; "
                         f"found {len(tables)}: {[p.name for p in tables]}. "
                         "Check the campaign path and branch configuration.")
    settings = read_tsv(tables[0])
    names = [row["rungroup"] for row in settings]
    if not settings or len(set(names)) != len(names):
        raise ValueError("Empty campaign or duplicate rungroup names")
    if len({row["optics_id"] for row in settings}) != len(settings):
        raise ValueError("Duplicate optics IDs in campaign table")
    # The available entry key is local to a merged replay. Shared original runs
    # could otherwise put the same physical event in fit and holdout files.
    seen_runs = set()
    seen_sources = set()
    for row in settings:
        original_runs = {value.strip() for value in row["runs"].split(",")}
        if not all(value.isdigit() for value in original_runs):
            raise ValueError("Campaign runs must be comma-separated original run numbers")
        if seen_runs & original_runs or row["rootfile"] in seen_sources:
            raise ValueError("Rungroups must have disjoint original runs and source files")
        seen_runs.update(original_runs)
        seen_sources.add(row["rootfile"])
    if set(cfg["pools"]) - set(names):
        raise ValueError("pools contains rungroups absent from the campaign table")
    sources = []
    for row in settings:
        name = row["rungroup"]
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
            raise ValueError(f"Invalid rungroup name: {name}")
        paths = [campaign / f"04{axis}_candidate_trees_{letter}/root/{col}Candidates_{name}.root"
                 for axis, letter, col in (("b", "x", "Xscol"), ("a", "y", "Yscol"))]
        sources.append((row, paths))
    missing = [str(p) for _, paths in sources for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing candidate files; run both candidate stages first:\n" + "\n".join(missing))
    print(f"Campaign: {campaign.name}; {len(settings)} rungroups; settings: {json.dumps(cfg)}", flush=True)
    if check:
        for row, paths in sources:
            data, audit = read_joint(*paths, row["optics_id"])
            print(f"OK {row['rungroup']}: {len(data['entry'])} joint events; {len(audit)} excluded", flush=True)
        print("Preflight complete; no outputs written.")
        return
    target = campaign / "05c_core_sample" / tag
    if target.exists():
        raise FileExistsError(f"Output exists: {target}. Use a new short TAG to preserve the earlier split.")
    target.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".core_", dir=target.parent))
    try:
        for sub in ("root", "plots", "tsv", "models", "code"):
            (stage / sub).mkdir()
        (stage / "config.json").write_text(json.dumps(cfg, indent=2) + "\n")
        (stage / "code/core_sample.json").write_text(json.dumps(cfg, indent=2) + "\n")
        shutil.copy2(PROJECT / "core_sample.py", stage / "code/core_sample.py")
        shutil.copy2(tables[0], stage / tables[0].name)
        files = [tables[0], PROJECT / "core_sample.py", PROJECT / "core_sample.json"]
        files += [p for _, paths in sources for p in paths]
        manifest = dict(config=cfg, inputs={str(p.resolve()): digest(p) for p in files},
                        versions={m: importlib.metadata.version(m) for m in
                                  ("numpy", "scipy", "uproot", "scikit-image", "matplotlib")},
                        python=sys.version, samples=SAMPLES,
                        quality={0: "unsupported", 1: "noncore", 2: "core"})
        (stage / "requirements.txt").write_text("".join(f"{m}=={v}\n" for m, v in manifest["versions"].items()))
        all_data, summaries, plots = {}, [], []
        for row, paths in sources:
            name = row["rungroup"]
            data, audit = read_joint(*paths, row["optics_id"])
            write_tsv(stage / "tsv" / f"core_{name}_excluded.tsv", audit)
            protected = reserve(data, cfg, name)
            data["sample"] = np.where(protected, 2, 0).astype(np.int32)
            data["quality"] = np.zeros(len(protected), dtype=np.int32)
            data["core_score"] = np.full(len(protected), np.nan)
            spacing = np.array(cfg["spacing"]) if cfg["spacing"] else spacing_from_labels(data, np.flatnonzero(~protected))
            print(f"Select {name}: N={len(protected)}, spacing={spacing.round(4)} cm", flush=True)
            for (foil, ndel), idx in groups(data, ["foil", "ndel"]).items():
                stem = f"core_{name}_foil{foil}_ndel{ndel}"
                if not np.any(~protected[idx]):
                    summaries.append(dict(rungroup=name, foil=foil, ndel=ndel, status="no_development"))
                    continue
                model, rows = select_slice(data, idx, cfg, protected, spacing)
                data["quality"][idx] = model.pop("quality")
                data["core_score"][idx] = model.pop("score")
                np.savez_compressed(stage / "models" / (stem + ".npz"), **model)
                for summary in rows:
                    summary.update(rungroup=name, foil=foil, ndel=ndel)
                    summaries.append(summary)
                plots.append((name, idx, stem, foil, ndel))
            all_data[name] = data
        allocation = allocate(all_data, cfg)
        write_tsv(stage / "tsv/allocation.tsv", allocation,
                  ["pool", "rungroup", "zfoil", "delta_low", "delta_high", "xscol", "yscol", "available", "fit", "surplus"])
        counts = []
        for name, data in all_data.items():
            data["core_keep"] = (data["quality"] == 2).astype(np.int32)
            with uproot.recreate(stage / "root" / f"CoreSample_{name}.root") as root:
                root.mktree("CoreSample", data)
            for key, idx in groups(data, ["foil", "ndel", "xscol", "yscol"]).items():
                counts.append(dict(rungroup=name, foil=key[0], ndel=key[1], xscol=key[2], yscol=key[3],
                                   total=len(idx), fit=int(np.sum(data["sample"][idx] == 1)),
                                   holdout=int(np.sum(data["sample"][idx] == 2)),
                                   holdout_core=int(np.sum((data["sample"][idx] == 2) & (data["quality"][idx] == 2))),
                                   surplus=int(np.sum(data["sample"][idx] == 3)),
                                   unsupported=int(np.sum(data["quality"][idx] == 0))))
        write_tsv(stage / "tsv/counts.tsv", counts)
        fields = ["rungroup", "foil", "ndel", "xscol", "yscol", "status", "n_dev", "n_core_dev",
                  "peak_x", "peak_y", "center_x", "center_y", "peak_density", "background", "shift"]
        write_tsv(stage / "tsv/regions.tsv", summaries, fields)
        for name, idx, stem, foil, ndel in plots:
            with np.load(stage / "models" / (stem + ".npz")) as model:
                plot_slice(stage / "plots" / (stem + ".png"), all_data[name], idx, model,
                           f"{name} | foil {foil}, delta slice {ndel}")
        manifest["totals"] = {key: sum(r[key] for r in counts) for key in
                              ("total", "fit", "holdout", "holdout_core", "surplus", "unsupported")}
        manifest["outputs"] = {str(p.relative_to(stage)): digest(p) for p in sorted(stage.rglob("*")) if p.is_file()}
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        stage.rename(target)
        print(f"Core sample complete: {target}\n{json.dumps(manifest['totals'])}")
        if manifest["totals"]["fit"] == 0:
            print("WARNING: no fit cores accepted; inspect regions.tsv and the plots.")
    except BaseException:
        shutil.rmtree(stage)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("tag", nargs="?", default="core")
    parser.add_argument("--check", action="store_true", help="validate inputs without writing outputs")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.tag):
        parser.error("TAG must be a short name containing letters, digits, underscores or hyphens")
    campaign = args.campaign if args.campaign.is_absolute() else PROJECT / args.campaign
    try:
        run(campaign.resolve(), args.tag, args.check)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
