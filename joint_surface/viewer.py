"""Inspect saved memberships without fitting or changing any event assignments.

Run ``python3 -m joint_surface.viewer CAMPAIGN`` for the native desktop view,
or add ``--html membership.html`` for a self-contained, offline Firefox view.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


COORDINATES = ("xfp", "xpfp", "yfp", "ypfp", "delta")
PALETTE = ("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b",
           "#17becf", "#bcbd22", "#e377c2", "#393b79", "#637939", "#843c39")
REJECT_PALETTE = ("#d62728", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")


@dataclass
class Membership:
    x: np.ndarray
    y: np.ndarray
    accepted_x: np.ndarray
    accepted_y: np.ndarray
    reason_x: np.ndarray
    reason_y: np.ndarray
    extras: dict[str, np.ndarray]


@dataclass
class Case:
    key: str
    label: str
    entry: np.ndarray
    coordinates: np.ndarray  # cm, mrad, cm, mrad, percent, for display only
    configurations: dict[str, Membership]
    delta_edges: np.ndarray
    metadata: dict[str, Any]
    shown: np.ndarray


def _read_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def _array(data, names, n=None, default=None, dtype=None):
    for name in names:
        if name in data:
            result = np.asarray(data[name], dtype=dtype)
            break
    else:
        if default is None:
            raise ValueError("Missing array: " + " or ".join(names))
        result = np.full(n, default, dtype=dtype)
    if result.ndim != 1 or (n is not None and len(result) != n):
        raise ValueError(f"{names[0]} must be a one-dimensional array of length {n}")
    return result


def find_manifest(campaign: Path) -> Path:
    """Accept the campaign directory, work directory, or explicit manifest."""
    campaign = Path(campaign)
    if campaign.is_file():
        return campaign
    for relative in ("00_joint_surface_work/viewer_manifest.json", "viewer_manifest.json"):
        candidate = campaign / relative
        if candidate.is_file():
            return candidate
    raise ValueError(f"No viewer_manifest.json found in {campaign}; run the campaign label stage first")


def sample_indices(case: Case, max_points: int, seed: int = 731) -> np.ndarray:
    """A deterministic, bounded display sample with representation of rejections.

    Each rejection category first receives up to one quarter of the budget,
    shared between categories. Remaining places are uniformly sampled from all
    remaining events. This deliberate sampling is not a density estimator.
    """
    n = len(case.entry)
    available = np.flatnonzero(np.all(np.isfinite(case.coordinates), axis=1))
    if max_points <= 0 or len(available) <= max_points:
        return available
    salt = int.from_bytes(hashlib.sha256(case.key.encode()).digest()[:8], "little")
    rng = np.random.default_rng(seed ^ salt)
    groups = {}
    for config, membership in case.configurations.items():
        for axis in ("x", "y"):
            accepted = getattr(membership, "accepted_" + axis)
            reasons = getattr(membership, "reason_" + axis)
            for reason in np.unique(reasons[~accepted]):
                groups[(config, axis, str(reason))] = available[(~accepted & (reasons == reason))[available]]
    keep = set()
    budget = max_points // 4
    each = max(1, budget // max(1, len(groups)))
    for indices in groups.values():
        take = min(len(indices), each, budget - len(keep))
        if take > 0:
            keep.update(rng.choice(indices, take, replace=False).tolist())
    remaining = np.setdiff1d(available, np.fromiter(keep, dtype=int), assume_unique=False)
    need = max_points - len(keep)
    keep.update(rng.choice(remaining, need, replace=False).tolist())
    return np.array(sorted(keep), dtype=int)


def load_campaign(campaign: Path, max_points: int = 12000) -> list[Case]:
    manifest_path = find_manifest(campaign).resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != 1 or not manifest.get("cases"):
        raise ValueError("Expected a schema_version 1 manifest with a nonempty cases list")
    base = manifest_path.parent
    cases, keys = [], set()
    for item in manifest["cases"]:
        key = str(item["key"])
        if key in keys:
            raise ValueError(f"Duplicate case key: {key}")
        keys.add(key)
        events_path = base / item.get("events", item.get("npz", ""))
        data = _read_npz(events_path)
        entry = _array(data, ("entry", "entry_id"), dtype=np.int64)
        n = len(entry)
        if len(np.unique(entry)) != n:
            raise ValueError(f"{key}: duplicate entry IDs")
        coords = np.column_stack([_array(data, (name,), n, dtype=float) for name in COORDINATES])
        coords[:, [1, 3]] *= 1000  # Stored slopes are radians; display slopes are mrad.
        configurations = {}
        predictions = item.get("predictions", {"main": item.get("events", item.get("npz"))})
        for name, relative in predictions.items():
            pred = _read_npz(base / relative)
            if "entry" in pred and not np.array_equal(entry, pred["entry"]):
                raise ValueError(f"{key}/{name}: prediction event order differs from event coordinates")
            x = _array(pred, ("x_label", "physical_x", "latent_x"), n, dtype=np.int32)
            y = _array(pred, ("y_label", "physical_y", "latent_y"), n, dtype=np.int32)
            ax = _array(pred, ("accepted_x",), n, default=False, dtype=bool) if "accepted_x" in pred else x >= 0
            ay = _array(pred, ("accepted_y",), n, default=False, dtype=bool) if "accepted_y" in pred else y >= 0
            if np.any(ax & (x < 0)) or np.any(ay & (y < 0)):
                raise ValueError(f"{key}/{name}: accepted events need nonnegative membership indices")
            rx = _array(pred, ("reason_x", "rejection_x"), n, default="unassigned", dtype=str)
            ry = _array(pred, ("reason_y", "rejection_y"), n, default="unassigned", dtype=str)
            extras = {k: v for k, v in pred.items() if k.startswith(("proposed_", "latent_", "train_", "reference_", "final_")) and v.ndim == 1 and len(v) == n}
            configurations[str(name)] = Membership(x, y, ax, ay, rx, ry, extras)
        edges = np.asarray(item.get("delta_edges", manifest.get("delta_edges", [-10, -8, -5, 0, 5, 10])), dtype=float)
        if len(edges) < 2 or np.any(~np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
            raise ValueError(f"{key}: delta_edges must be finite and strictly increasing")
        case = Case(key, str(item.get("label", key)), entry, coords, configurations, edges, item, np.arange(n))
        case.shown = sample_indices(case, max_points)
        cases.append(case)
    return cases


def interval_mask(case: Case, interval: int | None) -> np.ndarray:
    if interval is None:
        return np.ones(len(case.entry), dtype=bool)
    lo, hi = case.delta_edges[interval:interval + 2]
    return (case.coordinates[:, 4] >= lo) & (case.coordinates[:, 4] < hi)


def counts(case: Case, config: str, interval=None, shown=False) -> dict[str, int]:
    selected = interval_mask(case, interval)
    if shown:
        sample = np.zeros(len(case.entry), dtype=bool)
        sample[case.shown] = True
        selected &= sample
    m = case.configurations[config]
    return {"events": int(selected.sum()), "x": int((selected & m.accepted_x).sum()),
            "y": int((selected & m.accepted_y).sum()),
            "both": int((selected & m.accepted_x & m.accepted_y).sum())}


def trace_groups(case: Case, config: str, axis: str, color_mode="each", interval=None, event_filter="all"):
    """Return explicit identity/rejection groups shared by desktop and HTML."""
    m = case.configurations[config]
    dimension = axis if color_mode == "each" else color_mode
    accepted = getattr(m, "accepted_" + dimension)
    labels = getattr(m, dimension)
    reasons = getattr(m, "reason_" + dimension)
    mask = interval_mask(case, interval)
    if event_filter == "accepted":
        mask &= m.accepted_x & m.accepted_y
    elif event_filter == "rejected":
        mask &= ~(m.accepted_x & m.accepted_y)
    indices = case.shown[mask[case.shown]]
    result = []
    for identity in np.unique(labels[indices][accepted[indices]]):
        ii = indices[accepted[indices] & (labels[indices] == identity)]
        result.append((f"{dimension.upper()} {'row' if dimension == 'x' else 'column'} {identity}",
                       PALETTE[int(identity) % len(PALETTE)], "o", ii))
    for k, reason in enumerate(sorted(set(reasons[indices][~accepted[indices]]))):
        ii = indices[~accepted[indices] & (reasons[indices] == reason)]
        result.append((f"{dimension.upper()} rejected: {reason or 'unassigned'}", REJECT_PALETTE[k % len(REJECT_PALETTE)], "x", ii))
    return result


def event_text(case: Case, config: str, index: int) -> str:
    m = case.configurations[config]
    c = case.coordinates[index]
    def identity(axis):
        return str(getattr(m, axis)[index]) if getattr(m, "accepted_" + axis)[index] else "rejected (" + str(getattr(m, "reason_" + axis)[index]) + ")"
    parts = [f"Entry {case.entry[index]}", f"delta {c[4]:.4g}%", f"X {identity('x')}", f"Y {identity('y')}",
             f"xfp {c[0]:.4g} cm", f"xpfp {c[1]:.4g} mrad", f"yfp {c[2]:.4g} cm", f"ypfp {c[3]:.4g} mrad"]
    parts.extend(f"{key}={value[index]}" for key, value in sorted(m.extras.items()))
    return " · ".join(parts)


def html_payload(cases):
    result = []
    for c in cases:
        ii = c.shown
        result.append({"key": c.key, "label": c.label, "entry": c.entry[ii].tolist(),
                       "coordinates": c.coordinates[ii].tolist(), "delta_edges": c.delta_edges.tolist(),
                       "full": {name: [counts(c, name), *[counts(c, name, j) for j in range(len(c.delta_edges)-1)]] for name in c.configurations},
                       "configurations": {name: {"x": m.x[ii].tolist(), "y": m.y[ii].tolist(),
                           "accepted_x": m.accepted_x[ii].tolist(), "accepted_y": m.accepted_y[ii].tolist(),
                           "reason_x": m.reason_x[ii].tolist(), "reason_y": m.reason_y[ii].tolist(),
                           "extras": {k: v[ii].tolist() for k, v in m.extras.items()}}
                           for name, m in c.configurations.items()}})
    return {"schema_version": 1, "cases": result, "palette": PALETTE, "rejection_palette": REJECT_PALETTE}


def write_html(cases, output: Path, force=False, plotly_js: Path | None = None):
    """Embed Plotly's installed/local JavaScript; the resulting file is offline."""
    if output.exists() and not force:
        raise ValueError(f"Output exists: {output}; use --force to replace it")
    if plotly_js:
        js = Path(plotly_js).read_text()
    else:
        try:
            from plotly.offline.offline import get_plotlyjs
        except ImportError as exc:
            raise ValueError("HTML export needs the optional plotly Python package (pip install plotly), or --plotly-js /path/to/plotly.min.js. The exported HTML then needs no network.") from exc
        js = get_plotlyjs()
    template = Path(__file__).parents[1] / "tools/membership_viewer/viewer.html"
    data = json.dumps(html_payload(cases), separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    html = template.read_text().replace("__PLOTLY_JS__", js).replace("__MEMBERSHIP_DATA__", data)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w" if force else "x") as handle:
        handle.write(html)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path, help="Campaign folder or viewer_manifest.json")
    parser.add_argument("--desktop", action="store_true", help="Open native Matplotlib/Tk window (default)")
    parser.add_argument("--html", type=Path, help="Export a single offline HTML file instead of opening a window")
    parser.add_argument("--plotly-js", type=Path, help="Local Plotly bundle, as an alternative to the optional plotly package")
    parser.add_argument("--max-points", type=int, default=12000, help="Display budget per foil; 0 means all events (default: 12000)")
    parser.add_argument("--force", action="store_true", help="Replace an existing HTML output")
    args = parser.parse_args(argv)
    try:
        cases = load_campaign(args.campaign, args.max_points)
        if args.html:
            output = write_html(cases, args.html, args.force, args.plotly_js)
            print(f"Created {output}: self-contained; open locally with Firefox")
        else:
            from .desktop_viewer import launch
            launch(cases)
    except (ValueError, OSError, KeyError, ImportError) as exc:
        parser.exit(2, f"membership viewer: {exc}\n")
    return 0


if __name__ == "__main__":
    main()
