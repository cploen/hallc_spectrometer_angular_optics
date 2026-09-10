#!/usr/bin/env python3
"""Compare available, training and capped core counts at nominal sieve positions."""
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Circle
import numpy as np

from spectrometer_config import from_campaign


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", type=Path)
    parser.add_argument("rungroup", help="Name or unique prefix")
    parser.add_argument("--delta", type=int, nargs="+", default=[0, 1])
    parser.add_argument("--cap", type=float, default=2, help="Multiplier of the 25th percentile")
    parser.add_argument("--log", action="store_true", help="Use logarithmic colors; save a separate _log PNG")
    args = parser.parse_args()
    if not math.isfinite(args.cap) or args.cap <= 0:
        parser.error("cap must be positive")
    campaign = args.sample.parent.parent
    spec = from_campaign(campaign)
    mask = campaign / "config/sieve_mask.json"
    blocked = {tuple(p) for p in json.loads(mask.read_text())["blocked"]} if mask.exists() else set()
    with (args.sample / "tsv/counts.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    names = {r["rungroup"] for r in rows if r["rungroup"].startswith(args.rungroup)}
    if len(names) != 1:
        parser.error("Rungroup must match exactly one saved name")
    name = names.pop()
    rows = [r for r in rows if r["rungroup"] == name and int(r["ndel"]) in args.delta]
    if set(args.delta) - {int(r["ndel"]) for r in rows}:
        parser.error("A requested delta slice is absent")
    vmax = max((int(r["fit"]) + int(r["surplus"]) for r in rows
                if (int(r["xscol"]), int(r["yscol"])) not in blocked), default=1)
    norm = LogNorm(1, max(vmax, 2)) if args.log else Normalize(0, max(vmax, 1))
    cmap = plt.get_cmap("viridis")
    plt.rcParams.update({"font.size": 11})
    for foil, delta in sorted({(int(r["foil"]), int(r["ndel"])) for r in rows}):
        cells = {(int(r["xscol"]), int(r["yscol"])): r for r in rows
                 if int(r["foil"]) == foil and int(r["ndel"]) == delta}
        available = {k: int(r["fit"]) + int(r["surplus"]) for k, r in cells.items() if k not in blocked}
        training = {k: int(r["fit"]) for k, r in cells.items() if k not in blocked}
        positive = [n for n in available.values() if n > 0]
        if not positive:
            parser.error(f"No available cores for foil {foil}, delta {delta}")
        p25 = float(np.percentile(positive, 25))
        cap = math.floor(args.cap * p25)
        proposed = {k: min(n, cap) for k, n in available.items()}
        assert all(training[k] <= available[k] and proposed[k] <= available[k] for k in available)
        fig, axes = plt.subplots(1, 3, figsize=(17, 10))
        fig.subplots_adjust(left=.045, right=.91, bottom=.27, top=.83, wspace=.23)
        for ax, title, values in zip(axes, ["Available development core events", "Current training allocation",
                                          f"Proposed: {args.cap:g} × P25 cap = {cap}"],
                                     [available, training, proposed]):
            for ix in range(spec.nx):
                for iy in range(spec.ny):
                    pair = (ix, iy)
                    x, y = spec.ys(iy), spec.xs(ix)
                    if pair in blocked:
                        ax.text(x, y, "×", ha="center", va="center", fontsize=20, color="#b22525")
                    elif pair not in values:
                        ax.plot(x, y, ".", color="#b8bdc5", markersize=4)
                    else:
                        value = values[pair]
                        fill = cmap(norm(value)) if value else "#ffffff"
                        ax.add_patch(Circle((x, y), .62, facecolor=fill,
                                            edgecolor="#68717a", linewidth=.6))
                        rgb = cmap(norm(value))[:3]
                        light = .2126 * rgb[0] + .7152 * rgb[1] + .0722 * rgb[2]
                        ink = "white" if value and light < .5 else "#151515"
                        ax.text(x, y, str(value), ha="center", va="center", fontsize=9, color=ink)
            nonzero = [n for n in values.values() if n > 0]
            ratio = max(nonzero) / min(nonzero) if nonzero else float("nan")
            ax.set_title(title, fontsize=13, pad=16)
            ax.text(.5, -.095, f"Total: {sum(values.values()):,} events\n"
                    f"Largest / smallest nonzero: {ratio:.1f}×", transform=ax.transAxes,
                    ha="center", va="top", fontsize=12, linespacing=1.6)
            ax.set_xlabel("Y sieve (cm)")
            ax.set_ylabel("X sieve (cm)")
            ax.set_xlim(spec.ys(0) - 1, spec.ys(spec.ny - 1) + 1)
            ax.set_ylim(spec.xs(0) - 1.3, spec.xs(spec.nx - 1) + 1.3)
            ax.set_aspect("equal")
            ax.spines[["top", "right"]].set_visible(False)
        cax = fig.add_axes([.93, .28, .014, .44])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                     label=f"Core events per hole ({'log' if args.log else 'linear'} scale; shared across panels)")
        fig.suptitle(f"{campaign.name} · {name} · foil {foil} · delta slice {delta}\n"
                     f"Sieve population balance | {args.sample.name} | {'Log' if args.log else 'Linear'} colors", fontsize=18, y=.97)
        fig.text(.5, .065,
                 f"P25 = {p25:g}; proposed cap retains {sum(proposed.values()) / sum(available.values()):.0%} "
                 f"of available core events; {sum(n > cap for n in available.values())} holes capped. "
                 f"{sum(n == 0 for n in available.values())} additional labels have no eligible core events.\n"
                 "× = confirmed blocked position; hollow 0 = no eligible core events (outside color scale); gray dot = no saved label.\n"
                 "Nominal geometry, not reconstructed centroids. Protected holdout excluded. "
                 "Proposal replaces current allocation rules; no 80% rule or campaign budget applied.",
                 ha="center", va="center", fontsize=10, linespacing=1.6)
        if not mask.exists():
            fig.text(.5, .015, "No blocked-position mask supplied.", ha="center", color="#b22525")
        suffix = f"_cap{args.cap:g}".replace(".", "p") if args.cap != 2 else ""
        suffix += "_log" if args.log else ""
        out = args.sample / "plots" / f"balance_{name}_foil{foil}_ndel{delta}{suffix}.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=160, facecolor="white")
        plt.close(fig)
        print(f"{out}: available={sum(available.values())}, training={sum(training.values())}, "
              f"proposed={sum(proposed.values())}, cap={cap}")


if __name__ == "__main__":
    main()
