#!/usr/bin/env python3
"""Map historical fit inputs versus the log-verified SVD selection."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from spectrometer_config import from_campaign


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("campaign", type=Path)
    p.add_argument("data", type=Path, help="Verified svd_counts.py output directory")
    args = p.parse_args()
    spec = from_campaign(args.campaign)
    with (args.data / "counts.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    provenance = json.loads((args.data / "provenance.json").read_text())
    assert sum(int(r["training"]) for r in rows) == provenance["selected"]
    mask = args.campaign / "config/sieve_mask.json"
    blocked = {tuple(v) for v in json.loads(mask.read_text())["blocked"]} if mask.exists() else set()
    groups = defaultdict(list)
    for r in rows:
        groups[(r["rungroup"], int(r["foil"]), int(r["ndel"]))].append(r)
    norm = LogNorm(1, max(int(r["available"]) for r in rows))
    cmap = plt.get_cmap("viridis")
    out = args.data / "plots"
    out.mkdir(exist_ok=True)
    for (name, foil, delta), cells in sorted(groups.items()):
        fig, axes = plt.subplots(1, 2, figsize=(12, 10))
        fig.subplots_adjust(left=.07, right=.88, bottom=.25, top=.84, wspace=.25)
        for ax, field, title in zip(axes, ["available", "training"],
                                   ["Available labeled fit-input events", "Events admitted to SVD"]):
            values = {(int(r["xscol"]), int(r["yscol"])): int(r[field]) for r in cells}
            for ix in range(spec.nx):
                for iy in range(spec.ny):
                    pair = ix, iy
                    x, y = spec.ys(iy), spec.xs(ix)
                    value = values.get(pair, 0)
                    if pair in blocked and not value:
                        ax.text(x, y, "×", ha="center", va="center", color="#b22525", fontsize=20)
                    elif pair not in values:
                        ax.plot(x, y, ".", color="#bbc0c5", markersize=4)
                    else:
                        rgb = cmap(norm(value))[:3] if value else (1, 1, 1)
                        ax.add_patch(Circle((x, y), .62, facecolor=rgb,
                                            edgecolor="#b22525" if pair in blocked else "#68717a",
                                            linewidth=2 if pair in blocked else .6))
                        light = sum(a * b for a, b in zip(rgb, [.2126, .7152, .0722]))
                        ax.text(x, y, str(value), ha="center", va="center", fontsize=9,
                                color="white" if light < .5 else "#151515")
            positive = [v for k, v in values.items() if v and k not in blocked]
            ratio = f"{max(positive) / min(positive):.1f}×" if positive else "n/a"
            ax.set_title(title, fontsize=13, pad=15)
            ax.set(xlabel="Y sieve (cm)", ylabel="X sieve (cm)",
                   xlim=(spec.ys(0)-1, spec.ys(spec.ny-1)+1),
                   ylim=(spec.xs(0)-1.3, spec.xs(spec.nx-1)+1.3))
            ax.set_aspect("equal")
            ax.spines[["top", "right"]].set_visible(False)
            ax.text(.5, -.10, f"Total: {sum(values.values()):,} events\n"
                    f"Largest / smallest nonzero open hole: {ratio}", transform=ax.transAxes,
                    ha="center", va="top", fontsize=11, linespacing=1.6)
        cax = fig.add_axes([.91, .32, .017, .43])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                     label="Events per hole (shared log scale across campaign)")
        fig.suptitle(f"HMS 6.667 GeV · {name}\nFoil {foil} · delta slice {delta} · saved July 16, 2026 SVD fit",
                     fontsize=16, y=.97)
        av = sum(int(r["available"]) for r in cells)
        tr = sum(int(r["training"]) for r in cells)
        fig.text(.5, .065,
                 f"SVD retained {tr:,} / {av:,} events ({tr/av:.1%}) in this slice.\n"
                 "Selection reconstructed in input-tree order; campaign total and every logged Y-column count match.\n"
                 "Available = original GMM-cleaned TFit inputs, before SVD limits. No new core selection or allocation applied.\n"
                 "× = blocked and zero; red outline = events at a blocked position (included in totals); gray dot = no saved label.\n"
                 "Hollow 0 = zero selected events; nominal sieve positions. Original event-ID list was not saved.",
                 ha="center", va="center", fontsize=9, linespacing=1.6)
        fig.savefig(out / f"svd_{name}_foil{foil}_ndel{delta}.png", dpi=150, facecolor="white")
        plt.close(fig)
    print(f"Wrote {len(groups)} maps to {out}")


if __name__ == "__main__":
    main()
