#!/usr/bin/env python3
"""Plot per-hole core populations from a saved core-sample tag."""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", type=Path, help="Saved 05c_core_sample tag directory")
    parser.add_argument("rungroup", help="Rungroup name or unique prefix, e.g. rg01")
    args = parser.parse_args()
    mask_path = args.sample.parent.parent / "config/sieve_mask.json"
    blocked = set()
    if mask_path.exists():
        mask = json.loads(mask_path.read_text())
        for pair in mask["blocked"]:
            if len(pair) != 2 or any(type(i) is not int or i < 0 for i in pair):
                parser.error("Blocked positions must be nonnegative integer [xscol, yscol] pairs")
            blocked.add(tuple(pair))
    with (args.sample / "tsv/counts.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    names = sorted({r["rungroup"] for r in rows if r["rungroup"].startswith(args.rungroup)})
    if len(names) != 1:
        parser.error("Rungroup must match exactly one saved rungroup")
    name = names[0]
    groups = defaultdict(dict)
    for row in rows:
        if row["rungroup"] == name:
            groups[int(row["foil"])].setdefault(int(row["ndel"]), []).append(row)
    # The selector's independent region summary must agree with fit + surplus.
    with (args.sample / "tsv/regions.tsv").open() as stream:
        regions = {(r["foil"], r["ndel"], r["xscol"], r["yscol"]): r
                   for r in csv.DictReader(stream, delimiter="\t") if r["rungroup"] == name}
    colors = dict(q25="#0057d9", median="#9700b8", mean="#101010",
                  cap15="#168457", cap20="#d11f38")
    styles = dict(q25="-", median=(0, (7, 3)), mean="-.",
                  cap15=(0, (4, 2)), cap20=":")
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False})
    for foil, slices in sorted(groups.items()):
        nrows = math.ceil((len(slices) + 1) / 2)
        fig, axes = plt.subplots(nrows, 2, figsize=(18, 6 * nrows))
        axes = axes.ravel()
        for ax, (ndel, cells) in zip(axes, sorted(slices.items())):
            for r in cells:
                region = regions[(r["foil"], r["ndel"], r["xscol"], r["yscol"])]
                expected = int(region["n_core_dev"]) if region["status"] == "accepted" else 0
                assert int(r["fit"]) + int(r["surplus"]) == expected
            available = np.array([int(r["fit"]) + int(r["surplus"]) for r in cells])
            training = np.array([int(r["fit"]) for r in cells])
            is_blocked = np.array([(int(r["xscol"]), int(r["yscol"])) in blocked for r in cells])
            keep = (available > 0) & ~is_blocked
            excluded = int(np.sum((available == 0) & ~is_blocked))
            nblocked = int(np.sum(is_blocked))
            blocked_training = int(training[is_blocked].sum())
            available, training = available[keep], training[keep]
            ax.set_title(f"Delta slice {ndel}  |  Open holes: {len(available)} with cores, {excluded} without\n"
                         f"Blocked positions with labels: {nblocked} (excluded)",
                         loc="left", fontsize=12, pad=12)
            ax.set_xlabel("Events per sieve hole")
            ax.set_ylabel("Number of sieve holes")
            if not len(available):
                ax.text(.5, .5, "No eligible cores", transform=ax.transAxes, ha="center")
                continue
            bins = np.arange(0, (int(available.max()) // 100 + 2) * 100, 100)
            counts, _ = np.histogram(available, bins)
            fitted, _ = np.histogram(training, bins)
            assert counts.sum() == fitted.sum() == len(available)
            ax.hist(available, bins=bins, color="#79b4db", alpha=.65, edgecolor="white")
            ax.hist(training, bins=bins, histtype="step", color="#c46a12", linewidth=2.2)
            q25, median = np.percentile(available, [25, 50])
            mean = available.mean()
            cap15, cap20 = math.floor(1.5 * q25), math.floor(2 * q25)
            handles = [Patch(facecolor="#79b4db", alpha=.65, label="Available cores"),
                       Line2D([], [], color="#c46a12", lw=2.2, label="Training cores")]
            for value, key, label in [(q25, "q25", "25th percentile"), (median, "median", "Median"),
                                      (mean, "mean", "Mean"), (cap15, "cap15", "1.5 × P25 cap"),
                                      (cap20, "cap20", "2 × P25 cap")]:
                ax.axvline(value, color=colors[key], linestyle=styles[key], linewidth=2)
                handles.append(Line2D([], [], color=colors[key], ls=styles[key], lw=2,
                                      label=f"{label}: {value:,.1f}"))
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -.20),
                      frameon=False, fontsize=10, ncol=3, handlelength=3.5,
                      columnspacing=1.5, borderaxespad=0)
            lines = [f"Available {available.sum():,}  |  Actual training {training.sum():,}",
                     f"P25 {q25:,.1f}  |  Median {median:,.0f}  |  Mean {mean:,.0f}"]
            for label, cap in [("1.5", cap15), ("2", cap20)]:
                retained = np.minimum(available, cap).sum()
                lines.append(f"{label} × P25: cap {cap:,} → {retained:,} events ({retained / available.sum():.0%}); "
                             f"{np.sum(available > cap)} holes capped")
            ax.text(.98, .97, "\n".join(lines), transform=ax.transAxes, va="top", ha="right",
                    fontsize=10, linespacing=1.6,
                    bbox=dict(facecolor="white", edgecolor="#dddddd", alpha=.94, pad=7))
            ax.set_xlim(0, bins[-1])
            ax.set_ylim(0, max(counts.max(), fitted.max()) * 1.55 + 1)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.grid(axis="y", alpha=.15)
            ax.set_axisbelow(True)
            print(f"foil={foil} ndel={ndel}: holes={len(available)} available={available.sum()} "
                  f"training={training.sum()} cap1.5={cap15}/{np.minimum(available, cap15).sum()} "
                  f"cap2={cap20}/{np.minimum(available, cap20).sum()} blocked={nblocked} "
                  f"blocked_training={blocked_training}")
        key = axes[len(slices)]
        key.axis("off")
        key.text(.02, .95,
                 "Each histogram entry represents one open hole with cores.\n"
                 "Bin width: 100 events; horizontal ranges vary by slice.\n"
                 "Protected holdout is excluded. P25 uses positive core counts.\n"
                 "Cap totals apply only the proposed cap to available cores;\n"
                 "they do not apply the current 80% rule or campaign budget.\n"
                 "Caps are rounded down; no allocation has been changed.\n"
                 "Blocked positions are excluded using campaign config/sieve_mask.json.\n"
                 "Open = not listed as blocked; absent labels are not counted.\n"
                 + ("No mask supplied: blocked geometry is unverified." if not mask_path.exists() else
                    "Mask pairs are [X label, Y label], not sieve coordinates."),
                 fontsize=10, va="top", linespacing=1.4)
        for ax in axes[len(slices) + 1:]:
            ax.axis("off")
        fig.suptitle(f"{args.sample.parent.parent.name} · {name} · foil {foil}\nPer-hole core populations · {args.sample.name}",
                     fontsize=18, y=.99)
        fig.tight_layout(rect=(0, 0, 1, .95), h_pad=5, w_pad=3)
        out = args.sample / "plots" / f"counts_{name}_foil{foil}.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=150, facecolor="white")
        plt.close(fig)
        print(out)


if __name__ == "__main__":
    main()
