#!/usr/bin/env python3
"""Compare joint labels, GMM survivors and verified historical SVD event clouds."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import uproot

from core_sample import read_joint, digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("rungroup", help="Name or unique prefix, e.g. rg01")
    args = parser.parse_args()
    campaign = args.campaign
    table, = campaign.glob("config/rungroups_*_inputs.tsv")
    with table.open() as f:
        matches = [r for r in csv.DictReader(f, delimiter="\t") if r["rungroup"].startswith(args.rungroup)]
    if len(matches) != 1:
        parser.error("Rungroup must match one saved setting")
    name, run = matches[0]["rungroup"], int(matches[0]["optics_id"])
    xp = campaign / f"04b_candidate_trees_x/root/XscolCandidates_{name}.root"
    yp = campaign / f"04a_candidate_trees_y/root/YscolCandidates_{name}.root"
    fp = campaign / f"06a_fit_ntuple/root/Optics_{run}_-1_fit_tree_gmm.root"
    data = campaign / "06b_svd_fit/diagnostics/balance"
    provenance = json.loads((data / "provenance.json").read_text())
    expected, = [v for k, v in provenance["input_sha256"].items() if k.endswith(fp.name)]
    assert digest(fp) == expected, "Fit input differs from verified historical selection"
    original, audit = read_joint(xp, yp, run)
    with uproot.open(fp) as f:
        gmm = f["TFit"].arrays(library="np")
    events, oi, gi = np.intersect1d(original["entry"], gmm["entry"], return_indices=True)
    assert len(events) == len(gmm["entry"]), "GMM events are not a unique subset of original labels"
    for before, after in [("xsieve", "xs"), ("ysieve", "ys"), ("foil", "foil"),
                          ("ndel", "ndel"), ("xscol", "xscol"), ("yscol", "yscol"), ("delta", "delta")]:
        assert np.allclose(original[before][oi], gmm[after][gi], rtol=0, atol=1e-9), before
    with np.load(data / "selected.npz") as selection:
        indices = selection[f"{run}_tree_index"]
        assert np.array_equal(gmm["entry"][indices], selection[f"{run}_entry"])
    svd = {k: v[indices] for k, v in gmm.items()}
    sources = [(original, "ysieve", "xsieve", "Original joint-labeled events"),
               (gmm, "ys", "xs", "After GMM cleanup"),
               (svd, "ys", "xs", "Events admitted to SVD")]
    # Fixed bins and limits preserve spatial and density comparisons across all figures.
    yedges = np.linspace(-7, 7, 281)
    xedges = np.linspace(-13, 13, 326)
    groups = sorted(set(zip(original["foil"].tolist(), original["ndel"].tolist())))
    histograms = {}
    for foil, delta in groups:
        for stage, (a, y, x, _) in enumerate(sources):
            keep = (a["foil"] == foil) & (a["ndel"] == delta)
            hist = np.histogram2d(a[y][keep], a[x][keep], bins=[yedges, xedges])[0]
            histograms[(foil, delta, stage)] = hist, int(keep.sum())
    norm = LogNorm(1, max(2, max(h.max() for h, _ in histograms.values())))
    out = data / "clouds"
    out.mkdir(exist_ok=True)
    summary = []
    for foil, delta in groups:
        fig, axes = plt.subplots(1, 3, figsize=(16, 9))
        fig.subplots_adjust(left=.055, right=.92, top=.79, bottom=.19, wspace=.25)
        totals = []
        for stage, (ax, (_, _, _, title)) in enumerate(zip(axes, sources)):
            h, n = histograms[(foil, delta, stage)]
            outside = n - int(h.sum())
            totals.append(n)
            mesh = ax.pcolormesh(yedges, xedges, np.ma.masked_equal(h.T, 0),
                                 cmap="viridis", norm=norm, rasterized=True)
            ax.set_title(f"{title}\nN = {n:,}; outside frame = {outside}", fontsize=12, pad=13)
            ax.set(xlabel="Y sieve (cm)", ylabel="X sieve (cm)",
                   xlim=(yedges[0], yedges[-1]), ylim=(xedges[0], xedges[-1]))
            ax.set_aspect("equal")
        cb = fig.add_axes([.94, .24, .014, .44])
        fig.colorbar(mesh, cax=cb, label="Events per bin (shared log scale, all five slices)")
        selected = (original["foil"] == foil) & (original["ndel"] == delta)
        low, high = original["delta_low"][selected][0], original["delta_high"][selected][0]
        fig.suptitle(f"HMS 6.667 GeV · {name} · foil {foil}\n"
                     f"Delta slice {delta}: [{low:g}, {high:g})% · historical SVD sample", fontsize=17, y=.95)
        fig.text(.5, .055,
                 "Identical sieve coordinates, axes and bins across stages. Bin size: 0.05 cm (Y) × 0.08 cm (X).\n"
                 "Original = unique joint X/Y labels; GMM survivors verified by event ID, labels and coordinates.\n"
                 "SVD membership reconstructed from the saved July 16 fit inputs and rules; all logged counts match. No core-selection cuts applied.",
                 ha="center", va="center", fontsize=10, linespacing=1.6)
        path = out / f"clouds_{name}_foil{foil}_ndel{delta}.png"
        fig.savefig(path, dpi=170, facecolor="white")
        plt.close(fig)
        summary.append(dict(foil=foil, ndel=delta, original=totals[0], gmm=totals[1], svd=totals[2]))
        print(path.name, totals)
    (out / f"{name}_provenance.json").write_text(json.dumps(dict(
        input_sha256={str(p): digest(p) for p in [xp, yp, fp, data / "selected.npz"]},
        script_sha256=digest(Path(__file__)), joint_excluded=len(audit),
        gmm_subset_verified=True, counts=summary,
        color_min=1, color_max=float(norm.vmax)), indent=2) + "\n")


if __name__ == "__main__":
    main()
