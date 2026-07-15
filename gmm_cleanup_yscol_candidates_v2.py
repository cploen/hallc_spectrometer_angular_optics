#!/usr/bin/env python3
"""
GMM cleanup for YFP/YPFP yscol candidate events already extracted into TYCand.

This script does NOT reopen or reapply TCutG files.  The candidate tree is assumed
already selected by:
  baseline cuts + ytar/delta foil cut + YFP/YPFP yscol band cut.

It runs scikit-learn GaussianMixture in projected sieve space only:
  X = [xsieve, ysieve]

Default behavior for one rungroup/foil/delta slice:
  - derive the candidate-tree input from the campaign directory
  - process all yscol values present in the candidate tree
  - remove columns with N < min_events
  - choose GMM component count by BIC, bounded by max_components <= 9
  - reject lowest-density events using keep_frac
  - write ROOT/CSV/TSV/PDF outputs into 05a_gmm_cleanup_y subdirectories

Example:
  python3 gmm_cleanup_yscol_candidates_v2.py \
    --rungroup rg04_theta12p395_foilpm3 \
    --campaign HMS_5p878GeV \
    --foil 0 --ndel 1 \
    --min-events 30 --max-components 9 --keep-frac 0.95
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    from sklearn.mixture import GaussianMixture
except Exception as exc:
    raise RuntimeError(
        "scikit-learn is required. Try: python3 -m pip install --user scikit-learn"
    ) from exc

try:
    import ROOT
except Exception as exc:
    raise RuntimeError("PyROOT is required to read/write ROOT files.") from exc


@dataclass
class EventRecord:
    entry: int
    run: int
    foil: int
    ndel: int
    yscol: int
    xfp: float
    xpfp: float
    yfp: float
    ypfp: float
    xsieve: float
    ysieve: float
    delta: float = float("nan")
    delta_low: float = float("nan")
    delta_high: float = float("nan")
    ytar: float = float("nan")
    xptar: float = float("nan")
    yptar: float = float("nan")
    gmm_keep: int = 0
    gmm_score: float = float("nan")
    gmm_component: int = -1
    gmm_action: str = "unset"


def branch_names(tree) -> set[str]:
    return {br.GetName() for br in tree.GetListOfBranches()}


def find_branch(tree, aliases: Iterable[str], required: bool = True) -> str | None:
    names = branch_names(tree)
    for name in aliases:
        if name in names:
            return name
    if required:
        raise KeyError(f"None of these branches found: {list(aliases)}\nAvailable: {sorted(names)}")
    return None


def val(tree, branch: str | None, default=float("nan")):
    if branch is None:
        return default
    return getattr(tree, branch)


def read_candidate_tree(args) -> List[EventRecord]:
    fin = ROOT.TFile.Open(str(args.input_root), "READ")
    if not fin or fin.IsZombie():
        raise FileNotFoundError(f"Could not open input ROOT file: {args.input_root}")

    tree = fin.Get(args.tree)
    if not tree:
        raise KeyError(f"Could not find tree '{args.tree}' in {args.input_root}")

    br_entry = find_branch(tree, ["entry", "Entry$"], required=False)
    br_run   = find_branch(tree, ["run", "Run", "run_out"], required=False)
    br_foil  = find_branch(tree, ["foil", "nfoil", "foil_out"], required=True)
    br_ndel  = find_branch(tree, ["ndel", "ndel_out"], required=True)
    br_yscol = find_branch(tree, ["yscol", "assigned_yscol", "yscol_out"], required=True)

    br_xfp   = find_branch(tree, ["xfp", "H.dc.x_fp"], required=True)
    br_xpfp  = find_branch(tree, ["xpfp", "H.dc.xp_fp"], required=True)
    br_yfp   = find_branch(tree, ["yfp", "H.dc.y_fp"], required=True)
    br_ypfp  = find_branch(tree, ["ypfp", "H.dc.yp_fp"], required=True)
    br_xs    = find_branch(tree, ["xsieve", "xs", "H.extcor.xsieve"], required=True)
    br_ys    = find_branch(tree, ["ysieve", "ys", "H.extcor.ysieve"], required=True)

    br_delta      = find_branch(tree, ["delta", "H.gtr.dp"], required=False)
    br_delta_low  = find_branch(tree, ["delta_low"], required=True)
    br_delta_high = find_branch(tree, ["delta_high"], required=True)
    br_ytar       = find_branch(tree, ["ytar", "H.gtr.y"], required=False)
    br_xptar = find_branch(tree, ["xptar", "H.gtr.th"], required=False)
    br_yptar = find_branch(tree, ["yptar", "H.gtr.ph"], required=False)

    records: List[EventRecord] = []
    nentries = tree.GetEntries()
    for i in range(nentries):
        tree.GetEntry(i)
        foil = int(val(tree, br_foil))
        ndel = int(val(tree, br_ndel))
        yscol = int(val(tree, br_yscol))

        if foil != args.foil or ndel != args.ndel:
            continue
        if args.yscol is not None and yscol != args.yscol:
            continue

        run = int(val(tree, br_run, -1)) if br_run else -1

        entry = int(val(tree, br_entry, i)) if br_entry else i
        records.append(EventRecord(
            entry=entry,
            run=run,
            foil=foil,
            ndel=ndel,
            yscol=yscol,
            xfp=float(val(tree, br_xfp)),
            xpfp=float(val(tree, br_xpfp)),
            yfp=float(val(tree, br_yfp)),
            ypfp=float(val(tree, br_ypfp)),
            xsieve=float(val(tree, br_xs)),
            ysieve=float(val(tree, br_ys)),
            delta=float(val(tree, br_delta)),
            delta_low=float(val(tree, br_delta_low)),
            delta_high=float(val(tree, br_delta_high)),
            ytar=float(val(tree, br_ytar)),
            xptar=float(val(tree, br_xptar)),
            yptar=float(val(tree, br_yptar)),
        ))

    fin.Close()
    return records


def fit_one_yscol(records: List[EventRecord], args) -> Dict[str, float | int | str]:
    n = len(records)
    summary: Dict[str, float | int | str] = {
        "yscol": records[0].yscol if records else -999,
        "n_before": n,
        "n_after": 0,
        "keep_frac_actual": 0.0,
        "action": "empty",
        "n_components": 0,
        "bic": float("nan"),
        "score_threshold": float("nan"),
    }

    if n < args.min_events:
        for r in records:
            r.gmm_keep = 0
            r.gmm_action = "removed_low_stat"
        summary.update(action="removed_low_stat")
        return summary

    X = np.array([[r.xsieve, r.ysieve] for r in records], dtype=float)
    finite = np.isfinite(X).all(axis=1)
    if not np.all(finite):
        # Non-finite points are not safe SVD anchors.
        for r, ok in zip(records, finite):
            if not ok:
                r.gmm_keep = 0
                r.gmm_action = "removed_nonfinite"
        fit_records = [r for r, ok in zip(records, finite) if ok]
        X = X[finite]
    else:
        fit_records = records

    if len(fit_records) < args.min_events:
        for r in fit_records:
            r.gmm_keep = 0
            r.gmm_action = "removed_low_stat_after_finite"
        summary.update(action="removed_low_stat_after_finite")
        return summary

    # Upper bound: physically no more than 9 xsieve hole groups in one yscol projection.
    # Require several points per component to avoid one/two-point nuisance blobs.
    max_k = min(args.max_components, max(1, len(fit_records) // args.min_points_per_component))

    best_gmm = None
    best_bic = math.inf
    best_k = 0
    for k in range(1, max_k + 1):
        gmm = GaussianMixture(
            n_components=k,
            covariance_type="full",
            reg_covar=args.reg_covar,
            random_state=args.random_state,
            n_init=args.n_init,
        )
        gmm.fit(X)
        bic = float(gmm.bic(X))
        if bic < best_bic:
            best_bic = bic
            best_gmm = gmm
            best_k = k

    assert best_gmm is not None
    scores = best_gmm.score_samples(X)
    labels = best_gmm.predict(X)

    threshold = float(np.quantile(scores, 1.0 - args.keep_frac))
    keep = scores >= threshold

    for r, sc, lab, ok in zip(fit_records, scores, labels, keep):
        r.gmm_score = float(sc)
        r.gmm_component = int(lab)
        r.gmm_keep = int(ok)
        r.gmm_action = "gmm_cleaned"

    n_after = sum(r.gmm_keep for r in records)
    summary.update(
        action="gmm_cleaned",
        n_components=best_k,
        bic=best_bic,
        score_threshold=threshold,
        n_after=int(n_after),
        keep_frac_actual=(float(n_after) / float(n) if n else 0.0),
    )
    return summary


def process_records(records: List[EventRecord], args) -> List[Dict[str, float | int | str]]:
    summaries = []
    for yscol in sorted({r.yscol for r in records}):
        sub = [r for r in records if r.yscol == yscol]
        summaries.append(fit_one_yscol(sub, args))
    return summaries


def write_root(records: List[EventRecord], summaries: List[Dict[str, float | int | str]], root_path: Path, args) -> None:
    import array

    fout = ROOT.TFile(str(root_path), "RECREATE")
    t = ROOT.TTree("GMMClean", "GMM event mask for yscol candidates in one rungroup/foil/delta slice")
    ROOT.TNamed("rungroup", str(args.rungroup)).Write()

    fields_i = {
        "entry": array.array("i", [0]),
        "run": array.array("i", [0]),
        "foil": array.array("i", [0]),
        "ndel": array.array("i", [0]),
        "yscol": array.array("i", [0]),
        "gmm_keep": array.array("i", [0]),
        "gmm_component": array.array("i", [0]),
    }
    fields_d = {
        "xfp": array.array("d", [0.0]),
        "xpfp": array.array("d", [0.0]),
        "yfp": array.array("d", [0.0]),
        "ypfp": array.array("d", [0.0]),
        "xsieve": array.array("d", [0.0]),
        "ysieve": array.array("d", [0.0]),
        "delta": array.array("d", [0.0]),
        "ytar": array.array("d", [0.0]),
        "xptar": array.array("d", [0.0]),
        "yptar": array.array("d", [0.0]),
        "gmm_score": array.array("d", [0.0]),
    }

    for name, arr in fields_i.items():
        t.Branch(name, arr, f"{name}/I")
    for name, arr in fields_d.items():
        t.Branch(name, arr, f"{name}/D")

    for r in records:
        fields_i["entry"][0] = r.entry
        fields_i["run"][0] = r.run
        fields_i["foil"][0] = r.foil
        fields_i["ndel"][0] = r.ndel
        fields_i["yscol"][0] = r.yscol
        fields_i["gmm_keep"][0] = r.gmm_keep
        fields_i["gmm_component"][0] = r.gmm_component
        fields_d["xfp"][0] = r.xfp
        fields_d["xpfp"][0] = r.xpfp
        fields_d["yfp"][0] = r.yfp
        fields_d["ypfp"][0] = r.ypfp
        fields_d["xsieve"][0] = r.xsieve
        fields_d["ysieve"][0] = r.ysieve
        fields_d["delta"][0] = r.delta
        fields_d["ytar"][0] = r.ytar
        fields_d["xptar"][0] = r.xptar
        fields_d["yptar"][0] = r.yptar
        fields_d["gmm_score"][0] = r.gmm_score
        t.Fill()

    ts = ROOT.TTree("GMMSummary", "GMM cleanup summary by yscol")
    s_yscol = array.array("i", [0])
    s_n_before = array.array("i", [0])
    s_n_after = array.array("i", [0])
    s_n_components = array.array("i", [0])
    s_bic = array.array("d", [0.0])
    s_thr = array.array("d", [0.0])
    s_keep_frac = array.array("d", [0.0])
    for name, arr, spec in [
        ("yscol", s_yscol, "yscol/I"),
        ("n_before", s_n_before, "n_before/I"),
        ("n_after", s_n_after, "n_after/I"),
        ("n_components", s_n_components, "n_components/I"),
        ("bic", s_bic, "bic/D"),
        ("score_threshold", s_thr, "score_threshold/D"),
        ("keep_frac_actual", s_keep_frac, "keep_frac_actual/D"),
    ]:
        ts.Branch(name, arr, spec)

    for s in summaries:
        s_yscol[0] = int(s["yscol"])
        s_n_before[0] = int(s["n_before"])
        s_n_after[0] = int(s["n_after"])
        s_n_components[0] = int(s["n_components"])
        s_bic[0] = float(s["bic"])
        s_thr[0] = float(s["score_threshold"])
        s_keep_frac[0] = float(s["keep_frac_actual"])
        ts.Fill()

    fout.Write()
    fout.Close()


def write_csv(records: List[EventRecord], summaries: List[Dict[str, float | int | str]], csv_path: Path, summary_path: Path, args) -> None:
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "rungroup", "entry", "run", "foil", "ndel", "yscol",
            "xfp", "xpfp", "yfp", "ypfp", "xsieve", "ysieve",
            "delta", "ytar", "xptar", "yptar",
            "gmm_keep", "gmm_score", "gmm_component", "gmm_action",
        ])
        for r in records:
            w.writerow([
                args.rungroup, r.entry, r.run, r.foil, r.ndel, r.yscol,
                r.xfp, r.xpfp, r.yfp, r.ypfp, r.xsieve, r.ysieve,
                r.delta, r.ytar, r.xptar, r.yptar,
                r.gmm_keep, r.gmm_score, r.gmm_component, r.gmm_action,
            ])

    delta_low, delta_high = get_delta_bounds(records)
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "rungroup", "foil", "ndel", "delta_low", "delta_high",
            "yscol", "n_before", "n_after", "keep_frac_actual",
            "action", "n_components", "bic", "score_threshold",
        ])
        for s in summaries:
            w.writerow([
                args.rungroup, args.foil, args.ndel, delta_low, delta_high,
                s["yscol"], s["n_before"], s["n_after"], s["keep_frac_actual"],
                s["action"], s["n_components"], s["bic"], s["score_threshold"],
            ])


def get_delta_bounds(records: List[EventRecord]) -> Tuple[float, float]:
    lows = sorted({
        round(r.delta_low, 12)
        for r in records
        if math.isfinite(r.delta_low)
    })
    highs = sorted({
        round(r.delta_high, 12)
        for r in records
        if math.isfinite(r.delta_high)
    })

    if len(lows) != 1 or len(highs) != 1:
        raise RuntimeError(
            f"Expected one delta interval; found lows={lows}, highs={highs}"
        )

    return float(lows[0]), float(highs[0])


def format_delta_edge(value: float) -> str:
    magnitude = f"{abs(value):g}".replace(".", "p")
    return f"m{magnitude}" if value < 0 else magnitude


def delta_tag(records: List[EventRecord]) -> str:
    low, high = get_delta_bounds(records)
    return (
        f"delta_{format_delta_edge(low)}"
        f"_to_{format_delta_edge(high)}"
    )


def delta_label(records: List[EventRecord]) -> str:
    low, high = get_delta_bounds(records)
    return f"δ = {low:g}% to {high:g}%"


def draw_ysieve_guides(ax, args) -> None:
    """Draw ROOT/manual-optics ysieve guide lines.

    This must match plot_yfp_cuts.C / fit_opt_matrix.C:
        pos = (nys - 4) * 0.6 * 2.54
    Do not make this runtime-relative.
    """
    y_bottom = float(args.xsieve_min)
    y_top = float(args.xsieve_max)
    label_y = y_bottom - 0.035 * (y_top - y_bottom)

    for nys in range(9):
        ypos = (nys - 4) * 0.6 * 2.54
        ax.axvline(ypos, linewidth=0.8, alpha=0.85, color="red")
        ax.text(
            ypos, label_y, str(nys),
            ha="center", va="top", fontsize=7, color="red", clip_on=False,
        )


def style_sieve_axis(ax, args) -> None:
    """Keep sieve diagnostics in a ROOT-like, equal-scale view."""
    ax.set_xlabel("ysieve (cm)")
    ax.set_ylabel("xsieve (cm)")
    ax.set_xlim(float(args.ysieve_min), float(args.ysieve_max))
    ax.set_ylim(float(args.xsieve_min), float(args.xsieve_max))
    ax.set_aspect("equal", adjustable="box")
    ax.minorticks_on()
    draw_ysieve_guides(ax, args)


def make_pdf(records: List[EventRecord], summaries: List[Dict[str, float | int | str]], pdf_path: Path, args) -> None:
    with PdfPages(pdf_path) as pdf:
        # Page 1: all color-coded yscol projections, before and after.
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for yscol in sorted({r.yscol for r in records}):
            sub = [r for r in records if r.yscol == yscol]
            axes[0].scatter([r.ysieve for r in sub], [r.xsieve for r in sub], s=7, label=f"y{yscol}")
            kept = [r for r in sub if r.gmm_keep]
            if kept:
                axes[1].scatter([r.ysieve for r in kept], [r.xsieve for r in kept], s=7, label=f"y{yscol}")
        axes[0].set_title("All yscol projections BEFORE GMM")
        axes[1].set_title("All yscol projections AFTER GMM")
        for ax in axes:
            style_sieve_axis(ax, args)
            ax.legend(fontsize=7, loc="best")
        fig.suptitle(
            f"{args.rungroup} — foil {args.foil}, {delta_label(records)}"
        )
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Per-yscol diagnostic pages.
        for s in summaries:
            yscol = int(s["yscol"])
            sub = [r for r in records if r.yscol == yscol]
            keep = np.array([bool(r.gmm_keep) for r in sub])
            ypfp = np.array([r.ypfp for r in sub])
            yfp = np.array([r.yfp for r in sub])
            xs = np.array([r.xsieve for r in sub])
            ys = np.array([r.ysieve for r in sub])

            fig = plt.figure(figsize=(11, 8.5))
            ax1 = fig.add_subplot(2, 2, 1)
            ax2 = fig.add_subplot(2, 2, 2)
            ax3 = fig.add_subplot(2, 2, 3)
            ax4 = fig.add_subplot(2, 2, 4)
            ax4.axis("off")

            ax1.scatter(ypfp, yfp, s=8)
            ax1.set_title(f"YFP/YPFP candidate reference: yscol {yscol}")
            ax1.set_xlabel("ypfp")
            ax1.set_ylabel("yfp")

            ax2.scatter(ys, xs, s=8)
            ax2.set_title("BEFORE GMM: projected sieve space")
            style_sieve_axis(ax2, args)

            if len(sub) > 0:
                ax3.scatter(ys[~keep], xs[~keep], s=10, marker="x", label="rejected")
                ax3.scatter(ys[keep], xs[keep], s=8, label="kept")
            ax3.set_title("AFTER GMM cleanup")
            style_sieve_axis(ax3, args)
            ax3.legend(fontsize=8, loc="best")

            bic = float(s["bic"])
            thr = float(s["score_threshold"])
            txt = [
                f"rungroup = {args.rungroup}",
                f"foil = {args.foil}   ndel = {args.ndel}   {delta_label(sub)}",
                f"yscol = {yscol}",
                "",
                f"N before = {s['n_before']}",
                f"N after  = {s['n_after']}",
                f"keep frac actual = {float(s['keep_frac_actual']):.3f}",
                f"action = {s['action']}",
                "",
                f"min_events = {args.min_events}",
                f"max_components = {args.max_components}",
                f"chosen components = {s['n_components']}",
                f"BIC = {bic:.3g}" if math.isfinite(bic) else "BIC = n/a",
                f"score threshold = {thr:.3g}" if math.isfinite(thr) else "score threshold = n/a",
            ]
            ax4.text(0.02, 0.98, "\n".join(txt), va="top", ha="left", fontsize=11, family="monospace")
            fig.suptitle("yscol candidate → GMM sieve-space cleanup", fontsize=14)
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            pdf.savefig(fig)
            plt.close(fig)


def output_paths(args, records: List[EventRecord]) -> Tuple[Path, Path, Path, Path]:
    outdir = Path(args.outdir)

    root_dir = outdir / "root"
    csv_dir = outdir / "csv"
    tsv_dir = outdir / "tsv"
    plot_dir = outdir / "plots"

    for directory in (root_dir, csv_dir, tsv_dir, plot_dir):
        directory.mkdir(parents=True, exist_ok=True)

    selection = "allY" if args.yscol is None else f"yscol{args.yscol}"
    tag = (
        f"{args.rungroup}_foil{args.foil}_"
        f"{delta_tag(records)}_{selection}"
    )

    return (
        root_dir / f"gmm_clean_{tag}.root",
        csv_dir / f"gmm_clean_{tag}.csv",
        tsv_dir / f"gmm_clean_{tag}_summary.tsv",
        plot_dir / f"gmm_clean_{tag}.pdf",
    )


def main() -> None:
    p = argparse.ArgumentParser(description="GMM cleanup for TYCand yscol candidates in xsieve/ysieve space.")
    p.add_argument("--rungroup", required=True, help="Full rungroup tag used for campaign I/O and labels.")
    p.add_argument("--campaign", type=Path, required=True, help="Campaign path relative to the current project-root directory, e.g. HMS_5p878GeV.")
    p.add_argument("--foil", type=int, required=True)
    p.add_argument("--ndel", type=int, required=True)
    p.add_argument("--yscol", type=int, default=None, help="Optional: process only one yscol. Default processes all present yscols.")
    p.add_argument("--input-root", type=Path, default=None, help="Override candidate ROOT input.")
    p.add_argument("--tree", default="TYCand")
    p.add_argument("--outdir", type=Path, default=None, help="Override 05a_gmm_cleanup_y output directory.")

    p.add_argument("--min-events", type=int, default=30, help="Remove candidate columns with fewer than this many events.")
    p.add_argument("--max-components", type=int, default=9, help="Upper bound on GMM components per yscol; must be <= 9.")
    p.add_argument("--keep-frac", type=float, default=0.95, help="Fraction retained after GMM log-density cut.")
    p.add_argument("--min-points-per-component", type=int, default=5, help="Guard against tiny nuisance components.")

    p.add_argument("--reg-covar", type=float, default=1e-6)
    p.add_argument("--n-init", type=int, default=5)
    p.add_argument("--random-state", type=int, default=13)

    # ROOT-style sieve diagnostic view. These are plotting-only; they do not affect GMM.
    p.add_argument("--ysieve-min", type=float, default=-6.8)
    p.add_argument("--ysieve-max", type=float, default=6.8)
    p.add_argument("--xsieve-min", type=float, default=-12.5)
    p.add_argument("--xsieve-max", type=float, default=12.5)
    p.add_argument("--guide-spacing", type=float, default=1.524, help="ysieve column-guide spacing in cm")
    p.add_argument("--guide-origin", type=float, default=0.0, help="ysieve position of the center guide column")
    p.add_argument("--guide-center-index", type=int, default=4, help="column index located at guide-origin")
    p.add_argument("--guide-first", type=int, default=0)
    p.add_argument("--guide-last", type=int, default=8)
    p.add_argument("--no-guide-labels", action="store_false", dest="show_guide_labels")
    p.set_defaults(show_guide_labels=True)

    args = p.parse_args()

    project_root = Path.cwd().resolve()
    args.campaign_dir = (project_root / args.campaign).resolve()

    if not args.campaign_dir.is_dir():
        raise FileNotFoundError(
            f"Campaign directory does not exist: {args.campaign_dir}"
        )

    if args.input_root is None:
        args.input_root = (
            args.campaign_dir
            / "04a_candidate_trees_y"
            / "root"
            / f"YscolCandidates_{args.rungroup}.root"
        )

    if args.outdir is None:
        args.outdir = args.campaign_dir / "05a_gmm_cleanup_y"

    if not (0.0 < args.keep_frac <= 1.0):
        raise ValueError("--keep-frac must be in (0, 1].")
    if args.max_components > 9:
        raise ValueError("--max-components should not exceed 9.")
    if args.min_points_per_component < 2:
        raise ValueError("--min-points-per-component should be at least 2.")

    records = read_candidate_tree(args)
    if not records:
        raise RuntimeError(
            f"No candidate events found for rungroup={args.rungroup} "
            f"foil={args.foil} ndel={args.ndel} yscol={args.yscol}"
        )

    summaries = process_records(records, args)
    root_path, csv_path, summary_path, pdf_path = output_paths(args, records)
    write_root(records, summaries, root_path, args)
    write_csv(records, summaries, csv_path, summary_path, args)
    make_pdf(records, summaries, pdf_path, args)

    print("GMM cleanup complete")
    print(f"  input       : {args.input_root}")
    print(f"  rungroup     : {args.rungroup}")
    print(f"  foil/ndel    : {args.foil}/{args.ndel}")
    print(f"  delta slice  : {delta_label(records)}")
    print(f"  yscols      : {sorted({r.yscol for r in records})}")
    for s in summaries:
        print(f"  y{s['yscol']}: action={s['action']} N={s['n_before']} -> {s['n_after']} components={s['n_components']}")
    print(f"  ROOT output : {root_path}")
    print(f"  CSV output  : {csv_path}")
    print(f"  Summary TSV : {summary_path}")
    print(f"  PDF output  : {pdf_path}")


if __name__ == "__main__":
    main()
