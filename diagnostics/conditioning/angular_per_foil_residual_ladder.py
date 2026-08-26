#!/usr/bin/env python3
"""Delta-style angular residual ladders for every nominal foil.

For every prefix N=1..210, independently solve the unit-column design matrix
directly with SVD, transform the coefficients back to the HMS convention, and
measure residual RMS separately for selected training events and eligible
events excluded by the fit caps.  The current unscaled-XTX ladder is evaluated
from its saved coefficients for a method-matched comparison.
"""

from __future__ import annotations

import argparse
import csv
import gc
import math
from array import array
from pathlib import Path

import numpy as np
import ROOT

from angular_term_ladder import load_problem
from angular_solver_comparison import load_heldout_problem
from preliminary_angular_conditioning import (
    read_campaign_rows,
    read_optics_metadata,
    selected_indices,
)


TARGET_NAMES = ("xptar", "ytar", "yptar")
TARGET_LABELS = ("x'_{tar}", "y_{tar}", "y'_{tar}")
TARGET_UNITS = ("mrad", "cm", "mrad")
PHYSICAL_SCALES = np.asarray([1000.0, 100.0, 1000.0])


def selected_foil_labels(campaign, metadata_path, file_id, nfit_max):
    settings = read_campaign_rows(
        next((campaign / "config").glob("rungroups_*_inputs.tsv"))
    )
    metadata = read_optics_metadata(
        metadata_path, {setting[1] for setting in settings}
    )
    labels = []
    global_count = 0
    branches = ["delta", "ztarT", "ysT"]
    for _, optics_id in settings:
        path = (
            campaign
            / "06a_fit_ntuple/root"
            / f"Optics_{optics_id}_{file_id}_fit_tree_gmm.root"
        )
        arrays = ROOT.RDataFrame("TFit", str(path)).AsNumpy(branches)
        selected, _ = selected_indices(
            arrays, *metadata[optics_id], global_count, nfit_max
        )
        labels.append(arrays["ztarT"][selected])
        global_count += len(selected)
    return np.concatenate(labels)


def rescaled_direct_ladder(design, targets):
    norms = np.linalg.norm(design, axis=0)
    # A single QR factorization supplies every prefix: later columns do not
    # change the QR factorization of earlier columns.  Scale the leading R
    # block before its small SVD so this is explicitly the unit-column solve.
    q, r = np.linalg.qr(design, mode="reduced")
    q_transpose_targets = q.T @ targets
    del q
    gc.collect()
    eps = np.finfo(np.float64).eps
    maximum_terms = design.shape[1]
    solutions = np.zeros((maximum_terms, maximum_terms, targets.shape[1]))
    ranks = np.zeros(maximum_terms, dtype=np.int32)
    conditions = np.zeros(maximum_terms)
    for nterms in range(1, maximum_terms + 1):
        scaled_r = r[:nterms, :nterms] / norms[:nterms][None, :]
        u, singular, vh = np.linalg.svd(scaled_r, full_matrices=False)
        keep = singular > eps * singular[0]
        scaled_coordinates = vh.T @ (
            np.where(keep, 1.0 / singular, 0.0)[:, None]
            * (u.T @ q_transpose_targets[:nterms])
        )
        solutions[nterms - 1, :nterms] = (
            scaled_coordinates / norms[:nterms, None]
        )
        ranks[nterms - 1] = int(np.count_nonzero(keep))
        conditions[nterms - 1] = singular[0] / singular[-1]
    return solutions, ranks, conditions, norms


def sufficient_statistics(design, targets, foil_labels, foil_values):
    result = {}
    for foil in foil_values:
        mask = np.isclose(foil_labels, foil, atol=1e-8)
        x = design[mask]
        y = targets[mask]
        result[foil] = {
            "count": int(np.count_nonzero(mask)),
            "gram": x.T @ x,
            "rhs": x.T @ y,
            "target_norm2": np.sum(y * y, axis=0),
        }
    return result


def residual_rms(statistics, coefficients, nterms):
    d = coefficients[:nterms]
    residual_sum2 = (
        statistics["target_norm2"]
        - 2.0 * np.sum(d * statistics["rhs"][:nterms], axis=0)
        + np.einsum(
            "ik,ij,jk->k",
            d,
            statistics["gram"][:nterms, :nterms],
            d,
        )
    )
    return (
        np.sqrt(np.maximum(residual_sum2, 0.0) / statistics["count"])
        * PHYSICAL_SCALES
    )


def graph(x, y, color, marker, line_style=1):
    item = ROOT.TGraph(len(x), array("d", x), array("d", y))
    item.SetLineColor(color)
    item.SetMarkerColor(color)
    item.SetMarkerStyle(marker)
    item.SetMarkerSize(0.48)
    item.SetLineWidth(2)
    item.SetLineStyle(line_style)
    return item


def foil_name(value):
    if abs(value) < 1e-8:
        return "foil_0cm"
    return f"foil_{'p' if value > 0 else 'm'}{abs(value):g}cm"


def minimum_and_degradation_onset(values, minimum_n=35, fraction=0.05, window=5):
    start = minimum_n - 1
    minimum_index = start + int(np.argmin(values[start:]))
    threshold = values[minimum_index] * (1.0 + fraction)
    onset = None
    for index in range(minimum_index + 1, len(values) - window + 1):
        if np.all(values[index : index + window] > threshold):
            onset = index + 1
            break
    return minimum_index + 1, onset


def plot_ladder(
    output,
    method,
    foil,
    target_index,
    n_values,
    training,
    heldout,
    training_count,
    heldout_count,
    reference_n,
    plot_min_n=35,
):
    method_title = (
        "rescaled direct X" if method == "rescaled_direct_X" else "current unscaled X^{T}X"
    )
    target_label = TARGET_LABELS[target_index]
    units = TARGET_UNITS[target_index]
    minimum_n, degradation_n = minimum_and_degradation_onset(
        heldout, minimum_n=plot_min_n
    )
    minimum_index = minimum_n - 1
    reference = heldout[reference_n - 1]
    relative = 100.0 * (heldout / reference - 1.0)
    plotted = n_values >= plot_min_n

    canvas = ROOT.TCanvas(
        f"c_{method}_{foil_name(foil)}_{TARGET_NAMES[target_index]}",
        "Per-foil residual ladder",
        1100,
        900,
    )
    unique = f"{method}_{foil_name(foil)}_{TARGET_NAMES[target_index]}"
    top = ROOT.TPad(f"top_{unique}", "top", 0.0, 0.34, 1.0, 1.0)
    bottom = ROOT.TPad(f"bottom_{unique}", "bottom", 0.0, 0.0, 1.0, 0.34)
    top.SetBottomMargin(0.02)
    bottom.SetTopMargin(0.04)
    bottom.SetBottomMargin(0.26)
    top.Draw()
    bottom.Draw()

    top.cd()
    training_graph = graph(
        n_values[plotted], training[plotted], ROOT.kBlue + 1, 20
    )
    heldout_graph = graph(
        n_values[plotted], heldout[plotted], ROOT.kRed + 1, 20
    )
    multi = ROOT.TMultiGraph()
    multi.Add(training_graph, "LP")
    multi.Add(heldout_graph, "LP")
    multi.SetTitle(
        f"HMS {target_label} residuals by N: foil z={foil:g} cm, {method_title};"
        f"Number of fitted terms N;Residual RMS ({units})"
    )
    multi.Draw("A")
    multi.GetXaxis().SetLabelSize(0.0)
    ROOT.gPad.Update()
    y_min = ROOT.gPad.GetUymin()
    y_max = ROOT.gPad.GetUymax()
    minimum_line = ROOT.TLine(minimum_n, y_min, minimum_n, y_max)
    minimum_line.SetLineStyle(2)
    minimum_line.SetLineColor(ROOT.kGray + 2)
    minimum_line.SetLineWidth(2)
    minimum_line.Draw()
    degradation_line = None
    if degradation_n is not None:
        degradation_line = ROOT.TLine(degradation_n, y_min, degradation_n, y_max)
        degradation_line.SetLineStyle(3)
        degradation_line.SetLineColor(ROOT.kRed + 1)
        degradation_line.SetLineWidth(2)
        degradation_line.Draw()
    legend = ROOT.TLegend(0.13, 0.72, 0.52, 0.87)
    legend.AddEntry(training_graph, f"Selected training events ({training_count})", "lp")
    legend.AddEntry(heldout_graph, f"Eligible but not selected ({heldout_count})", "lp")
    legend.Draw()
    annotation = ROOT.TPaveText(0.57, 0.67, 0.94, 0.88, "NDC")
    annotation.SetFillColor(ROOT.kWhite)
    annotation.SetBorderSize(1)
    annotation.SetTextAlign(12)
    annotation.SetTextSize(0.030)
    annotation.AddText(
        f"Held-out minimum for N #geq {plot_min_n}: N={minimum_n}"
    )
    annotation.AddText(
        f"RMS change vs N={reference_n}: {relative[minimum_index]:+.3g}%"
    )
    annotation.AddText("Independent fit at every N")
    if degradation_n is not None:
        annotation.AddText(
            f"First 5-N run >5% above minimum: N={degradation_n}",
        )
    else:
        annotation.AddText("No 5-N run >5% above minimum by N=210")
    annotation.Draw()

    bottom.cd()
    subset = n_values >= plot_min_n
    relative_graph = graph(
        n_values[subset], relative[subset], ROOT.kRed + 1, 20
    )
    relative_graph.SetTitle(
        f"Held-out change relative to N={reference_n};Number of fitted terms N;Held-out RMS change (%)"
    )
    relative_graph.Draw("ALP")
    ROOT.gPad.Update()
    zero = ROOT.TLine(
        float(n_values[subset][0]),
        0.0,
        float(n_values[subset][-1]),
        0.0,
    )
    zero.SetLineStyle(3)
    zero.SetLineColor(ROOT.kGray + 1)
    zero.Draw()
    min_bottom = ROOT.TLine(
        minimum_n,
        ROOT.gPad.GetUymin(),
        minimum_n,
        ROOT.gPad.GetUymax(),
    )
    min_bottom.SetLineStyle(2)
    min_bottom.SetLineColor(ROOT.kGray + 2)
    min_bottom.Draw()
    degradation_bottom = None
    if degradation_n is not None:
        degradation_bottom = ROOT.TLine(
            degradation_n,
            ROOT.gPad.GetUymin(),
            degradation_n,
            ROOT.gPad.GetUymax(),
        )
        degradation_bottom.SetLineStyle(3)
        degradation_bottom.SetLineColor(ROOT.kRed + 1)
        degradation_bottom.Draw()

    target_directory = output / method / foil_name(foil)
    target_directory.mkdir(parents=True, exist_ok=True)
    stem = target_directory / f"{TARGET_NAMES[target_index]}_residual_ladder"
    canvas.SaveAs(str(stem.with_suffix(".png")))
    canvas.SaveAs(str(stem.with_suffix(".pdf")))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--current-coefficients", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-n", type=int, default=126)
    parser.add_argument("--file-id", type=int, default=-1)
    parser.add_argument("--nfit-max", type=int, default=200000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    (
        training_design,
        training_targets,
        _,
        terms,
        _,
        _,
        _,
    ) = load_problem(args.campaign, args.metadata, args.file_id, args.nfit_max)
    training_foils = selected_foil_labels(
        args.campaign, args.metadata, args.file_id, args.nfit_max
    )
    (
        heldout_design,
        heldout_targets,
        _,
        _,
        heldout_auxiliary,
    ) = load_heldout_problem(
        args.campaign, args.metadata, terms, args.file_id, args.nfit_max
    )
    heldout_foils = heldout_auxiliary["true_z"]
    if len(training_foils) != len(training_targets):
        raise ValueError("Training foil labels do not match the selected sample")

    foil_values = sorted(
        set(float(value) for value in np.unique(training_foils)).union(
            float(value) for value in np.unique(heldout_foils)
        )
    )
    print(f"Nominal foils: {foil_values}", flush=True)
    print("Computing independent rescaled direct-X solution at every N...", flush=True)
    direct_solutions, direct_ranks, direct_conditions, norms = rescaled_direct_ladder(
        training_design, training_targets
    )
    current_archive = np.load(args.current_coefficients)
    current_solutions = current_archive["coefficients"]
    if current_solutions.shape != direct_solutions.shape:
        raise ValueError(
            f"Current coefficient ladder shape {current_solutions.shape} does not match {direct_solutions.shape}"
        )

    print("Forming per-foil sufficient statistics...", flush=True)
    training_stats = sufficient_statistics(
        training_design, training_targets, training_foils, foil_values
    )
    heldout_stats = sufficient_statistics(
        heldout_design, heldout_targets, heldout_foils, foil_values
    )

    methods = {
        "rescaled_direct_X": direct_solutions,
        "current_XTX": current_solutions,
    }
    maximum_terms = training_design.shape[1]
    n_values = np.arange(1, maximum_terms + 1, dtype=np.float64)
    records = []
    onset_records = []
    curves = {}
    for method, solutions in methods.items():
        for foil in foil_values:
            training_curve = np.empty((maximum_terms, 3))
            heldout_curve = np.empty((maximum_terms, 3))
            for nterms in range(1, maximum_terms + 1):
                coefficients = solutions[nterms - 1, :nterms]
                training_curve[nterms - 1] = residual_rms(
                    training_stats[foil], coefficients, nterms
                )
                heldout_curve[nterms - 1] = residual_rms(
                    heldout_stats[foil], coefficients, nterms
                )
                for target, (name, units) in enumerate(
                    zip(TARGET_NAMES, TARGET_UNITS)
                ):
                    records.append(
                        {
                            "method": method,
                            "foil_z_cm": foil,
                            "N": nterms,
                            "target": name,
                            "units": units,
                            "training_count": training_stats[foil]["count"],
                            "heldout_count": heldout_stats[foil]["count"],
                            "training_residual_rms": training_curve[nterms - 1, target],
                            "heldout_residual_rms": heldout_curve[nterms - 1, target],
                            "heldout_change_from_reference_percent": math.nan,
                        }
                    )
            curves[(method, foil)] = (training_curve, heldout_curve)
            for target, (name, units) in enumerate(
                zip(TARGET_NAMES, TARGET_UNITS)
            ):
                best_n, onset_n = minimum_and_degradation_onset(
                    heldout_curve[:, target], minimum_n=35
                )
                onset_records.append(
                    {
                        "method": method,
                        "foil_z_cm": foil,
                        "target": name,
                        "units": units,
                        "training_count": training_stats[foil]["count"],
                        "heldout_count": heldout_stats[foil]["count"],
                        "best_N_ge_35": best_n,
                        "best_heldout_residual_rms": heldout_curve[
                            best_n - 1, target
                        ],
                        "first_5N_run_above_minimum_by_5pct_N": (
                            onset_n if onset_n is not None else ""
                        ),
                    }
                )

    reference_index = args.reference_n - 1
    for record in records:
        _, heldout_curve = curves[(record["method"], record["foil_z_cm"])]
        target = TARGET_NAMES.index(record["target"])
        reference = heldout_curve[reference_index, target]
        record["heldout_change_from_reference_percent"] = 100.0 * (
            record["heldout_residual_rms"] / reference - 1.0
        )

    table_path = args.output / "angular_per_foil_residual_ladder.tsv"
    with table_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(records[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(records)
    onset_path = args.output / "angular_per_foil_overfitting_onsets.tsv"
    with onset_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(onset_records[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(onset_records)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    for method in methods:
        for foil in foil_values:
            training_curve, heldout_curve = curves[(method, foil)]
            for target in range(3):
                plot_ladder(
                    args.output,
                    method,
                    foil,
                    target,
                    n_values,
                    training_curve[:, target],
                    heldout_curve[:, target],
                    training_stats[foil]["count"],
                    heldout_stats[foil]["count"],
                    args.reference_n,
                )

    summary_lines = [
        "Angular per-foil training/held-out residual ladder",
        f"training events = {len(training_targets)}",
        f"held-out eligible events = {len(heldout_targets)}",
        f"foils = {foil_values}",
        f"evaluated N = 1..{maximum_terms}",
        f"relative-change reference = N={args.reference_n}",
        "primary method = direct SVD of unit-column X, independently refit at every N",
        "comparison method = current unscaled-XTX ROOT-SVD coefficient ladder",
        "",
    ]
    for method in methods:
        summary_lines.append(method)
        for foil in foil_values:
            _, heldout_curve = curves[(method, foil)]
            parts = []
            for target, name in enumerate(TARGET_NAMES):
                best, onset = minimum_and_degradation_onset(
                    heldout_curve[:, target], minimum_n=35
                )
                parts.append(
                    f"{name}: best N>={35} is {best}, RMS={heldout_curve[best - 1, target]:.9g} {TARGET_UNITS[target]}, "
                    f"first 5-N run >5% above minimum={onset if onset is not None else 'none through 210'}"
                )
            summary_lines.append(
                f"  foil {foil:g} cm, training={training_stats[foil]['count']}, heldout={heldout_stats[foil]['count']}: "
                + "; ".join(parts)
            )
        summary_lines.append("")
    summary_lines.append(f"Table: {table_path}")
    summary_lines.append(f"Onset table: {onset_path}")
    summary_path = args.output / "angular_per_foil_residual_ladder_summary.txt"
    summary_path.write_text("\n".join(summary_lines) + "\n")
    np.savez_compressed(
        args.output / "angular_rescaled_direct_ladder_coefficients.npz",
        coefficients=direct_solutions,
        ranks=direct_ranks,
        conditions=direct_conditions,
        column_norms=norms,
        term_codes=np.asarray(["".join(map(str, term)) for term in terms]),
    )
    print(summary_path.read_text())


if __name__ == "__main__":
    main()
