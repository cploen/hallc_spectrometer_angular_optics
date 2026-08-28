#!/usr/bin/env python3
"""Compare angular-optics least-squares formulations on the same selected data.

The comparison mirrors the earlier delta study:

* the saved/current unscaled normal-equation solve, X^T X;
* a direct SVD solve of unscaled X;
* a direct SVD solve after giving every column of X unit length;
* an SVD solve of the normal equations after the same column scaling.

Scaled solutions are transformed back to the original transport-matrix
coefficient convention before predictions or coefficients are compared.
"""

from __future__ import annotations

import argparse
import csv
import math
from array import array
from pathlib import Path

import numpy as np
import ROOT

from angular_term_ladder import (
    COEFF_RE,
    evaluate_terms,
    load_problem,
    read_clean_matrix_rows,
    root_svd,
)
from preliminary_angular_conditioning import (
    build_design,
    read_campaign_rows,
    read_optics_metadata,
    selected_indices,
)


TARGET_NAMES = ("xptar", "ytar", "yptar")
TARGET_LABELS = ("x'_{tar}", "y_{tar}", "y'_{tar}")
TARGET_UNITS = ("mrad", "cm", "mrad")
PHYSICAL_SCALES = np.asarray([1000.0, 100.0, 1000.0])


def read_saved_fit(path: Path, terms):
    rows = []
    for line in path.read_text().splitlines():
        match = COEFF_RE.match(line)
        if not match:
            continue
        coefficients = tuple(float(match.group(index)) for index in range(1, 5))
        exponents = tuple(int(digit) for digit in match.group(5))
        rows.append((coefficients, exponents))
    fit_rows = [row for row in rows if row[1][4] == 0]
    fixed_rows = [row for row in rows if row[1][4] != 0]
    if len(fit_rows) != len(terms):
        raise ValueError(
            f"Expected {len(terms)} fitted rows in {path}, found {len(fit_rows)}"
        )
    for column, (row, term) in enumerate(zip(fit_rows, terms)):
        if row[1] != term:
            raise ValueError(
                f"Term mismatch at column {column}: matrix={row[1]}, design={term}"
            )
    coefficients = np.asarray([row[0][:3] for row in fit_rows], dtype=np.float64)
    return coefficients, len(rows), len(fit_rows), len(fixed_rows)


def rms_columns(values):
    return np.sqrt(np.mean(values * values, axis=0))


def solution_metrics(design, targets, coefficients):
    predictions = design @ coefficients
    residuals = predictions - targets
    return predictions, rms_columns(residuals) * PHYSICAL_SCALES


def compare_solutions(first, second, design, column_norms):
    difference = first - second
    first_raw_norm = np.linalg.norm(first, axis=0)
    second_raw_norm = np.linalg.norm(second, axis=0)
    raw_relative = np.linalg.norm(difference, axis=0) / np.maximum(
        np.maximum(first_raw_norm, second_raw_norm), np.finfo(float).tiny
    )
    scaled_difference = column_norms[:, None] * difference
    first_scaled = column_norms[:, None] * first
    second_scaled = column_norms[:, None] * second
    scaled_relative = np.linalg.norm(scaled_difference, axis=0) / np.maximum(
        np.maximum(np.linalg.norm(first_scaled, axis=0), np.linalg.norm(second_scaled, axis=0)),
        np.finfo(float).tiny,
    )
    prediction_difference = design @ difference
    prediction_rms = rms_columns(prediction_difference) * PHYSICAL_SCALES
    prediction_max = np.max(np.abs(prediction_difference), axis=0) * PHYSICAL_SCALES
    return raw_relative, scaled_relative, prediction_rms, prediction_max


def direct_prefix_spectra(design, norms):
    # For a prefix-ordered basis, X=QR gives X_N=Q_N R_N. Q_N has orthonormal
    # columns, so the singular values of each large X_N are exactly those of
    # the small leading R_N block. One QR therefore supplies the full ladder.
    r_full = np.linalg.qr(design, mode="r")
    raw = []
    scaled = []
    eps = np.finfo(np.float64).eps
    for nterms in range(1, design.shape[1] + 1):
        r_prefix = r_full[:nterms, :nterms]
        raw_singular = np.linalg.svd(r_prefix, compute_uv=False)
        scaled_singular = np.linalg.svd(
            r_prefix / norms[:nterms][None, :], compute_uv=False
        )
        raw.append(
            (
                raw_singular[0] / raw_singular[-1],
                int(np.count_nonzero(raw_singular > eps * raw_singular[0])),
            )
        )
        scaled.append(
            (
                scaled_singular[0] / scaled_singular[-1],
                int(np.count_nonzero(scaled_singular > eps * scaled_singular[0])),
            )
        )
    return raw, scaled


def eligible_indices(arrays, foils, delta_edges):
    delta = arrays["delta"]
    foil_ok = np.zeros(len(delta), dtype=bool)
    for center in foils:
        foil_ok |= np.abs(arrays["ztarT"] - center) < 2.5
    delta_bin_ok = np.zeros(len(delta), dtype=bool)
    for lower, upper in zip(delta_edges[:-1], delta_edges[1:]):
        delta_bin_ok |= (delta >= lower) & (delta < upper)
    ysieve_ok = np.zeros(len(delta), dtype=bool)
    for center in (np.arange(9) - 4) * 0.6 * 2.54:
        ysieve_ok |= np.abs(arrays["ysT"] - center) < 0.5
    return np.flatnonzero(
        (np.abs(delta) < 100.0)
        & (delta > -15.0)
        & (delta < 30.0)
        & foil_ok
        & delta_bin_ok
        & ysieve_ok
    )


def read_central_angles(path, wanted_ids):
    result = {}
    for line in path.read_text().splitlines():
        fields = [field.strip() for field in line.split(",")]
        if fields and fields[0].isdigit() and int(fields[0]) in wanted_ids:
            result[int(fields[0])] = float(fields[2])
    missing = wanted_ids - result.keys()
    if missing:
        raise ValueError(f"Missing central angles for {sorted(missing)}")
    return result


def load_heldout_problem(campaign, metadata_path, terms, file_id, nfit_max):
    settings = read_campaign_rows(next((campaign / "config").glob("rungroups_*_inputs.tsv")))
    metadata = read_optics_metadata(metadata_path, {setting[1] for setting in settings})
    angles_deg = read_central_angles(
        metadata_path, {setting[1] for setting in settings}
    )
    matrix_rows = read_clean_matrix_rows(campaign / "config/oldfit.dat")
    xtar_rows = [row for row in matrix_rows if row[1][4] != 0]
    xtar_coefficients = np.asarray([row[0][:3] for row in xtar_rows])
    branches = [
        "delta", "ztarT", "ysT", "xfp", "xpfp", "yfp", "ypfp",
        "xtar", "xptarT", "yptarT", "ytarT",
    ]
    design_parts = []
    target_parts = []
    setting_parts = []
    xtar_parts = []
    true_z_parts = []
    angle_parts = []
    ymis_parts = []
    xbeam_parts = []
    summary = []
    global_count = 0
    for setting_index, (rungroup, optics_id) in enumerate(settings):
        path = campaign / "06a_fit_ntuple/root" / f"Optics_{optics_id}_{file_id}_fit_tree_gmm.root"
        arrays = ROOT.RDataFrame("TFit", str(path)).AsNumpy(branches)
        selected, _ = selected_indices(
            arrays, *metadata[optics_id], global_count, nfit_max
        )
        eligible = eligible_indices(arrays, *metadata[optics_id])
        heldout = np.setdiff1d(eligible, selected, assume_unique=True)
        global_count += len(selected)
        if not len(heldout):
            summary.append((rungroup, 0))
            continue
        design_parts.append(build_design(arrays, heldout, terms))
        xtar_values = evaluate_terms(arrays, heldout, xtar_rows)
        xtar_prediction = xtar_values @ xtar_coefficients
        angle_deg = angles_deg[optics_id]
        angle = math.radians(angle_deg)
        ymis = 0.1 * (0.52 - 0.012 * abs(angle_deg) + 0.002 * abs(angle_deg) ** 2)
        true_z = arrays["ztarT"][heldout]
        true_yptar = arrays["yptarT"][heldout]
        true_ytar = arrays["ytarT"][heldout]
        xbeam = (
            true_ytar
            + ymis
            - true_z * (math.sin(angle) - true_yptar * math.cos(angle))
        ) / (math.cos(angle) + true_yptar * math.sin(angle))
        target_parts.append(
            np.column_stack(
                (
                    arrays["xptarT"][heldout] - xtar_prediction[:, 0],
                    (arrays["ytarT"][heldout] - 100.0 * xtar_prediction[:, 1]) / 100.0,
                    arrays["yptarT"][heldout] - xtar_prediction[:, 2],
                )
            )
        )
        setting_parts.append(np.full(len(heldout), setting_index, dtype=np.int16))
        xtar_parts.append(xtar_prediction)
        true_z_parts.append(true_z)
        angle_parts.append(np.full(len(heldout), angle))
        ymis_parts.append(np.full(len(heldout), ymis))
        xbeam_parts.append(xbeam)
        summary.append((rungroup, len(heldout)))
    return (
        np.vstack(design_parts),
        np.vstack(target_parts),
        np.concatenate(setting_parts),
        summary,
        {
            "xtar_prediction": np.vstack(xtar_parts),
            "true_z": np.concatenate(true_z_parts),
            "angle": np.concatenate(angle_parts),
            "ymis": np.concatenate(ymis_parts),
            "xbeam": np.concatenate(xbeam_parts),
        },
    )


def reconstruct_vertex(clean_prediction, auxiliary):
    total_ytar_cm = 100.0 * (
        clean_prediction[:, 1] + auxiliary["xtar_prediction"][:, 1]
    )
    total_yptar = clean_prediction[:, 2] + auxiliary["xtar_prediction"][:, 2]
    sine = np.sin(auxiliary["angle"])
    cosine = np.cos(auxiliary["angle"])
    return (
        total_ytar_cm
        + auxiliary["ymis"]
        - auxiliary["xbeam"] * (cosine + total_yptar * sine)
    ) / (sine - total_yptar * cosine)


def load_gram_ladder(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def graph(x, y, color, marker):
    result = ROOT.TGraph(len(x), array("d", x), array("d", y))
    result.SetLineColor(color)
    result.SetMarkerColor(color)
    result.SetMarkerStyle(marker)
    result.SetMarkerSize(0.45)
    result.SetLineWidth(2)
    return result


def save_canvas(canvas, path):
    canvas.SaveAs(str(path.with_suffix(".png")))
    canvas.SaveAs(str(path.with_suffix(".pdf")))


def condition_plot(output, direct_raw, direct_scaled, gram_rows):
    n_values = [float(index) for index in range(1, len(direct_raw) + 1)]
    curves = (
        ("Current: unscaled X^{T}X", [float(row["kappa_Ay"]) for row in gram_rows], ROOT.kBlue + 1, 20),
        ("Direct SVD of unscaled X", [row[0] for row in direct_raw], ROOT.kGreen + 2, 21),
        ("Direct SVD of rescaled X", [row[0] for row in direct_scaled], ROOT.kRed + 1, 22),
        ("SVD of rescaled X^{T}X", [float(row["kappa_scaled_Gram"]) for row in gram_rows], ROOT.kOrange + 7, 23),
    )
    canvas = ROOT.TCanvas("c_solver_condition", "Condition by solve method", 1120, 800)
    canvas.SetLogy()
    multi = ROOT.TMultiGraph()
    graphs = []
    for _, values, color, marker in curves:
        item = graph(n_values, values, color, marker)
        multi.Add(item, "LP")
        graphs.append(item)
    multi.SetTitle(
        "Angular fit conditioning by numerical method;Number of adjusted coefficients N;Condition number (smaller is better)"
    )
    multi.Draw("A")
    legend = ROOT.TLegend(0.13, 0.68, 0.53, 0.86)
    for (label, _, _, _), item in zip(curves, graphs):
        legend.AddEntry(item, label, "lp")
    legend.Draw()
    save_canvas(canvas, output / "condition_by_solve_method_vs_N")


def retained_plot(output, direct_raw, direct_scaled, gram_rows):
    n_values = [float(index) for index in range(1, len(direct_raw) + 1)]
    dropped_current = [float(row["truncated"]) for row in gram_rows]
    dropped_direct = [float(index - item[1]) for index, item in enumerate(direct_raw, 1)]
    dropped_scaled_direct = [
        float(index - item[1]) for index, item in enumerate(direct_scaled, 1)
    ]
    dropped_scaled_gram = [float(row["truncated_scaled"]) for row in gram_rows]
    curves = (
        ("Current: unscaled X^{T}X", dropped_current, ROOT.kBlue + 1, 20),
        ("Direct SVD of unscaled X", dropped_direct, ROOT.kGreen + 2, 21),
        ("Direct SVD of rescaled X", dropped_scaled_direct, ROOT.kRed + 1, 22),
        ("SVD of rescaled X^{T}X", dropped_scaled_gram, ROOT.kOrange + 7, 23),
    )
    canvas = ROOT.TCanvas("c_solver_dropped", "Discarded directions", 1120, 800)
    multi = ROOT.TMultiGraph()
    graphs = []
    for _, values, color, marker in curves:
        item = graph(n_values, values, color, marker)
        multi.Add(item, "LP")
        graphs.append(item)
    multi.SetTitle(
        "Directions suppressed by each solver's numerical cutoff;Number of adjusted coefficients N;Number of suppressed coefficient combinations"
    )
    multi.Draw("A")
    multi.SetMinimum(0.0)
    legend = ROOT.TLegend(0.13, 0.68, 0.53, 0.86)
    for (label, _, _, _), item in zip(curves, graphs):
        legend.AddEntry(item, label, "lp")
    legend.Draw()
    save_canvas(canvas, output / "discarded_directions_by_solve_method_vs_N")


def prediction_plot(output, comparisons, suffix="", sample_label="calibration sample"):
    methods = list(comparisons)
    canvas = ROOT.TCanvas("c_prediction_differences", "Prediction differences", 1120, 900)
    canvas.Divide(1, 3)
    histograms = []
    colors = (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7)
    for target, (label, units) in enumerate(zip(TARGET_LABELS, TARGET_UNITS), 1):
        canvas.cd(target)
        ROOT.gPad.SetLogy()
        histogram = ROOT.TH1D(
            f"h_prediction_difference_{target}",
            f"Difference from rescaled direct-X solution on {sample_label}: {label};Solve method;RMS difference ({units})",
            len(methods),
            0.0,
            float(len(methods)),
        )
        histogram.SetStats(False)
        histogram.SetFillColor(colors[target - 1])
        histogram.SetLineColor(colors[target - 1])
        for index, method in enumerate(methods, 1):
            value = max(float(comparisons[method][2][target - 1]), 1e-15)
            histogram.SetBinContent(index, value)
            histogram.GetXaxis().SetBinLabel(index, method)
        histogram.LabelsOption("v", "X")
        histogram.SetMinimum(1e-6 if units == "mrad" else 1e-7)
        histogram.Draw("HIST")
        histograms.append(histogram)
    save_canvas(canvas, output / f"prediction_difference_from_rescaled_direct_X{suffix}")


def candidate_residual_plot(output, records, sample):
    n_values = sorted({int(record["N"]) for record in records})
    canvas = ROOT.TCanvas(
        f"c_candidate_residual_{sample}",
        f"Candidate residuals on {sample}",
        1120,
        900,
    )
    canvas.Divide(1, 3)
    graphs = []
    for target, (label, units) in enumerate(zip(TARGET_LABELS, TARGET_UNITS), 1):
        canvas.cd(target)
        multi = ROOT.TMultiGraph()
        current_values = [
            next(
                record[f"{sample}_{TARGET_NAMES[target - 1]}_rms"]
                for record in records
                if record["N"] == nterms and record["method"] == "current XTX"
            )
            for nterms in n_values
        ]
        scaled_values = [
            next(
                record[f"{sample}_{TARGET_NAMES[target - 1]}_rms"]
                for record in records
                if record["N"] == nterms and record["method"] == "rescaled direct X"
            )
            for nterms in n_values
        ]
        current_graph = graph(
            [float(value) for value in n_values], current_values, ROOT.kBlue + 1, 20
        )
        scaled_graph = graph(
            [float(value) for value in n_values], scaled_values, ROOT.kRed + 1, 22
        )
        multi.Add(current_graph, "LP")
        multi.Add(scaled_graph, "LP")
        multi.SetTitle(
            f"Candidate complexity on {sample.replace('_', ' ')}: {label};N;Residual RMS ({units})"
        )
        multi.Draw("A")
        legend = ROOT.TLegend(0.13, 0.72, 0.43, 0.86)
        legend.AddEntry(current_graph, "Current unscaled X^{T}X", "lp")
        legend.AddEntry(scaled_graph, "Rescaled direct X", "lp")
        legend.Draw()
        graphs.extend((current_graph, scaled_graph, multi, legend))
    save_canvas(canvas, output / f"candidate_residuals_{sample}")


def candidate_training_zoom_plot(output, records):
    """Show the small training-RMS differences without the low-N scale."""
    n_values = sorted({int(record["N"]) for record in records if int(record["N"]) >= 35})
    canvas = ROOT.TCanvas(
        "c_candidate_training_zoom",
        "Training residuals zoom",
        1200,
        900,
    )
    canvas.Divide(1, 3)
    objects = []
    for target, (name, label, units) in enumerate(
        zip(TARGET_NAMES, TARGET_LABELS, TARGET_UNITS), 1
    ):
        canvas.cd(target)
        ROOT.gPad.SetRightMargin(0.30)
        current_values = [
            next(
                record[f"fit_{name}_rms"]
                for record in records
                if record["N"] == nterms and record["method"] == "current XTX"
            )
            for nterms in n_values
        ]
        scaled_values = [
            next(
                record[f"fit_{name}_rms"]
                for record in records
                if record["N"] == nterms
                and record["method"] == "rescaled direct X"
            )
            for nterms in n_values
        ]
        current_graph = graph(
            [float(value) for value in n_values], current_values, ROOT.kBlue + 1, 20
        )
        scaled_graph = graph(
            [float(value) for value in n_values], scaled_values, ROOT.kRed + 1, 22
        )
        current_graph.SetMarkerSize(0.8)
        scaled_graph.SetMarkerSize(0.8)
        multi = ROOT.TMultiGraph()
        multi.Add(current_graph, "LP")
        multi.Add(scaled_graph, "LP")
        all_values = current_values + scaled_values
        span = max(all_values) - min(all_values)
        padding = 0.08 * span
        multi.SetMinimum(min(all_values) - padding)
        multi.SetMaximum(max(all_values) + padding)
        multi.SetTitle(
            f"Training residuals, N #geq 35: {label};N;Residual RMS ({units})"
        )
        multi.Draw("A")
        objects.extend((multi, current_graph, scaled_graph))
        if target == 1:
            legend = ROOT.TLegend(0.72, 0.70, 0.99, 0.86)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.AddEntry(current_graph, "Current: unscaled X^{T}X", "lp")
            legend.AddEntry(scaled_graph, "Direct: scaled X", "lp")
            legend.Draw()
            objects.append(legend)
    save_canvas(canvas, output / "candidate_residuals_fit_zoom")


def candidate_heldout_zoom_plot(output, records):
    """Show held-out minima before the large high-N rise dominates the scale."""
    n_values = sorted(
        {
            int(record["N"])
            for record in records
            if 35 <= int(record["N"]) <= 180
        }
    )
    canvas = ROOT.TCanvas(
        "c_candidate_heldout_zoom",
        "Held-out residuals zoom",
        1200,
        900,
    )
    canvas.Divide(1, 3)
    objects = []
    for target, (name, label, units) in enumerate(
        zip(TARGET_NAMES, TARGET_LABELS, TARGET_UNITS), 1
    ):
        canvas.cd(target)
        ROOT.gPad.SetRightMargin(0.30)
        current_values = [
            next(
                record[f"heldout_{name}_rms"]
                for record in records
                if record["N"] == nterms and record["method"] == "current XTX"
            )
            for nterms in n_values
        ]
        scaled_values = [
            next(
                record[f"heldout_{name}_rms"]
                for record in records
                if record["N"] == nterms
                and record["method"] == "rescaled direct X"
            )
            for nterms in n_values
        ]
        current_graph = graph(
            [float(value) for value in n_values], current_values, ROOT.kBlue + 1, 20
        )
        scaled_graph = graph(
            [float(value) for value in n_values], scaled_values, ROOT.kRed + 1, 22
        )
        current_graph.SetMarkerSize(0.8)
        scaled_graph.SetMarkerSize(0.8)
        multi = ROOT.TMultiGraph()
        multi.Add(current_graph, "LP")
        multi.Add(scaled_graph, "LP")
        all_values = current_values + scaled_values
        span = max(all_values) - min(all_values)
        padding = 0.08 * span
        multi.SetMinimum(min(all_values) - padding)
        multi.SetMaximum(max(all_values) + padding)
        multi.SetTitle(
            f"Held-out residuals, 35 #leq N #leq 180: {label};N;Residual RMS ({units})"
        )
        multi.Draw("A")

        minimum_index = int(np.argmin(scaled_values))
        best_n = n_values[minimum_index]
        scaled_minimum = scaled_values[minimum_index]
        current_at_best_n = current_values[minimum_index]
        percent = 100.0 * (scaled_minimum - current_at_best_n) / current_at_best_n
        direction = "above" if percent >= 0.0 else "below"
        note = ROOT.TLatex()
        note.SetNDC(True)
        note.SetTextSize(0.055)
        note.DrawLatex(0.72, 0.52, f"Scaled minimum: N={best_n}")
        note.DrawLatex(
            0.72,
            0.43,
            f"{abs(percent):.1f}% {direction} current",
        )
        objects.extend((multi, current_graph, scaled_graph, note))
        if target == 1:
            legend = ROOT.TLegend(0.72, 0.70, 0.99, 0.86)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.AddEntry(current_graph, "Current: unscaled X^{T}X", "lp")
            legend.AddEntry(scaled_graph, "Direct: scaled X", "lp")
            legend.Draw()
            objects.append(legend)
    save_canvas(canvas, output / "candidate_residuals_heldout_zoom")


def candidate_vertex_plot(output, records):
    n_values = sorted({int(record["N"]) for record in records})
    current_values = [
        next(
            record["heldout_vertex_rms_cm"]
            for record in records
            if record["N"] == nterms and record["method"] == "current XTX"
        )
        for nterms in n_values
    ]
    scaled_values = [
        next(
            record["heldout_vertex_rms_cm"]
            for record in records
            if record["N"] == nterms and record["method"] == "rescaled direct X"
        )
        for nterms in n_values
    ]
    canvas = ROOT.TCanvas("c_candidate_vertex", "Candidate vertex residual", 1120, 800)
    canvas.SetLogy()
    current_graph = graph(
        [float(value) for value in n_values], current_values, ROOT.kBlue + 1, 20
    )
    scaled_graph = graph(
        [float(value) for value in n_values], scaled_values, ROOT.kRed + 1, 22
    )
    multi = ROOT.TMultiGraph()
    multi.Add(current_graph, "LP")
    multi.Add(scaled_graph, "LP")
    multi.SetTitle(
        "Reconstructed HMS vertex on eligible events excluded from the fit;N;Vertex z residual RMS (cm)"
    )
    multi.Draw("A")
    legend = ROOT.TLegend(0.13, 0.72, 0.43, 0.86)
    legend.AddEntry(current_graph, "Current unscaled X^{T}X", "lp")
    legend.AddEntry(scaled_graph, "Rescaled direct X", "lp")
    legend.Draw()
    save_canvas(canvas, output / "candidate_vertex_residuals_heldout")


def training_heldout_plot(output, records):
    n_values = sorted({int(record["N"]) for record in records})
    canvas = ROOT.TCanvas(
        "c_training_heldout", "Training and held-out residuals", 1120, 900
    )
    canvas.Divide(1, 3)
    objects = []
    for target, (name, label, units) in enumerate(
        zip(TARGET_NAMES, TARGET_LABELS, TARGET_UNITS), 1
    ):
        canvas.cd(target)
        multi = ROOT.TMultiGraph()
        specifications = (
            ("current XTX", "fit", ROOT.kBlue + 1, 20, 2),
            ("current XTX", "heldout", ROOT.kBlue + 1, 24, 1),
            ("rescaled direct X", "fit", ROOT.kRed + 1, 22, 2),
            ("rescaled direct X", "heldout", ROOT.kRed + 1, 26, 1),
        )
        graphs = []
        for method, sample, color, marker, line_style in specifications:
            values = [
                next(
                    record[f"{sample}_{name}_rms"]
                    for record in records
                    if record["N"] == nterms and record["method"] == method
                )
                for nterms in n_values
            ]
            item = graph(
                [float(value) for value in n_values], values, color, marker
            )
            item.SetLineStyle(line_style)
            multi.Add(item, "LP")
            graphs.append(item)
        multi.SetTitle(
            f"Training improvement versus held-out behavior: {label};N;Residual RMS ({units})"
        )
        multi.Draw("A")
        legend = ROOT.TLegend(0.12, 0.62, 0.48, 0.86)
        legend.AddEntry(graphs[0], "Current X^{T}X: training", "lp")
        legend.AddEntry(graphs[1], "Current X^{T}X: held out", "lp")
        legend.AddEntry(graphs[2], "Rescaled direct X: training", "lp")
        legend.AddEntry(graphs[3], "Rescaled direct X: held out", "lp")
        legend.Draw()
        objects.extend((multi, legend, *graphs))
    save_canvas(canvas, output / "training_vs_heldout_residuals_by_N")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--saved-matrix", type=Path, required=True)
    parser.add_argument("--gram-ladder", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--file-id", type=int, default=-1)
    parser.add_argument("--nfit-max", type=int, default=200000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    design, targets, _, terms, selected_summary, valid_seed_rows, fixed_seed_rows = load_problem(
        args.campaign, args.metadata, args.file_id, args.nfit_max
    )
    events, columns = design.shape
    saved, saved_rows, saved_fit_rows, saved_fixed_rows = read_saved_fit(
        args.saved_matrix, terms
    )
    norms = np.linalg.norm(design, axis=0)
    scaled_design = design / norms
    eps = np.finfo(np.float64).eps

    print("Solving unscaled X directly with SVD...", flush=True)
    direct, _, direct_rank, direct_singular = np.linalg.lstsq(
        design, targets, rcond=eps
    )
    print("Solving unit-column X directly with SVD...", flush=True)
    scaled_coordinates, _, scaled_direct_rank, scaled_direct_singular = np.linalg.lstsq(
        scaled_design, targets, rcond=eps
    )
    scaled_direct = scaled_coordinates / norms[:, None]

    print("Solving unit-column normal equations with ROOT SVD...", flush=True)
    scaled_gram_solution, scaled_gram_singular, scaled_gram_threshold, scaled_gram_rank, _ = root_svd(
        scaled_design.T @ scaled_design, scaled_design.T @ targets
    )
    scaled_gram = scaled_gram_solution / norms[:, None]

    print("Reforming the current unscaled normal equations for comparison...", flush=True)
    current_reformed, current_singular, current_threshold, current_rank, _ = root_svd(
        design.T @ design, design.T @ targets
    )

    solutions = {
        "saved current XTX": saved,
        "reformed current XTX": current_reformed,
        "unscaled direct X": direct,
        "rescaled XTX": scaled_gram,
        "rescaled direct X": scaled_direct,
    }
    metrics = {}
    for name, coefficients in solutions.items():
        prediction, residual = solution_metrics(design, targets, coefficients)
        metrics[name] = {"prediction": prediction, "residual": residual}

    reference = solutions["rescaled direct X"]
    comparisons = {}
    for name, coefficients in solutions.items():
        if name == "rescaled direct X":
            continue
        comparisons[name] = compare_solutions(
            coefficients, reference, design, norms
        )

    print("Computing direct-X condition ladders from one QR factorization...", flush=True)
    direct_raw_ladder, direct_scaled_ladder = direct_prefix_spectra(design, norms)
    gram_rows = load_gram_ladder(args.gram_ladder)
    if len(gram_rows) != columns:
        raise ValueError(f"Expected {columns} Gram ladder rows, found {len(gram_rows)}")

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    condition_plot(args.output, direct_raw_ladder, direct_scaled_ladder, gram_rows)
    retained_plot(args.output, direct_raw_ladder, direct_scaled_ladder, gram_rows)
    prediction_plot(args.output, comparisons)

    print("Loading eligible events excluded by the fit-sample caps...", flush=True)
    (
        heldout_design,
        heldout_targets,
        heldout_setting,
        heldout_summary,
        heldout_auxiliary,
    ) = load_heldout_problem(args.campaign, args.metadata, terms, args.file_id, args.nfit_max)
    heldout_metrics = {}
    heldout_comparisons = {}
    for name, coefficients in solutions.items():
        prediction, residual = solution_metrics(
            heldout_design, heldout_targets, coefficients
        )
        heldout_metrics[name] = {"prediction": prediction, "residual": residual}
        heldout_metrics[name]["vertex"] = reconstruct_vertex(
            prediction, heldout_auxiliary
        )
        if name != "rescaled direct X":
            heldout_comparisons[name] = compare_solutions(
                coefficients, reference, heldout_design, norms
            )
    prediction_plot(
        args.output,
        heldout_comparisons,
        suffix="_heldout",
        sample_label="eligible events excluded from the fit",
    )

    candidate_ns = tuple(
        sorted(
            set(
                (5, 15, 35, 39, 40, 53, 54, 70, 126, 210)
                + tuple(range(50, 211, 10))
            )
        )
    )
    ladder_coefficients = np.load(
        args.gram_ladder.parent / "angular_term_ladder_coefficients.npz"
    )["coefficients"]
    candidate_records = []
    for nterms in candidate_ns:
        scaled_prefix = scaled_design[:, :nterms]
        scaled_c, _, scaled_rank, _ = np.linalg.lstsq(
            scaled_prefix, targets, rcond=eps
        )
        scaled_d = scaled_c / norms[:nterms, None]
        current_d = ladder_coefficients[nterms - 1, :nterms]
        for method, coefficients in (
            ("current XTX", current_d),
            ("rescaled direct X", scaled_d),
        ):
            _, fit_rms = solution_metrics(
                design[:, :nterms], targets, coefficients
            )
            heldout_prediction, heldout_rms = solution_metrics(
                heldout_design[:, :nterms], heldout_targets, coefficients
            )
            heldout_vertex = reconstruct_vertex(
                heldout_prediction, heldout_auxiliary
            )
            row = next(item for item in gram_rows if int(item["N"]) == nterms)
            candidate_records.append(
                {
                    "N": nterms,
                    "method": method,
                    "rank": int(row["rank"]) if method == "current XTX" else scaled_rank,
                    "condition": float(row["kappa_Ay"])
                    if method == "current XTX"
                    else direct_scaled_ladder[nterms - 1][0],
                    "fit_xptar_rms": fit_rms[0],
                    "fit_ytar_rms": fit_rms[1],
                    "fit_yptar_rms": fit_rms[2],
                    "heldout_xptar_rms": heldout_rms[0],
                    "heldout_ytar_rms": heldout_rms[1],
                    "heldout_yptar_rms": heldout_rms[2],
                    "heldout_vertex_rms_cm": float(
                        np.sqrt(
                            np.mean(
                                (heldout_vertex - heldout_auxiliary["true_z"]) ** 2
                            )
                        )
                    ),
                }
            )
    candidate_table_path = args.output / "angular_candidate_N_comparison.tsv"
    with candidate_table_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(candidate_records[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(candidate_records)
    candidate_residual_plot(args.output, candidate_records, "fit")
    candidate_training_zoom_plot(args.output, candidate_records)
    candidate_residual_plot(args.output, candidate_records, "heldout")
    candidate_heldout_zoom_plot(args.output, candidate_records)
    candidate_vertex_plot(args.output, candidate_records)
    training_heldout_plot(args.output, candidate_records)

    table_path = args.output / "angular_solver_comparison.tsv"
    with table_path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(
            [
                "method",
                "target",
                "residual_rms",
                "units",
                "raw_coefficient_relative_difference_from_rescaled_direct",
                "scaled_coordinate_relative_difference_from_rescaled_direct",
                "prediction_rms_difference_from_rescaled_direct",
                "prediction_max_difference_from_rescaled_direct",
            ]
        )
        for name in solutions:
            comparison = comparisons.get(name)
            for target, (target_name, units) in enumerate(zip(TARGET_NAMES, TARGET_UNITS)):
                writer.writerow(
                    [
                        name,
                        target_name,
                        metrics[name]["residual"][target],
                        units,
                        comparison[0][target] if comparison else 0.0,
                        comparison[1][target] if comparison else 0.0,
                        comparison[2][target] if comparison else 0.0,
                        comparison[3][target] if comparison else 0.0,
                    ]
                )

    heldout_table_path = args.output / "angular_solver_comparison_heldout.tsv"
    with heldout_table_path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(
            [
                "method", "target", "heldout_residual_rms", "units",
                "prediction_rms_difference_from_rescaled_direct",
                "prediction_max_difference_from_rescaled_direct",
            ]
        )
        for name in solutions:
            comparison = heldout_comparisons.get(name)
            for target, (target_name, units) in enumerate(zip(TARGET_NAMES, TARGET_UNITS)):
                writer.writerow(
                    [
                        name,
                        target_name,
                        heldout_metrics[name]["residual"][target],
                        units,
                        comparison[2][target] if comparison else 0.0,
                        comparison[3][target] if comparison else 0.0,
                    ]
                )

    full_current = gram_rows[-1]
    report_lines = [
        "Preliminary angular-optics solver comparison",
        f"events = {events}",
        f"adjusted coefficients per reconstructed quantity = {columns}",
        f"valid seed-matrix rows = {valid_seed_rows}",
        f"seed rows carried fixed because they depend on xtar = {fixed_seed_rows}",
        f"saved clean output rows = {saved_rows} ({saved_fit_rows} adjusted, {saved_fixed_rows} fixed)",
        "scaling = divide each design-matrix column by its Euclidean length, then transform coefficients back",
        f"relative singular-value cutoff used for comparisons = {eps:.17g}",
        "",
        "Full N=210 condition numbers (smaller is better):",
        f"  current unscaled XTX = {float(full_current['kappa_Ay']):.9g}",
        f"  unscaled direct X = {direct_singular[0] / direct_singular[-1]:.9g}",
        f"  rescaled direct X = {scaled_direct_singular[0] / scaled_direct_singular[-1]:.9g}",
        f"  rescaled XTX = {scaled_gram_singular[0] / scaled_gram_singular[-1]:.9g}",
        "",
        "Full N=210 retained directions using the same relative cutoff:",
        f"  current unscaled XTX = {full_current['rank']}/{columns}",
        f"  unscaled direct X = {direct_rank}/{columns}",
        f"  rescaled direct X = {scaled_direct_rank}/{columns}",
        f"  rescaled XTX = {scaled_gram_rank}/{columns}",
        "",
        "Residual RMS on the selected fit sample:",
    ]
    for name in solutions:
        values = metrics[name]["residual"]
        report_lines.append(
            f"  {name}: xptar={values[0]:.9g} mrad, ytar={values[1]:.9g} cm, yptar={values[2]:.9g} mrad"
        )
    report_lines.extend(["", "Difference from the rescaled direct-X solution:"])
    for name, comparison in comparisons.items():
        report_lines.append(
            f"  {name}: prediction RMS xptar={comparison[2][0]:.9g} mrad, "
            f"ytar={comparison[2][1]:.9g} cm, yptar={comparison[2][2]:.9g} mrad"
        )
        report_lines.append(
            f"    raw coefficient relative differences = {comparison[0][0]:.9g}, "
            f"{comparison[0][1]:.9g}, {comparison[0][2]:.9g}"
        )
        report_lines.append(
            f"    scaled-coordinate relative differences = {comparison[1][0]:.9g}, "
            f"{comparison[1][1]:.9g}, {comparison[1][2]:.9g}"
        )
    report_lines.extend(
        [
            "",
            f"Eligible events excluded by the fit-sample caps = {len(heldout_targets)}",
            "Held-out events by setting:",
        ]
    )
    for setting, count in heldout_summary:
        report_lines.append(f"  {setting}: {count}")
    report_lines.extend(["", "Residual RMS on eligible events excluded from the fit:"])
    for name in solutions:
        values = heldout_metrics[name]["residual"]
        report_lines.append(
            f"  {name}: xptar={values[0]:.9g} mrad, ytar={values[1]:.9g} cm, yptar={values[2]:.9g} mrad"
        )
    report_lines.extend(
        ["", "Held-out prediction difference from the rescaled direct-X solution:"]
    )
    for name, comparison in heldout_comparisons.items():
        report_lines.append(
            f"  {name}: RMS xptar={comparison[2][0]:.9g} mrad, "
            f"ytar={comparison[2][1]:.9g} cm, yptar={comparison[2][2]:.9g} mrad; "
            f"max xptar={comparison[3][0]:.9g} mrad, ytar={comparison[3][1]:.9g} cm, "
            f"yptar={comparison[3][2]:.9g} mrad"
        )
    report_lines.extend(["", "Reconstructed vertex z on held-out events:"])
    for name in solutions:
        vertex_residual = (
            heldout_metrics[name]["vertex"] - heldout_auxiliary["true_z"]
        )
        finite = np.isfinite(vertex_residual)
        absolute = np.abs(vertex_residual[finite])
        report_lines.append(
            f"  {name}: RMS={np.sqrt(np.mean(vertex_residual[finite] ** 2)):.9g} cm, "
            f"median abs={np.median(absolute):.9g} cm, p99 abs={np.percentile(absolute, 99.0):.9g} cm, "
            f"max abs={np.max(absolute):.9g} cm"
        )
    vertex_difference = (
        heldout_metrics["saved current XTX"]["vertex"]
        - heldout_metrics["rescaled direct X"]["vertex"]
    )
    finite_vertex_difference = vertex_difference[np.isfinite(vertex_difference)]
    absolute_vertex_difference = np.abs(finite_vertex_difference)
    report_lines.append(
        "  current-versus-rescaled vertex difference: "
        f"RMS={np.sqrt(np.mean(finite_vertex_difference ** 2)):.9g} cm, "
        f"median abs={np.median(absolute_vertex_difference):.9g} cm, "
        f"p90 abs={np.percentile(absolute_vertex_difference, 90.0):.9g} cm, "
        f"p99 abs={np.percentile(absolute_vertex_difference, 99.0):.9g} cm, "
        f"max abs={np.max(absolute_vertex_difference):.9g} cm"
    )
    report_lines.extend(
        [
            "",
            "Absolute current-matrix versus rescaled-direct prediction difference percentiles on held-out events:",
        ]
    )
    absolute_difference = np.abs(
        heldout_metrics["saved current XTX"]["prediction"]
        - heldout_metrics["rescaled direct X"]["prediction"]
    ) * PHYSICAL_SCALES
    for target, (name, units) in enumerate(zip(TARGET_NAMES, TARGET_UNITS)):
        percentiles = np.percentile(
            absolute_difference[:, target], [50.0, 90.0, 95.0, 99.0, 99.9, 100.0]
        )
        report_lines.append(
            f"  {name} ({units}): median={percentiles[0]:.9g}, p90={percentiles[1]:.9g}, "
            f"p95={percentiles[2]:.9g}, p99={percentiles[3]:.9g}, "
            f"p99.9={percentiles[4]:.9g}, max={percentiles[5]:.9g}"
        )
    report_lines.extend(
        [
            "",
            "Per-setting held-out RMS: saved current XTX versus rescaled direct X, followed by their prediction difference",
        ]
    )
    for setting_index, (setting, count) in enumerate(heldout_summary):
        if not count:
            continue
        mask = heldout_setting == setting_index
        current_residual = rms_columns(
            heldout_metrics["saved current XTX"]["prediction"][mask]
            - heldout_targets[mask]
        ) * PHYSICAL_SCALES
        scaled_residual = rms_columns(
            heldout_metrics["rescaled direct X"]["prediction"][mask]
            - heldout_targets[mask]
        ) * PHYSICAL_SCALES
        difference = rms_columns(
            heldout_metrics["saved current XTX"]["prediction"][mask]
            - heldout_metrics["rescaled direct X"]["prediction"][mask]
        ) * PHYSICAL_SCALES
        report_lines.append(
            f"  {setting} (n={count}): current=({current_residual[0]:.6g} mrad, {current_residual[1]:.6g} cm, {current_residual[2]:.6g} mrad); "
            f"rescaled=({scaled_residual[0]:.6g} mrad, {scaled_residual[1]:.6g} cm, {scaled_residual[2]:.6g} mrad); "
            f"difference=({difference[0]:.6g} mrad, {difference[1]:.6g} cm, {difference[2]:.6g} mrad)"
        )
    report_lines.extend(["", "Candidate-N fit/held-out residual summary:"])
    for nterms in candidate_ns:
        current_row = next(
            row for row in candidate_records if row["N"] == nterms and row["method"] == "current XTX"
        )
        scaled_row = next(
            row for row in candidate_records if row["N"] == nterms and row["method"] == "rescaled direct X"
        )
        report_lines.append(
            f"  N={nterms}: current rank={current_row['rank']}/{nterms}, scaled rank={scaled_row['rank']}/{nterms}; "
            f"held-out current=({current_row['heldout_xptar_rms']:.6g} mrad, {current_row['heldout_ytar_rms']:.6g} cm, {current_row['heldout_yptar_rms']:.6g} mrad); "
            f"scaled=({scaled_row['heldout_xptar_rms']:.6g} mrad, {scaled_row['heldout_ytar_rms']:.6g} cm, {scaled_row['heldout_yptar_rms']:.6g} mrad); "
            f"vertex current={current_row['heldout_vertex_rms_cm']:.6g} cm, scaled={scaled_row['heldout_vertex_rms_cm']:.6g} cm"
        )
    report_lines.extend(
        [
            "",
            f"Current reformed XTX weakest singular value / cutoff = {current_singular[-1] / current_threshold:.9g}",
            f"Rescaled XTX weakest singular value / cutoff = {scaled_gram_singular[-1] / scaled_gram_threshold:.9g}",
            "",
            f"Table: {table_path}",
            f"Held-out table: {heldout_table_path}",
            f"Candidate-N table: {candidate_table_path}",
        ]
    )
    report_path = args.output / "angular_solver_comparison.txt"
    report_path.write_text("\n".join(report_lines) + "\n")
    np.savez_compressed(
        args.output / "angular_solver_coefficients.npz",
        terms=np.asarray(["".join(map(str, term)) for term in terms]),
        column_norms=norms,
        **{name.replace(" ", "_"): value for name, value in solutions.items()},
    )
    print(report_path.read_text())


if __name__ == "__main__":
    main()
