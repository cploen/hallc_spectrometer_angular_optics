#!/usr/bin/env python3
"""Nested prefix ladder for the corrected HMS angular optics basis.

The production solve factorizes G_N = X_N^T X_N with ROOT's default relative
SVD tolerance.  Since the absolute cutoff changes with sigma_max(G_N), this
script reports and plots the solve-specific signed margin

    log10(sigma_min(G_N)) - log10(tolerance * sigma_max(G_N)).

Zero therefore means "at this solve's own truncation threshold".  A negative
value means at least the weakest mode is below it.  No common absolute
threshold is assumed across ladder rungs.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from array import array
from pathlib import Path

import numpy as np
import ROOT

from preliminary_angular_conditioning import (
    build_design,
    read_campaign_rows,
    read_optics_metadata,
    selected_indices,
)


COEFF_RE = re.compile(
    r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"(\d{5})\s*$"
)


def read_clean_matrix_rows(path: Path):
    rows = []
    in_data = False
    for line in path.read_text().splitlines():
        if not in_data:
            if line.startswith(" --------------------------------"):
                in_data = True
            continue
        if line.startswith(" --------------------------------"):
            break
        match = COEFF_RE.match(line)
        if not match:
            continue
        coefficients = tuple(float(match.group(i)) for i in range(1, 5))
        exponents = tuple(int(digit) for digit in match.group(5))
        rows.append((coefficients, exponents))
    return rows


def evaluate_terms(arrays, indices, rows):
    variables = (
        arrays["xfp"][indices] / 100.0,
        arrays["xpfp"][indices],
        arrays["yfp"][indices] / 100.0,
        arrays["ypfp"][indices],
        arrays["xtar"][indices] / 100.0,
    )
    values = np.empty((len(indices), len(rows)), dtype=np.float64)
    for column, (_, exponents) in enumerate(rows):
        term = np.ones(len(indices), dtype=np.float64)
        for variable, exponent in zip(variables, exponents):
            if exponent:
                term *= variable**exponent
        values[:, column] = term
    return values


def load_problem(campaign: Path, metadata_path: Path, file_id: int, nfit_max: int):
    config = campaign / "config"
    rungroup_files = list(config.glob("rungroups_*_inputs.tsv"))
    if len(rungroup_files) != 1:
        raise ValueError(f"Expected one rungroup TSV, found {rungroup_files}")
    settings = read_campaign_rows(rungroup_files[0])
    metadata = read_optics_metadata(metadata_path, {setting[1] for setting in settings})

    matrix_rows = read_clean_matrix_rows(config / "oldfit.dat")
    fit_rows = [row for row in matrix_rows if row[1][4] == 0]
    xtar_rows = [row for row in matrix_rows if row[1][4] != 0]
    fit_terms = [(0, 0, 0, 0, 0)] + [row[1] for row in fit_rows]
    if len(fit_terms) != 210:
        raise ValueError(f"Expected 210 corrected fit terms, found {len(fit_terms)}")

    branches = [
        "delta", "ztarT", "ysT", "xfp", "xpfp", "yfp", "ypfp",
        "xtar", "xptarT", "yptarT", "ytarT",
    ]
    design_parts = []
    target_parts = []
    setting_parts = []
    selected_summary = []
    global_count = 0

    xtar_coefficients = np.asarray([row[0][:3] for row in xtar_rows])
    for setting_index, (rungroup, optics_id) in enumerate(settings):
        tree_path = (
            campaign / "06a_fit_ntuple/root"
            / f"Optics_{optics_id}_{file_id}_fit_tree_gmm.root"
        )
        arrays = ROOT.RDataFrame("TFit", str(tree_path)).AsNumpy(branches)
        indices, per_foil = selected_indices(
            arrays, *metadata[optics_id], global_count, nfit_max
        )
        design = build_design(arrays, indices, fit_terms)
        xtar_values = evaluate_terms(arrays, indices, xtar_rows)
        xtar_prediction = xtar_values @ xtar_coefficients

        targets = np.column_stack(
            (
                arrays["xptarT"][indices] - xtar_prediction[:, 0],
                (arrays["ytarT"][indices] - 100.0 * xtar_prediction[:, 1]) / 100.0,
                arrays["yptarT"][indices] - xtar_prediction[:, 2],
            )
        )
        design_parts.append(design)
        target_parts.append(targets)
        setting_parts.append(np.full(len(indices), setting_index, dtype=np.int16))
        selected_summary.append((rungroup, len(indices), per_foil))
        global_count += len(indices)

    return (
        np.vstack(design_parts),
        np.vstack(target_parts),
        np.concatenate(setting_parts),
        fit_terms,
        selected_summary,
        len(matrix_rows),
        len(xtar_rows),
    )


def root_svd(gram, rhs=None):
    """Run the same ROOT TDecompSVD class used by the production fit."""
    nterms = gram.shape[0]
    contiguous = np.ascontiguousarray(gram, dtype=np.float64)
    matrix = ROOT.TMatrixD(nterms, nterms, contiguous.ravel())
    solver = ROOT.TDecompSVD(matrix)
    singular = np.asarray(
        [solver.GetSig()[index] for index in range(nterms)], dtype=np.float64
    )
    tolerance = float(solver.GetTol())
    threshold = tolerance * singular[0]
    rank = int(np.count_nonzero(singular > threshold))
    if rhs is None:
        return None, singular, threshold, rank, tolerance

    solution = np.empty_like(rhs)
    for target in range(rhs.shape[1]):
        vector_data = np.ascontiguousarray(rhs[:, target], dtype=np.float64)
        vector = ROOT.TVectorD(nterms, vector_data)
        if not solver.Solve(vector):
            raise RuntimeError(f"ROOT TDecompSVD solve failed at N={nterms}")
        solution[:, target] = [vector[index] for index in range(nterms)]
    return solution, singular, threshold, rank, tolerance


def signed_margin(value, threshold):
    if value <= 0.0 or threshold <= 0.0:
        return -math.inf
    return math.log10(value / threshold)


def graph(x, y, color, marker, title=""):
    result = ROOT.TGraph(len(x), array("d", x), array("d", y))
    result.SetLineColor(color)
    result.SetMarkerColor(color)
    result.SetMarkerStyle(marker)
    result.SetMarkerSize(0.55)
    result.SetLineWidth(2)
    if title:
        result.SetTitle(title)
    return result


def save_canvas(canvas, output: Path):
    canvas.SaveAs(str(output.with_suffix(".png")))
    canvas.SaveAs(str(output.with_suffix(".pdf")))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--file-id", type=int, default=-1)
    parser.add_argument("--nfit-max", type=int, default=200000)
    parser.add_argument("--step", type=int, default=1)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    design, targets, setting_index, terms, selected_summary, valid_rows, xtar_rows = load_problem(
        args.campaign, args.metadata, args.file_id, args.nfit_max
    )
    event_count, maximum_terms = design.shape
    gram_full = design.T @ design
    rhs_full = design.T @ targets
    target_norm2 = np.sum(targets * targets, axis=0)
    column_norms = np.sqrt(np.diag(gram_full))
    tolerance = np.finfo(np.float64).eps

    requested = list(range(1, maximum_terms + 1, max(args.step, 1)))
    if requested[-1] != maximum_terms:
        requested.append(maximum_terms)

    records = []
    solutions = []
    previous_prediction_coefficients = None
    physical_scales = np.asarray([1000.0, 100.0, 1000.0])

    for nterms in requested:
        gram = gram_full[:nterms, :nterms]
        rhs = rhs_full[:nterms]
        solution, singular, threshold, rank, root_tolerance = root_svd(gram, rhs)
        if root_tolerance != tolerance:
            raise RuntimeError(
                f"Unexpected ROOT tolerance at N={nterms}: {root_tolerance}"
            )

        norms = column_norms[:nterms]
        scaled_gram = gram / np.outer(norms, norms)
        _, scaled_singular, scaled_threshold, scaled_rank, scaled_tolerance = root_svd(
            scaled_gram
        )
        if scaled_tolerance != tolerance:
            raise RuntimeError(
                f"Unexpected scaled ROOT tolerance at N={nterms}: {scaled_tolerance}"
            )

        residual_variance = (
            target_norm2
            - 2.0 * np.sum(solution * rhs, axis=0)
            + np.einsum("ik,ij,jk->k", solution, gram, solution)
        ) / event_count
        residual_rms = np.sqrt(np.maximum(residual_variance, 0.0)) * physical_scales

        embedded = np.zeros((maximum_terms, 3), dtype=np.float64)
        embedded[:nterms] = solution
        if previous_prediction_coefficients is None:
            step_rms = np.full(3, np.nan)
        else:
            difference = embedded - previous_prediction_coefficients
            step_variance = np.einsum(
                "ik,ij,jk->k", difference, gram_full, difference
            ) / event_count
            step_rms = np.sqrt(np.maximum(step_variance, 0.0)) * physical_scales
        previous_prediction_coefficients = embedded
        solutions.append(embedded)

        smallest_retained = singular[rank - 1] if rank else math.nan
        code = "".join(str(value) for value in terms[nterms - 1])
        records.append(
            {
                "N": nterms,
                "last_code": code,
                "degree": sum(terms[nterms - 1]),
                "rank": rank,
                "truncated": nterms - rank,
                "sigma_max_Ay": singular[0],
                "sigma_min_Ay": singular[-1],
                "threshold_Ay": threshold,
                "sigma_min_minus_threshold_Ay": singular[-1] - threshold,
                "log10_min_over_threshold_Ay": signed_margin(singular[-1], threshold),
                "log10_min_retained_over_threshold_Ay": signed_margin(
                    smallest_retained, threshold
                ),
                "kappa_Ay": singular[0] / singular[-1] if singular[-1] > 0 else math.inf,
                "rank_scaled": scaled_rank,
                "truncated_scaled": nterms - scaled_rank,
                "threshold_scaled_Gram": scaled_threshold,
                "log10_min_over_threshold_scaled_Gram": signed_margin(
                    scaled_singular[-1], scaled_threshold
                ),
                "kappa_scaled_Gram": (
                    scaled_singular[0] / scaled_singular[-1]
                    if scaled_singular[-1] > 0 else math.inf
                ),
                "xptar_residual_rms_mrad": residual_rms[0],
                "ytar_residual_rms_cm": residual_rms[1],
                "yptar_residual_rms_mrad": residual_rms[2],
                "xptar_step_prediction_rms_mrad": step_rms[0],
                "ytar_step_prediction_rms_cm": step_rms[1],
                "yptar_step_prediction_rms_mrad": step_rms[2],
            }
        )

    full_solution = solutions[-1]
    for record, solution in zip(records, solutions):
        difference = solution - full_solution
        variance = np.einsum("ik,ij,jk->k", difference, gram_full, difference) / event_count
        rms = np.sqrt(np.maximum(variance, 0.0)) * physical_scales
        record["xptar_prediction_diff_from_full_rms_mrad"] = rms[0]
        record["ytar_prediction_diff_from_full_rms_cm"] = rms[1]
        record["yptar_prediction_diff_from_full_rms_mrad"] = rms[2]

    table_path = args.output / "angular_term_ladder.tsv"
    fields = list(records[0].keys())
    with table_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    coefficient_path = args.output / "angular_term_ladder_coefficients.npz"
    np.savez_compressed(
        coefficient_path,
        N=np.asarray(requested, dtype=np.int32),
        coefficients=np.stack(solutions),
        term_codes=np.asarray(["".join(map(str, term)) for term in terms]),
    )

    n_values = [float(record["N"]) for record in records]
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    canvas = ROOT.TCanvas("c_margin", "Solve-specific truncation margin", 1050, 760)
    margin_ay = graph(
        n_values,
        [record["log10_min_over_threshold_Ay"] for record in records],
        ROOT.kBlue + 1,
        20,
    )
    margin_scaled = graph(
        n_values,
        [record["log10_min_over_threshold_scaled_Gram"] for record in records],
        ROOT.kRed + 1,
        24,
    )
    multi = ROOT.TMultiGraph()
    multi.SetTitle(
        "Distance from each solve's own SVD truncation threshold;"
        "Number of fitted terms N;log_{10}(#sigma_{min}) - log_{10}(#tau_{N})  [decades]"
    )
    multi.Add(margin_ay, "LP")
    multi.Add(margin_scaled, "LP")
    multi.Draw("A")
    zero = ROOT.TLine(min(n_values), 0.0, max(n_values), 0.0)
    zero.SetLineStyle(2)
    zero.SetLineColor(ROOT.kGray + 2)
    zero.Draw()
    legend = ROOT.TLegend(0.14, 0.15, 0.55, 0.27)
    legend.AddEntry(margin_ay, "Production G_{N}=X_{N}^{T}X_{N}", "lp")
    legend.AddEntry(margin_scaled, "Unit-column Gram matrix", "lp")
    legend.AddEntry(zero, "zero margin (solve-specific threshold)", "l")
    legend.Draw()
    save_canvas(canvas, args.output / "truncation_margin_vs_N")

    low_records = [record for record in records if record["N"] <= 70]
    low_n = [float(record["N"]) for record in low_records]
    canvas_margin_low = ROOT.TCanvas(
        "c_margin_low", "Low-N solve-specific truncation margin", 1050, 760
    )
    margin_ay_low = graph(
        low_n,
        [record["log10_min_over_threshold_Ay"] for record in low_records],
        ROOT.kBlue + 1,
        20,
    )
    margin_scaled_low = graph(
        low_n,
        [record["log10_min_over_threshold_scaled_Gram"] for record in low_records],
        ROOT.kRed + 1,
        24,
    )
    multi_margin_low = ROOT.TMultiGraph()
    multi_margin_low.SetTitle(
        "Low-N distance from each solve's own SVD threshold;"
        "Number of fitted terms N;log_{10}(#sigma_{min}) - log_{10}(#tau_{N})  [decades]"
    )
    multi_margin_low.Add(margin_ay_low, "LP")
    multi_margin_low.Add(margin_scaled_low, "LP")
    multi_margin_low.Draw("A")
    zero_low = ROOT.TLine(min(low_n), 0.0, max(low_n), 0.0)
    zero_low.SetLineStyle(2)
    zero_low.SetLineColor(ROOT.kGray + 2)
    zero_low.Draw()
    legend_margin_low = ROOT.TLegend(0.14, 0.15, 0.55, 0.27)
    legend_margin_low.AddEntry(margin_ay_low, "Production G_{N}=X_{N}^{T}X_{N}", "lp")
    legend_margin_low.AddEntry(margin_scaled_low, "Unit-column Gram matrix", "lp")
    legend_margin_low.AddEntry(zero_low, "zero margin (solve-specific threshold)", "l")
    legend_margin_low.Draw()
    save_canvas(canvas_margin_low, args.output / "truncation_margin_lowN")

    canvas_condition = ROOT.TCanvas("c_condition", "Condition ladder", 1050, 760)
    canvas_condition.SetLogy()
    condition_raw = graph(
        n_values, [record["kappa_Ay"] for record in records], ROOT.kBlue + 1, 20
    )
    condition_scaled = graph(
        n_values,
        [record["kappa_scaled_Gram"] for record in records],
        ROOT.kRed + 1,
        24,
    )
    condition_multi = ROOT.TMultiGraph()
    condition_multi.SetTitle(
        "Angular term ladder conditioning;Number of fitted terms N;Condition number"
    )
    condition_multi.Add(condition_raw, "LP")
    condition_multi.Add(condition_scaled, "LP")
    condition_multi.Draw("A")
    condition_legend = ROOT.TLegend(0.14, 0.72, 0.48, 0.84)
    condition_legend.AddEntry(condition_raw, "Production Gram matrix", "lp")
    condition_legend.AddEntry(condition_scaled, "Unit-column Gram matrix", "lp")
    condition_legend.Draw()
    save_canvas(canvas_condition, args.output / "condition_vs_N")

    canvas_condition_low = ROOT.TCanvas(
        "c_condition_low", "Low-N condition ladder", 1050, 760
    )
    canvas_condition_low.SetLogy()
    condition_raw_low = graph(
        low_n, [record["kappa_Ay"] for record in low_records], ROOT.kBlue + 1, 20
    )
    condition_scaled_low = graph(
        low_n,
        [record["kappa_scaled_Gram"] for record in low_records],
        ROOT.kRed + 1,
        24,
    )
    condition_multi_low = ROOT.TMultiGraph()
    condition_multi_low.SetTitle(
        "Low-N angular conditioning;Number of fitted terms N;Condition number"
    )
    condition_multi_low.Add(condition_raw_low, "LP")
    condition_multi_low.Add(condition_scaled_low, "LP")
    condition_multi_low.Draw("A")
    condition_legend_low = ROOT.TLegend(0.14, 0.72, 0.48, 0.84)
    condition_legend_low.AddEntry(condition_raw_low, "Production Gram matrix", "lp")
    condition_legend_low.AddEntry(condition_scaled_low, "Unit-column Gram matrix", "lp")
    condition_legend_low.Draw()
    save_canvas(canvas_condition_low, args.output / "condition_lowN")

    canvas_step = ROOT.TCanvas("c_step", "Prediction stability", 1050, 900)
    canvas_step.Divide(1, 3)
    step_specs = (
        ("xptar_step_prediction_rms_mrad", "x'_{tar}", "mrad", ROOT.kBlue + 1),
        ("ytar_step_prediction_rms_cm", "y_{tar}", "cm", ROOT.kGreen + 2),
        ("yptar_step_prediction_rms_mrad", "y'_{tar}", "mrad", ROOT.kRed + 1),
    )
    step_graphs = []
    for pad, (field, label, units, color) in enumerate(step_specs, start=1):
        canvas_step.cd(pad)
        ROOT.gPad.SetLogy()
        values = [
            max(float(record[field]), 1e-15) if math.isfinite(float(record[field])) else 1e-15
            for record in records[1:]
        ]
        item = graph(
            n_values[1:],
            values,
            color,
            20,
            f"Change from previous ladder solve: {label};N;prediction RMS difference ({units})",
        )
        item.Draw("ALP")
        step_graphs.append(item)
    save_canvas(canvas_step, args.output / "prediction_step_vs_N")

    canvas_residual = ROOT.TCanvas("c_residual", "Fit residual ladder", 1050, 900)
    canvas_residual.Divide(1, 3)
    residual_specs = (
        ("xptar_residual_rms_mrad", "x'_{tar}", "mrad", ROOT.kBlue + 1),
        ("ytar_residual_rms_cm", "y_{tar}", "cm", ROOT.kGreen + 2),
        ("yptar_residual_rms_mrad", "y'_{tar}", "mrad", ROOT.kRed + 1),
    )
    residual_graphs = []
    for pad, (field, label, units, color) in enumerate(residual_specs, start=1):
        canvas_residual.cd(pad)
        item = graph(
            n_values,
            [float(record[field]) for record in records],
            color,
            20,
            f"Selected-sample residual: {label};N;residual RMS ({units})",
        )
        item.Draw("ALP")
        residual_graphs.append(item)
    save_canvas(canvas_residual, args.output / "residual_rms_vs_N")

    first_truncated = next((record for record in records if record["truncated"]), None)
    first_scaled_truncated = next(
        (record for record in records if record["truncated_scaled"]), None
    )
    report_lines = [
        "Preliminary corrected angular term ladder",
        f"events = {event_count}",
        f"maximum terms = {maximum_terms}",
        f"evaluated rungs = {len(records)} (step={args.step})",
        f"valid seed-matrix rows = {valid_rows}",
        f"retained xtar rows after malformed-line rejection = {xtar_rows}",
        f"ROOT-compatible relative SVD tolerance = {tolerance:.17g}",
        "retention rule = singular_value > tolerance * largest_singular_value",
        "truncation metric = log10(sigma_min) - log10(solve-specific threshold)",
        "zero margin means the weakest mode equals that rung's own threshold",
        "",
        "Malformed-row ambiguity retained for provenance:",
        "  The original fit accepted a blank seed-matrix line as an uninitialized xtar term.",
        "  This ladder excludes malformed rows and uses 461 valid seed rows / 252 xtar rows.",
        "",
        "Selected event summary:",
    ]
    for label, count, per_foil in selected_summary:
        report_lines.append(f"  {label}: {count}, per foil {per_foil}")
    report_lines.extend(
        [
            "",
            "First production-Gram truncation: "
            + (
                f"N={first_truncated['N']} ({first_truncated['truncated']} modes)"
                if first_truncated else "none"
            ),
            "First unit-column-Gram truncation: "
            + (
                f"N={first_scaled_truncated['N']} ({first_scaled_truncated['truncated_scaled']} modes)"
                if first_scaled_truncated else "none"
            ),
            "",
            f"Full N={maximum_terms} production rank = {records[-1]['rank']}/{maximum_terms}",
            f"Full N={maximum_terms} production truncation margin = "
            f"{records[-1]['log10_min_over_threshold_Ay']:.9g} decades",
            f"Full N={maximum_terms} unit-column rank = "
            f"{records[-1]['rank_scaled']}/{maximum_terms}",
            f"Full N={maximum_terms} unit-column truncation margin = "
            f"{records[-1]['log10_min_over_threshold_scaled_Gram']:.9g} decades",
            "",
            f"Table: {table_path}",
            f"Coefficients: {coefficient_path}",
        ]
    )
    report_path = args.output / "angular_term_ladder_summary.txt"
    report_path.write_text("\n".join(report_lines) + "\n")
    print(report_path.read_text())


if __name__ == "__main__":
    main()
