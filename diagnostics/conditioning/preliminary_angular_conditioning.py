#!/usr/bin/env python3
"""Preliminary conditioning study for the selected HMS angular-fit sample.

This intentionally mirrors the event-order caps in fit_opt_matrix_gmm.C while
ignoring blank or malformed matrix rows. It does not refit or alter a matrix.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np
import ROOT


COEFF_RE = re.compile(
    r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"(\d{5})\s*$"
)


def read_fit_terms(matrix_path: Path) -> list[tuple[int, int, int, int, int]]:
    terms: list[tuple[int, int, int, int, int]] = [(0, 0, 0, 0, 0)]
    in_data = False
    for line in matrix_path.read_text().splitlines():
        if not in_data:
            if line.startswith(" --------------------------------"):
                in_data = True
            continue
        if line.startswith(" --------------------------------"):
            break
        match = COEFF_RE.match(line)
        if not match:
            continue
        exponents = tuple(int(digit) for digit in match.group(5))
        if exponents[4] == 0:
            terms.append(exponents)
    return terms


def read_output_fit_coefficients(matrix_path: Path, expected_terms):
    coefficients = []
    terms = []
    for line in matrix_path.read_text().splitlines():
        match = COEFF_RE.match(line)
        if not match:
            if coefficients:
                break
            continue
        exponents = tuple(int(digit) for digit in match.group(5))
        if exponents[4] != 0:
            break
        coefficients.append([float(match.group(i)) for i in range(1, 4)])
        terms.append(exponents)
    if terms != expected_terms:
        raise ValueError(
            f"Fit terms in {matrix_path} do not match the intended basis: "
            f"{len(terms)} versus {len(expected_terms)}"
        )
    return np.asarray(coefficients)


def read_campaign_rows(tsv_path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with tsv_path.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            rows.append((row["rungroup"], int(row["optics_id"])))
    return rows


def read_optics_metadata(path: Path, wanted_ids: set[int]):
    lines = path.read_text().splitlines()
    result = {}
    for index, line in enumerate(lines):
        fields = [field.strip() for field in line.split(",")]
        if not fields or not fields[0].isdigit():
            continue
        optics_id = int(fields[0])
        if optics_id not in wanted_ids:
            continue
        nfoils = int(fields[3])
        ndel_edges = int(fields[5])
        foils = [float(value) for value in lines[index + 1].split(",")]
        delta_edges = [float(value) for value in lines[index + 2].split(",")]
        if len(foils) != nfoils or len(delta_edges) != ndel_edges:
            raise ValueError(f"Malformed metadata for optics ID {optics_id}")
        result[optics_id] = (foils, delta_edges)
    missing = wanted_ids - result.keys()
    if missing:
        raise ValueError(f"Missing optics metadata for {sorted(missing)}")
    return result


def selected_indices(arrays, foils, delta_edges, global_before, nfit_max):
    delta = arrays["delta"]
    ztar = arrays["ztarT"]
    ysieve = arrays["ysT"]
    nentries = len(delta)

    foil_index = np.full(nentries, -1, dtype=np.int16)
    for index, center in enumerate(foils):
        foil_index[np.abs(ztar - center) < 2.5] = index

    delta_index = np.full(nentries, -1, dtype=np.int16)
    for index in range(len(delta_edges) - 1):
        mask = (delta >= delta_edges[index]) & (delta < delta_edges[index + 1])
        delta_index[mask] = index

    ys_centers = (np.arange(9) - 4) * 0.6 * 2.54
    ys_index = np.full(nentries, -1, dtype=np.int16)
    for index, center in enumerate(ys_centers):
        ys_index[np.abs(ysieve - center) < 0.5] = index

    eligible = (
        (np.abs(delta) < 100.0)
        & (delta > -15.0)
        & (delta < 30.0)
        & (foil_index >= 0)
        & (delta_index >= 0)
        & (ys_index >= 0)
    )

    per_cell: dict[tuple[int, int, int], int] = {}
    per_foil = [0] * len(foils)
    selected: list[int] = []
    for event in np.flatnonzero(eligible):
        if global_before + len(selected) >= nfit_max:
            break
        foil = int(foil_index[event])
        cell = (foil, int(delta_index[event]), int(ys_index[event]))
        if per_cell.get(cell, 0) >= 1000 or per_foil[foil] >= 15000:
            continue
        per_cell[cell] = per_cell.get(cell, 0) + 1
        per_foil[foil] += 1
        selected.append(int(event))
    return np.asarray(selected, dtype=np.int64), per_foil


def build_design(arrays, indices, terms):
    variables = (
        arrays["xfp"][indices] / 100.0,
        arrays["xpfp"][indices],
        arrays["yfp"][indices] / 100.0,
        arrays["ypfp"][indices],
    )
    design = np.empty((len(indices), len(terms)), dtype=np.float64)
    for column, exponents in enumerate(terms):
        values = np.ones(len(indices), dtype=np.float64)
        for variable, exponent in zip(variables, exponents[:4]):
            if exponent:
                values *= variable**exponent
        design[:, column] = values
    return design


def spectrum_metrics(design, direct=True):
    norms = np.linalg.norm(design, axis=0)
    scaled = design / norms
    if direct:
        singular = np.linalg.svd(design, full_matrices=False, compute_uv=False)
        singular_scaled = np.linalg.svd(scaled, full_matrices=False, compute_uv=False)
    else:
        raw_eigenvalues = np.linalg.eigvalsh(design.T @ design)
        scaled_eigenvalues = np.linalg.eigvalsh(scaled.T @ scaled)
        singular = np.sqrt(np.maximum(raw_eigenvalues[::-1], 0.0))
        singular_scaled = np.sqrt(np.maximum(scaled_eigenvalues[::-1], 0.0))
    tolerance = max(design.shape) * np.finfo(np.float64).eps
    rank = int(np.count_nonzero(singular > tolerance * singular[0]))
    rank_scaled = int(
        np.count_nonzero(singular_scaled > tolerance * singular_scaled[0])
    )
    gaps = singular_scaled[:-1] / singular_scaled[1:]

    gram_scaled = scaled.T @ scaled
    eigenvalues, eigenvectors = np.linalg.eigh(gram_scaled)
    weakest = eigenvectors[:, 0]
    correlations = gram_scaled.copy()
    np.fill_diagonal(correlations, 0.0)
    max_corr = float(np.max(np.abs(correlations)))

    return {
        "norms": norms,
        "singular": singular,
        "singular_scaled": singular_scaled,
        "rank": rank,
        "rank_scaled": rank_scaled,
        "largest_scaled_gap": float(np.max(gaps)),
        "largest_scaled_gap_after_mode": int(np.argmax(gaps)),
        "weakest_scaled_vector": weakest,
        "max_abs_column_correlation": max_corr,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--file-id", type=int, default=-1)
    parser.add_argument("--nfit-max", type=int, default=200000)
    parser.add_argument("--matrix-a", type=Path)
    parser.add_argument("--matrix-b", type=Path)
    args = parser.parse_args()

    config = args.campaign / "config"
    rungroup_files = list(config.glob("rungroups_*_inputs.tsv"))
    if len(rungroup_files) != 1:
        raise SystemExit(f"Expected one rungroup input TSV, found {rungroup_files}")
    rows = read_campaign_rows(rungroup_files[0])
    metadata = read_optics_metadata(args.metadata, {row[1] for row in rows})
    terms = read_fit_terms(config / "oldfit.dat")

    branches = ["delta", "ztarT", "ysT", "xfp", "xpfp", "yfp", "ypfp"]
    designs = []
    setting_labels = []
    selected_counts = []
    global_count = 0
    for rungroup, optics_id in rows:
        tree_path = (
            args.campaign
            / "06a_fit_ntuple/root"
            / f"Optics_{optics_id}_{args.file_id}_fit_tree_gmm.root"
        )
        arrays = ROOT.RDataFrame("TFit", str(tree_path)).AsNumpy(branches)
        indices, per_foil = selected_indices(
            arrays, *metadata[optics_id], global_count, args.nfit_max
        )
        design = build_design(arrays, indices, terms)
        designs.append(design)
        setting_labels.append(rungroup)
        selected_counts.append((len(indices), per_foil))
        global_count += len(indices)

    design = np.vstack(designs)
    metrics = spectrum_metrics(design)

    loo_rows = []
    for index, label in enumerate(setting_labels):
        keep = np.vstack([part for j, part in enumerate(designs) if j != index])
        loo = spectrum_metrics(keep, direct=False)
        loo_rows.append(
            (
                label,
                len(keep),
                loo["singular"][0] / loo["singular"][-1],
                loo["singular_scaled"][0] / loo["singular_scaled"][-1],
                loo["rank"],
                loo["rank_scaled"],
            )
        )

    args.output.mkdir(parents=True, exist_ok=True)
    report = args.output / "preliminary_angular_conditioning.txt"
    spectrum = args.output / "angular_singular_values.tsv"
    loo_file = args.output / "leave_one_setting_out_conditioning.tsv"

    with spectrum.open("w") as stream:
        stream.write("mode\traw_singular_value\tunit_column_singular_value\n")
        for index, (raw, scaled) in enumerate(
            zip(metrics["singular"], metrics["singular_scaled"])
        ):
            stream.write(f"{index}\t{raw:.17g}\t{scaled:.17g}\n")

    with loo_file.open("w") as stream:
        stream.write(
            "omitted_setting\tn_events\traw_condition\tunit_column_condition"
            "\traw_rank\tunit_column_rank\n"
        )
        for row in loo_rows:
            stream.write("\t".join(map(str, row)) + "\n")

    weakest = metrics["weakest_scaled_vector"]
    top = np.argsort(np.abs(weakest))[::-1][:12]
    lines = [
        "Preliminary HMS angular design-matrix conditioning",
        f"events = {len(design)}",
        f"columns = {design.shape[1]}",
        f"raw condition number = {metrics['singular'][0] / metrics['singular'][-1]:.9g}",
        "unit-column condition number = "
        f"{metrics['singular_scaled'][0] / metrics['singular_scaled'][-1]:.9g}",
        f"raw numerical rank = {metrics['rank']}",
        f"unit-column numerical rank = {metrics['rank_scaled']}",
        f"column norm ratio = {metrics['norms'].max() / metrics['norms'].min():.9g}",
        f"max absolute unit-column correlation = {metrics['max_abs_column_correlation']:.9g}",
        "largest adjacent gap in unit-column spectrum = "
        f"{metrics['largest_scaled_gap']:.9g} after mode "
        f"{metrics['largest_scaled_gap_after_mode']}",
        "",
        "Selected events by setting:",
    ]
    for label, (count, per_foil) in zip(setting_labels, selected_counts):
        lines.append(f"  {label}: {count}, per foil {per_foil}")
    lines.extend(["", "Largest components of weakest unit-column mode:"])
    for index in top:
        code = "".join(str(value) for value in terms[index])
        lines.append(f"  column {index:3d} term {code}: {weakest[index]: .9g}")
    lines.extend(["", "Leave-one-setting-out conditioning:"])
    for row in loo_rows:
        lines.append(
            f"  omit {row[0]}: events={row[1]}, raw={row[2]:.9g}, "
            f"unit-column={row[3]:.9g}, ranks={row[4]}/{row[5]}"
        )
    if args.matrix_a and args.matrix_b:
        coefficients_a = read_output_fit_coefficients(args.matrix_a, terms)
        coefficients_b = read_output_fit_coefficients(args.matrix_b, terms)
        difference = coefficients_b - coefficients_a
        lines.extend(
            [
                "",
                "Matrix A/B comparison on the common selected sample:",
                f"  A = {args.matrix_a}",
                f"  B = {args.matrix_b}",
            ]
        )
        for column, (name, scale, units) in enumerate(
            (("xptar", 1000.0, "mrad"), ("ytar", 100.0, "cm"), ("yptar", 1000.0, "mrad"))
        ):
            prediction_difference = design @ difference[:, column] * scale
            lines.append(
                f"  {name}: coefficient relative L2="
                f"{np.linalg.norm(difference[:, column]) / np.linalg.norm(coefficients_a[:, column]):.9g}, "
                f"prediction difference RMS={np.sqrt(np.mean(prediction_difference**2)):.9g} {units}, "
                f"max={np.max(np.abs(prediction_difference)):.9g} {units}"
            )
    report.write_text("\n".join(lines) + "\n")
    print(report.read_text())


if __name__ == "__main__":
    main()
