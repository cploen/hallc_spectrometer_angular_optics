#!/usr/bin/env python3

import argparse
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# Presentation-readable defaults.
plt.rcParams.update({
    "font.size": 17,
    "axes.titlesize": 21,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 13,
    "figure.titlesize": 21,
})


def classify_arrangement(foils):
    vals = sorted(round(float(z), 3) for z in foils)

    if vals == [0.0]:
        return "0 cm"
    if vals == [-8.0, 8.0]:
        return "±8 cm"
    if vals == [-3.0, 3.0]:
        return "±3 cm"

    return ",".join(str(v) for v in vals)


def add_labels(ax, xvals, yvals, labels, dx=0.03, dy=0.0):
    for x, y, label in zip(xvals, yvals, labels):
        ax.text(x + dx, y + dy, str(label), fontsize=12)


def add_labels_staggered(ax, xvals, yvals, labels, dx=0.18, dy_step=0.18):
    """
    Stagger run-number labels aggressively so same-foil points stay readable.
    """
    x_offsets = [dx, dx, -1.15 * dx, 1.8 * dx, -1.8 * dx, 2.5 * dx, -2.5 * dx]
    y_offsets = [0.0, dy_step, -dy_step, 2.0 * dy_step, -2.0 * dy_step, 3.0 * dy_step, -3.0 * dy_step]

    for i, (x, y, label) in enumerate(zip(xvals, yvals, labels)):
        ax.text(
            x + x_offsets[i % len(x_offsets)],
            y + y_offsets[i % len(y_offsets)],
            str(label),
            fontsize=14,
            fontweight="bold",
        )


def add_hms_angle(df):
    angle_map = {
        # 6.667 GeV/c, HMS angle 12.490
        1537: 12.490,
        1538: 12.490,
        1539: 12.490,
        1540: 12.490,
        1541: 12.490,
        1542: 12.490,
        1543: 12.490,
        1544: 12.490,

        # 6.667 GeV/c, HMS angle 15.195
        51600: 15.195,
        5161: 15.195,
        5162: 15.195,
        5163: 15.195,
        51640: 15.195,

        # 6.667 GeV/c, HMS angle 12.500
        69820: 12.500,
        69840: 12.500,
        6985: 12.500,
    }

    df["hms_angle_deg"] = df["run"].map(angle_map)
    return df


def savefig(fig, stem, combined_pdf):
    # Larger padding prevents big axis labels from being clipped in the PDF.
    fig.tight_layout(pad=2.0)
    fig.subplots_adjust(left=0.18, bottom=0.18)
    combined_pdf.savefig(fig)
    print(f"Added {stem} to combined PDF")


def add_mean_error_column(df):
    """
    Add statistical error on the fitted mean.

    Preferred input, if available:
      fit_mean_err_cm

    Fallback:
      fit_sigma_cm / sqrt(entries)

    This keeps fitted peak width / RMS separate from uncertainty on the centroid.
    """
    if "fit_mean_err_cm" in df.columns:
        return df

    required = {"fit_sigma_cm", "entries"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Cannot compute error on mean; missing columns: {sorted(missing)}")

    df = df.copy()
    entries = df["entries"].clip(lower=1)
    df["fit_mean_err_cm"] = df["fit_sigma_cm"] / np.sqrt(entries)
    return df


def tilted_errorbar(ax, x, y, yerr, dx=0.035, **kwargs):
    """
    Draw compact diagonal error bars.

    The vertical span is still y +/- yerr, but the segment is drawn with
    a slight x-offset so it reads as a slash instead of a vertical whisker.
    This reduces clutter in crowded comparison plots.
    """
    segments = []
    for xi, yi, ei in zip(x, y, yerr):
        if not np.isfinite(xi) or not np.isfinite(yi) or not np.isfinite(ei):
            continue
        segments.append([(xi - dx, yi - ei), (xi + dx, yi + ei)])

    if not segments:
        return None

    lc = LineCollection(segments, **kwargs)
    ax.add_collection(lc)
    return lc



def set_clean_foil_ticks(ax, lo=-8, hi=8, step=2):
    """
    Use physically meaningful foil-position ticks instead of matplotlib's
    default half-step ticks. This keeps ±8 on labeled ticks and places ±3
    cleanly between ±2 and ±4.
    """
    ticks = np.arange(lo, hi + 0.5 * step, step)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)



def make_axes_readable(ax):
    """Increase tick size/weight and make plots easier to read in slides."""
    ax.tick_params(axis="both", which="major", labelsize=16, width=1.2, length=6)
    ax.tick_params(axis="both", which="minor", width=1.0, length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def main():
    parser = argparse.ArgumentParser(
        description="Plot ztar peak stability from TSV produced by fit_ztar_peaks_from_replay.C"
    )
    parser.add_argument("tsv", help="Input TSV file")
    parser.add_argument(
        "--outdir",
        default="plots/ztar_stability",
        help="Output directory [default: plots/ztar_stability]",
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="Optional output tag. Default uses TSV filename stem.",
    )
    args = parser.parse_args()

    tsv_path = Path(args.tsv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tag = args.tag if args.tag else tsv_path.stem
    combined_pdf_path = outdir / f"{tag}_ALL_PLOTS.pdf"
    pair_tsv = outdir / f"{tag}_pair_summary.tsv"

    df = pd.read_csv(tsv_path, sep="\t")
    df = add_mean_error_column(df)

    # Write a compact table separating centroid uncertainty from peak width.
    df["fit_sigma_mm"] = 10.0 * df["fit_sigma_cm"]
    df["fit_mean_err_mm"] = 10.0 * df["fit_mean_err_cm"]
    df["mean_minus_nominal_mm"] = 10.0 * df["mean_minus_nominal_cm"]

    summary_cols = [
        "run",
        "opticsID",
        "nominal_foil_z_cm",
        "fit_mean_cm",
        "mean_minus_nominal_cm",
        "mean_minus_nominal_mm",
        "fit_sigma_cm",
        "fit_sigma_mm",
        "fit_mean_err_cm",
        "fit_mean_err_mm",
        "entries",
    ]

    summary_tsv = outdir / f"{tag}_fit_summary_errors.tsv"
    df[summary_cols].to_csv(summary_tsv, sep="\t", index=False)
    print(f"Wrote {summary_tsv}")

    required = {
        "run",
        "opticsID",
        "nominal_foil_z_cm",
        "fit_mean_cm",
        "fit_sigma_cm",
        "mean_minus_nominal_cm",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    print("\n=== DataFrame preview ===")
    print(df.head())

    print("\n=== Columns ===")
    print(df.columns.tolist())

    print("\n=== Shape ===")
    print(df.shape)

    print("\n=== Runs / foil positions ===")
    print(
        df[
            [
                "run",
                "opticsID",
                "nominal_foil_z_cm",
                "fit_mean_cm",
                "fit_sigma_cm",
                "mean_minus_nominal_cm",
            ]
        ]
    )

    # Derived dataframe columns.
    df["run_label"] = df["run"].astype(str)

    arrangement_map = (
        df.groupby("run")["nominal_foil_z_cm"]
        .apply(lambda s: classify_arrangement(s.tolist()))
        .to_dict()
    )
    df["arrangement"] = df["run"].map(arrangement_map)
    df["group_label"] = df["run_label"] + " (" + df["arrangement"] + ")"
    df = add_hms_angle(df)

    arrangements = ["0 cm", "±3 cm", "±8 cm"]
    markers = {"0 cm": "o", "±3 cm": "s", "±8 cm": "^"}

    with PdfPages(combined_pdf_path) as combined_pdf:

        # 1. Residual vs nominal foil position.
        fig, ax = plt.subplots(figsize=(8, 5))

        for arr in arrangements:
            sub = df[df["arrangement"] == arr]
            if sub.empty:
                continue

            ax.scatter(
                sub["nominal_foil_z_cm"],
                sub["mean_minus_nominal_cm"],
                marker="s",
                s=120,
                color="red",
            )
    
            add_labels(
                ax,
                sub["nominal_foil_z_cm"],
                sub["mean_minus_nominal_cm"],
                sub["run_label"],
                dx=0.08,
            )

        ax.axhline(0, linewidth=1.6)
        ax.set_xlabel("Nominal foil z [cm]")
        ax.set_ylabel("Fit mean - nominal z [cm]")
        ax.set_title("Ztar residual vs nominal foil position")
        ax.grid(True, alpha=0.3)

        make_axes_readable(ax)
        savefig(fig, f"{tag}_residual_vs_foil_position", combined_pdf)
        plt.close(fig)

        # 2. Fitted mean vs nominal foil position.
        fig, ax = plt.subplots(figsize=(8, 5))

        xmin = min(df["nominal_foil_z_cm"].min(), df["fit_mean_cm"].min()) - 0.5
        xmax = max(df["nominal_foil_z_cm"].max(), df["fit_mean_cm"].max()) + 0.5

        ax.plot(
            [xmin, xmax],
            [xmin, xmax],
            linestyle="--",
            linewidth=1,
            
        )

        for arr in arrangements:
            sub = df[df["arrangement"] == arr]
            if sub.empty:
                continue

            ax.scatter(
                sub["nominal_foil_z_cm"],
                sub["fit_mean_cm"],
                marker="s",
                s=120,
                color="red",
            )
            tilted_errorbar(
                ax,
                sub["nominal_foil_z_cm"].to_numpy(),
                sub["fit_mean_cm"].to_numpy(),
                sub["fit_mean_err_cm"].to_numpy(),
                dx=0.035,
                linewidths=1.4,
                alpha=0.9,
                colors="red",
            )

            add_labels_staggered(
                ax,
                sub["nominal_foil_z_cm"],
                sub["fit_mean_cm"],
                sub["run_label"],
                dx=0.24,
                dy_step=0.22,
            )

        set_clean_foil_ticks(ax)
        ax.set_xlabel("Nominal foil z [cm]")
        ax.set_ylabel("Fitted H.react.z peak mean [cm]")
        ax.set_title("Fitted z peak mean vs nominal foil position")
        ax.grid(True, alpha=0.3)

        make_axes_readable(ax)
        savefig(fig, f"{tag}_mean_vs_nominal", combined_pdf)
        plt.close(fig)

        # 3. Residual by run / foil.
        fig, ax = plt.subplots(figsize=(10, 5))

        df_plot = df.copy()
        df_plot["x_label"] = (
            df_plot["run"].astype(str)
            + " z="
            + df_plot["nominal_foil_z_cm"].astype(str)
        )

        x = range(len(df_plot))
        ax.axhline(0, linewidth=1.6)

        for arr in arrangements:
            sub = df_plot[df_plot["arrangement"] == arr]
            if sub.empty:
                continue

            xpos = [df_plot.index.get_loc(i) for i in sub.index]

            ax.scatter(
                xpos,
                sub["mean_minus_nominal_cm"],
                marker="s",
                s=120,
                color="red",
            )
            tilted_errorbar(
                ax,
                xpos,
                sub["mean_minus_nominal_cm"].to_numpy(),
                sub["fit_mean_err_cm"].to_numpy(),
                dx=0.035,
                linewidths=1.4,
                alpha=0.9,
                colors="red",
            )

        ax.set_xticks(list(x))
        ax.set_xticklabels(df_plot["x_label"], rotation=60, ha="right", fontsize=10)
        ax.set_ylabel("Fit mean - nominal z [cm]")
        ax.set_title("Ztar residual by run and foil")
        ax.grid(True, alpha=0.3)

        make_axes_readable(ax)
        savefig(fig, f"{tag}_residual_by_run_and_foil", combined_pdf)
        plt.close(fig)

        # Pair-derived table: two-foil runs.
        pair_rows = []

        for run, sub in df.groupby("run"):
            sub = sub.sort_values("nominal_foil_z_cm")
            foils = sub["nominal_foil_z_cm"].tolist()
            arr = classify_arrangement(foils)

            if len(sub) != 2:
                continue

            minus = sub.iloc[0]
            plus = sub.iloc[1]

            mean_minus = float(minus["fit_mean_cm"])
            mean_plus = float(plus["fit_mean_cm"])

            pair_center = 0.5 * (mean_plus + mean_minus)
            pair_sep = mean_plus - mean_minus
            expected_sep = float(
                plus["nominal_foil_z_cm"] - minus["nominal_foil_z_cm"]
            )

            pair_rows.append(
                {
                    "run": int(run),
                    "opticsID": str(minus["opticsID"]),
                    "arrangement": arr,
                    "nominal_minus_cm": float(minus["nominal_foil_z_cm"]),
                    "nominal_plus_cm": float(plus["nominal_foil_z_cm"]),
                    "fit_mean_minus_cm": mean_minus,
                    "fit_mean_plus_cm": mean_plus,
                    "sigma_minus_cm": float(minus["fit_sigma_cm"]),
                    "sigma_plus_cm": float(plus["fit_sigma_cm"]),
                    "pair_center_cm": pair_center,
                    "pair_separation_cm": pair_sep,
                    "expected_separation_cm": expected_sep,
                    "separation_residual_cm": pair_sep - expected_sep,
                    "foil_abs_cm": 0.5 * expected_sep,
                    "hms_angle_deg": (
                        float(minus["hms_angle_deg"])
                        if pd.notna(minus["hms_angle_deg"])
                        else float("nan")
                    ),
                }
            )

        pair_df = pd.DataFrame(pair_rows)

        # Add single-foil z=0 runs as "center" points.
        zero_rows = []
        zero_df = df[df["nominal_foil_z_cm"].abs() < 1e-6].copy()

        for _, row in zero_df.iterrows():
            zero_rows.append(
                {
                    "run": int(row["run"]),
                    "opticsID": str(row["opticsID"]),
                    "arrangement": "0 cm",
                    "nominal_minus_cm": 0.0,
                    "nominal_plus_cm": 0.0,
                    "fit_mean_minus_cm": float(row["fit_mean_cm"]),
                    "fit_mean_plus_cm": float(row["fit_mean_cm"]),
                    "sigma_minus_cm": float(row["fit_sigma_cm"]),
                    "sigma_plus_cm": float(row["fit_sigma_cm"]),
                    "pair_center_cm": float(row["fit_mean_cm"]),
                    "pair_separation_cm": 0.0,
                    "expected_separation_cm": 0.0,
                    "separation_residual_cm": 0.0,
                    "foil_abs_cm": 0.0,
                    "hms_angle_deg": (
                        float(row["hms_angle_deg"])
                        if pd.notna(row["hms_angle_deg"])
                        else float("nan")
                    ),
                }
            )

        zero_pair_df = pd.DataFrame(zero_rows)

        if not zero_pair_df.empty:
            pair_df = pd.concat([pair_df, zero_pair_df], ignore_index=True)

        if not pair_df.empty:
            pair_df["run_label"] = pair_df["run"].astype(str)
            pair_df["x_label"] = (
                pair_df["run_label"] + " " + pair_df["arrangement"]
            )

        pair_df.to_csv(pair_tsv, sep="\t", index=False)
        print(f"Wrote {pair_tsv}")

        # 4. Normalized pair separation residual.
        # Exclude 0 cm single-foil rows because expected separation is zero.
        if not pair_df.empty:
            sep_df = pair_df[pair_df["expected_separation_cm"] > 0].copy()

            if not sep_df.empty:
                fig, ax = plt.subplots(figsize=(10.5, 6.2))
                x = range(len(sep_df))

                sep_df["separation_residual_percent"] = (
                    100.0
                    * sep_df["separation_residual_cm"]
                    / sep_df["expected_separation_cm"]
                )

                ax.scatter(
                    list(x),
                    sep_df["separation_residual_percent"],
                    marker="s",
                    s=120,
                    color="red",
                )

                for i, row in sep_df.reset_index(drop=True).iterrows():
                    ax.text(
                        i + 0.03,
                        row["separation_residual_percent"],
                        f'{row["separation_residual_percent"]:.2f}%',
                        fontsize=12,
                    )

                ax.axhline(0, linestyle="--", linewidth=1)
                ax.set_xticks(list(x))
                ax.set_xticklabels(sep_df["x_label"], rotation=35, ha="right")
                ax.set_ylabel("Sep. residual [%]", fontsize=15)
                ax.set_title("Normalized foil-pair separation residual")
                ax.grid(True, alpha=0.3)

                savefig(
                    fig,
                    f"{tag}_pair_separation_residual_percent",
                    combined_pdf,
                )
                plt.close(fig)

            # 5. Pair/single center vs absolute foil distance.
            fig, ax = plt.subplots(figsize=(10.5, 6.2))

            markers_by_foil = {
                0.0: "o",
                3.0: "s",
                8.0: "^",
            }

            for foil_abs in sorted(pair_df["foil_abs_cm"].dropna().unique()):
                sub = pair_df[pair_df["foil_abs_cm"] == foil_abs]

                ax.scatter(
                    sub["foil_abs_cm"],
                    sub["pair_center_cm"],
                    marker="s",
                    s=120,
                    color="red",
                )

                for _, row in sub.iterrows():
                    ax.text(
                        row["foil_abs_cm"] + 0.08,
                        row["pair_center_cm"],
                        str(int(row["run"])),
                        fontsize=14,
                        fontweight="bold",
                    )

            ax.axhline(0, linestyle="--", linewidth=1.6)
            ax.set_xticks([0, 3, 8])
            ax.set_xlabel("Absolute foil distance from nominal center [cm]")
            ax.set_ylabel("Fitted center [cm]")
            ax.set_title("Ztar center stability vs foil distance")
            ax.grid(True, alpha=0.3)

            make_axes_readable(ax)
            savefig(fig, f"{tag}_center_vs_foil_distance", combined_pdf)
            plt.close(fig)

            # 6. Pair/single center vs HMS angle.
            fig, ax = plt.subplots(figsize=(10.5, 6.2))

            for foil_abs in sorted(pair_df["foil_abs_cm"].dropna().unique()):
                sub = pair_df[pair_df["foil_abs_cm"] == foil_abs]
                if sub.empty:
                    continue

                ax.scatter(
                    sub["hms_angle_deg"],
                    sub["pair_center_cm"],
                    marker="s",
                    s=120,
                    color="red",
                )

                for j, (_, row) in enumerate(sub.iterrows()):
                    # Larger offsets are intentional: several runs share nearly identical HMS angle.
                    x_offsets = [0.030, 0.115, -0.105, 0.200, -0.190, 0.285, -0.275]
                    y_offsets = [0.000, 0.090, -0.090, 0.180, -0.180, 0.270, -0.270]
                    ax.text(
                        row["hms_angle_deg"] + x_offsets[j % len(x_offsets)],
                        row["pair_center_cm"] + y_offsets[j % len(y_offsets)],
                        str(int(row["run"])),
                        fontsize=14,
                        fontweight="bold",
                    )

            ax.axhline(0, linestyle="--", linewidth=1.6)
            ax.set_xlim(pair_df["hms_angle_deg"].min() - 0.45, pair_df["hms_angle_deg"].max() + 0.45)
            ax.set_xlabel("HMS angle [deg]")
            ax.set_ylabel("Fitted center [cm]")
            ax.set_title("Ztar center stability vs HMS angle")
            ax.grid(True, alpha=0.3)

            make_axes_readable(ax)
            savefig(fig, f"{tag}_center_vs_hms_angle", combined_pdf)
            plt.close(fig)

    print(f"Wrote combined PDF: {combined_pdf_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
