#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import ROOT

EDGES = [-10, -8, -5, 0, 5, 10]
MACRO = "assign_xfp_xpfp_angleScanBands_split.C"
OPTICS_DAT = Path("DATfiles/list_of_optics_run.dat")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply stored XFP/XPFP theta values from one reference foil "
            "to every foil in a target rungroup."
        )
    )
    parser.add_argument("campaign")
    parser.add_argument("target_rungroup")
    parser.add_argument("reference_rungroup")
    parser.add_argument("reference_foil", type=int)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def find_rungroup_config(campaign: Path) -> Path:
    matches = sorted((campaign / "config").glob("rungroups_*_inputs.tsv"))

    if len(matches) != 1:
        raise RuntimeError(
            f"expected exactly one rungroups_*_inputs.tsv under "
            f"{campaign}/config; found {len(matches)}"
        )

    return matches[0]


def load_rungroups(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as stream:
        return {
            row["rungroup"]: row
            for row in csv.DictReader(stream, delimiter="\t")
        }


def load_num_foils(optics_id: int) -> int:
    for raw in OPTICS_DAT.read_text().splitlines():
        fields = [field.strip() for field in raw.split(",")]

        if len(fields) < 7:
            continue

        if fields[0] != str(optics_id):
            continue

        try:
            return int(fields[3])
        except ValueError as exc:
            raise RuntimeError(
                f"invalid NumFoil for optics ID {optics_id}: {fields[3]!r}"
            ) from exc

    raise RuntimeError(
        f"optics ID {optics_id} not found in {OPTICS_DAT}"
    )


def load_reference_thetas(
    theta_tsv: Path,
    reference_rungroup: str,
    reference_foil: int,
) -> dict[tuple[int, str], float]:
    selected: dict[tuple[int, str], float] = {}

    with theta_tsv.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["tag"] != reference_rungroup:
                continue

            if int(row["foil"]) != reference_foil:
                continue

            key = (int(row["ndel"]), row["zone"])
            selected[key] = float(row["thetaDeg"])

    missing = [
        (ndel, zone)
        for ndel in range(5)
        for zone in ("low", "high")
        if (ndel, zone) not in selected
    ]

    if missing:
        raise RuntimeError(
            f"missing X theta rows for reference={reference_rungroup}, "
            f"foil={reference_foil}: {missing}"
        )

    return selected


def quantile(values: list[float], q: float) -> float | None:
    values = sorted(values)

    if not values:
        return None

    position = q * (len(values) - 1)
    index = int(position)
    fraction = position - index

    if index + 1 < len(values):
        return (
            values[index] * (1.0 - fraction)
            + values[index + 1] * fraction
        )

    return values[index]


def compute_split(
    rootfile: Path,
    ytar_cut_file: Path,
    foil: int,
    delta_min: float,
    delta_max: float,
) -> tuple[float | None, int]:
    input_file = ROOT.TFile.Open(str(rootfile), "READ")

    if not input_file or input_file.IsZombie():
        raise RuntimeError(f"cannot open {rootfile}")

    tree = input_file.Get("T") or input_file.Get("Tout")

    if not tree:
        raise RuntimeError(f"cannot find T or Tout in {rootfile}")

    cut_file = ROOT.TFile.Open(str(ytar_cut_file), "READ")

    if not cut_file or cut_file.IsZombie():
        raise RuntimeError(f"cannot open {ytar_cut_file}")

    ytar_cut = cut_file.Get(f"delta_vs_ytar_cut_foil{foil}")

    if not ytar_cut:
        raise RuntimeError(
            f"cannot find delta_vs_ytar_cut_foil{foil} "
            f"in {ytar_cut_file}"
        )

    values: list[float] = []

    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)

        cer = float(getattr(tree, "H.cer.npeSum"))
        cal = float(getattr(tree, "H.cal.etottracknorm"))
        delta = float(getattr(tree, "H.gtr.dp"))
        ytar = float(getattr(tree, "H.gtr.y"))
        xpfp = float(getattr(tree, "H.dc.xp_fp"))

        if cer <= 2.0 or cal <= 0.65:
            continue

        if not ytar_cut.IsInside(ytar, delta):
            continue

        if delta_min <= delta < delta_max:
            values.append(xpfp)

    input_file.Close()
    cut_file.Close()

    if len(values) < 1000:
        return None, len(values)

    q05 = quantile(values, 0.05)
    q95 = quantile(values, 0.95)

    if q05 is None or q95 is None:
        return None, len(values)

    return 0.5 * (q05 + q95), len(values)


def run_expression(expression: str, dry_run: bool) -> None:
    command = ["hcana", "-b", "-l", "-q", expression]

    if dry_run:
        print(subprocess.list2cmdline(command))
        return

    subprocess.run(command, check=True)


def main() -> None:
    args = parse_args()

    campaign = Path(args.campaign)
    config = find_rungroup_config(campaign)
    rungroups = load_rungroups(config)

    if args.target_rungroup not in rungroups:
        raise RuntimeError(
            f"target rungroup not found in {config}: "
            f"{args.target_rungroup}"
        )

    if args.reference_rungroup not in rungroups:
        raise RuntimeError(
            f"reference rungroup not found in {config}: "
            f"{args.reference_rungroup}"
        )

    target = rungroups[args.target_rungroup]
    target_id = int(target["optics_id"])
    rootfile = Path(target["rootfile"])

    if not rootfile.is_file():
        raise RuntimeError(f"target ROOT file not found: {rootfile}")

    num_foils = load_num_foils(target_id)

    theta_tsv = (
        campaign
        / "02b_angle_scan_x"
        / "tsv"
        / "xfp_xpfp_selected_thetas.tsv"
    )

    if not theta_tsv.is_file():
        raise RuntimeError(f"X theta table not found: {theta_tsv}")

    ytar_cut_file = (
        campaign
        / "01_ytar_cuts"
        / "cuts"
        / f"ytar_ridge_cut_{args.target_rungroup}.root"
    )

    if not ytar_cut_file.is_file():
        raise RuntimeError(f"Ytar cut file not found: {ytar_cut_file}")

    theta = load_reference_thetas(
        theta_tsv,
        args.reference_rungroup,
        args.reference_foil,
    )

    print(f"Campaign:          {campaign}")
    print(f"Target rungroup:   {args.target_rungroup}")
    print(f"Target optics ID:  {target_id}")
    print(f"Target foils:      {num_foils}")
    print(f"Reference group:   {args.reference_rungroup}")
    print(f"Reference foil:    {args.reference_foil}")
    print(f"Theta table:       {theta_tsv}")
    print(f"Ytar cuts:         {ytar_cut_file}")

    for foil in range(num_foils):
        for ndel in range(5):
            delta_min = EDGES[ndel]
            delta_max = EDGES[ndel + 1]

            split, event_count = compute_split(
                rootfile,
                ytar_cut_file,
                foil,
                delta_min,
                delta_max,
            )

            if split is None:
                print(
                    f"SKIP target={args.target_rungroup} "
                    f"foil={foil} ndel={ndel}: "
                    f"only N={event_count}"
                )
                continue

            print(
                f"\nTARGET foil={foil} ndel={ndel} "
                f"delta=[{delta_min},{delta_max}) "
                f"split={split:.8g} N={event_count}"
            )

            zones = (
                ("low", -999.0, split),
                ("high", split, 999.0),
            )

            for zone, xpfp_min, xpfp_max in zones:
                fixed_theta = theta[(ndel, zone)]

                print(
                    f"  zone={zone} theta={fixed_theta} "
                    f"xpfp=[{xpfp_min},{xpfp_max}]"
                )

                expression = (
                    f'{MACRO}('
                    f'{target_id},{delta_min},{delta_max},'
                    f'"{args.target_rungroup}",'
                    f'9,1.0,0.12,0.06,0.18,2,'
                    f'0.003,0.025,0.30,'
                    f'true,{foil},-1,-1,-999,false,'
                    f'{fixed_theta},{xpfp_min},{xpfp_max},true,'
                    f'"{campaign}","{rootfile}")'
                )

                run_expression(expression, args.dry_run)


if __name__ == "__main__":
    main()
