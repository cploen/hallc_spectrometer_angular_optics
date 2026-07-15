#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

EDGES = [-10, -8, -5, 0, 5, 10]
MACRO = "assign_yfp_ypfp_angleScanBands.C"
OPTICS_DAT = Path("DATfiles/list_of_optics_run.dat")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply stored YFP/YPFP theta values from one reference foil to every foil in a target rungroup."
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
            f"expected exactly one rungroups_*_inputs.tsv under {campaign}/config; found {len(matches)}"
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
) -> dict[int, float]:
    selected: dict[int, float] = {}

    with theta_tsv.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["tag"] != reference_rungroup:
                continue
            if int(row["foil"]) != reference_foil:
                continue

            selected[int(row["ndel"])] = float(row["thetaDeg"])

    missing = [ndel for ndel in range(5) if ndel not in selected]
    if missing:
        raise RuntimeError(
            f"missing Y theta rows for reference={reference_rungroup}, "
            f"foil={reference_foil}, ndel={missing}"
        )

    return selected


def run_command(expression: str, dry_run: bool) -> None:
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
            f"target rungroup not found in {config}: {args.target_rungroup}"
        )

    if args.reference_rungroup not in rungroups:
        raise RuntimeError(
            f"reference rungroup not found in {config}: {args.reference_rungroup}"
        )

    target = rungroups[args.target_rungroup]
    target_id = int(target["optics_id"])
    rootfile = Path(target["rootfile"])

    if not rootfile.exists():
        raise RuntimeError(f"target ROOT file does not exist: {rootfile}")

    num_foils = load_num_foils(target_id)

    theta_tsv = (
        campaign
        / "02a_angle_scan_y"
        / "tsv"
        / "yfp_ypfp_selected_thetas.tsv"
    )

    if not theta_tsv.is_file():
        raise RuntimeError(f"Y theta table not found: {theta_tsv}")

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
    print()

    for foil in range(num_foils):
        for ndel in range(5):
            delta_min = EDGES[ndel]
            delta_max = EDGES[ndel + 1]
            fixed_theta = theta[ndel]

            print(
                f"RUN target={args.target_rungroup} "
                f"foil={foil} ndel={ndel} "
                f"delta=[{delta_min},{delta_max}) "
                f"theta={fixed_theta}"
            )

            expression = (
                f'{MACRO}('
                f'{target_id},{delta_min},{delta_max},'
                f'"{args.target_rungroup}",'
                f'9,1.0,0.18,0.08,0.25,1,0.005,0.08,1.0,'
                f'true,{foil},-1,-1,-999,false,{fixed_theta},'
                f'"{campaign}","{rootfile}")'
            )

            run_command(expression, args.dry_run)


if __name__ == "__main__":
    main()
