#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import subprocess
import os
import socket
import time
from pathlib import Path
from spectrometer_config import from_campaign, run_metadata

import ROOT

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
    parser.add_argument("--verbose", action="store_true",
                        help="Report host/PID and split-reading progress every five seconds.")
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
    n_slices: int,
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
        for ndel in range(n_slices)
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


def compute_splits(
    rootfile: Path,
    ytar_cut_file: Path,
    num_foils: int,
    edges: list[float],
    spec,
    verbose: bool = False,
) -> dict[tuple[int, int], tuple[float | None, int]]:
    """Collect all foil/slice samples in one pass, preserving the split definition."""
    started = time.monotonic()
    if verbose:
        print(f"SPLIT START all foils/slices opening {rootfile}", flush=True)
    input_file = ROOT.TFile.Open(str(rootfile), "READ")
    if not input_file or input_file.IsZombie():
        raise RuntimeError(f"cannot open {rootfile}")

    cut_file = None
    try:
        tree = input_file.Get("T") or input_file.Get("Tout")
        if not tree:
            raise RuntimeError(f"cannot find T or Tout in {rootfile}")

        cut_file = ROOT.TFile.Open(str(ytar_cut_file), "READ")
        if not cut_file or cut_file.IsZombie():
            raise RuntimeError(f"cannot open {ytar_cut_file}")

        cuts = []
        for foil in range(num_foils):
            name = f"delta_vs_ytar_cut_foil{foil}"
            cut = cut_file.Get(name)
            if not cut:
                raise RuntimeError(f"cannot find {name} in {ytar_cut_file}")
            cuts.append(cut)

        branches = (
            spec.cherenkov_branch,
            spec.branch("cal.etottracknorm"),
            spec.branch("gtr.dp"),
            spec.branch("gtr.y"),
            spec.branch("dc.xp_fp"),
        )
        tree.SetBranchStatus("*", 0)
        for name in branches:
            if not tree.GetBranch(name):
                raise RuntimeError(f"missing required branch {name} in {rootfile}")
            tree.SetBranchStatus(name, 1)

        intervals = list(zip(edges[:-1], edges[1:]))
        samples = {
            (foil, ndel): []
            for foil in range(num_foils)
            for ndel in range(len(intervals))
        }
        total = tree.GetEntries()
        selected = 0
        last_report = time.monotonic()
        if verbose:
            print(f"SPLIT READ total={total} foils={num_foils} "
                  f"slices={len(intervals)} branches={len(branches)} passes=1",
                  flush=True)
        for entry in range(total):
            if verbose and time.monotonic() - last_report >= 5.0:
                elapsed = time.monotonic() - started
                print(f"SPLIT PROGRESS all foils/slices read={entry}/{total} "
                      f"selected_assignments={selected} elapsed={elapsed:.1f}s "
                      f"rate={entry / elapsed:.0f} events/s", flush=True)
                last_report = time.monotonic()
            if tree.GetEntry(entry) <= 0:
                raise RuntimeError(f"cannot read entry {entry} in {rootfile}")

            cer = float(getattr(tree, branches[0]))
            cal = float(getattr(tree, branches[1]))
            delta = float(getattr(tree, branches[2]))
            ytar = float(getattr(tree, branches[3]))
            xpfp = float(getattr(tree, branches[4]))

            if spec.name == "SHMS" and not spec.delta_min < delta < spec.delta_max:
                continue
            if cer <= 2.0 or cal <= 0.65:
                continue

            # Preserve half-open intervals and independent foil membership.
            # Do not assign an event to just one foil if supplied cuts overlap.
            matching_slices = [
                ndel for ndel, (lo, hi) in enumerate(intervals)
                if lo <= delta < hi
            ]
            if not matching_slices:
                continue
            for foil, cut in enumerate(cuts):
                if not cut.IsInside(ytar, delta):
                    continue
                for ndel in matching_slices:
                    samples[(foil, ndel)].append(xpfp)
                    selected += 1

        if verbose:
            print(f"SPLIT DONE read={total}/{total} "
                  f"selected_assignments={selected} "
                  f"elapsed={time.monotonic() - started:.1f}s", flush=True)
    finally:
        if cut_file:
            cut_file.Close()
        input_file.Close()

    splits = {}
    for key, values in samples.items():
        count = len(values)
        split = None
        if count >= 1000:
            q05 = quantile(values, 0.05)
            q95 = quantile(values, 0.95)
            if q05 is not None and q95 is not None:
                split = 0.5 * (q05 + q95)
        splits[key] = (split, count)
    return splits


def run_expression(expression: str, dry_run: bool) -> None:
    command = ["hcana", "-b", "-l", "-q", expression]

    if dry_run:
        print(subprocess.list2cmdline(command))
        return

    subprocess.run(command, check=True)


def main() -> None:
    args = parse_args()
    if args.verbose:
        print(f"PROCESS host={socket.gethostname()} pid={os.getpid()}", flush=True)

    campaign = Path(args.campaign)
    spec = from_campaign(campaign)
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
    target_meta=run_metadata(target_id)
    reference_meta=run_metadata(int(rungroups[args.reference_rungroup]["optics_id"]))
    edges=target_meta["edges"]
    if edges!=reference_meta["edges"]:
        raise ValueError("Target and reference delta boundaries differ; saved angles cannot be reused by index")
    if spec.name=="SHMS" and (target_meta["sieve_flag"]!=1 or reference_meta["sieve_flag"]!=1):
        raise ValueError("This version assumes the centered SHMS sieve")
    n_slices=len(edges)-1
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
        n_slices,
    )

    print(f"Campaign:          {campaign}")
    print(f"Target rungroup:   {args.target_rungroup}")
    print(f"Target optics ID:  {target_id}")
    print(f"Target foils:      {num_foils}")
    print(f"Reference group:   {args.reference_rungroup}")
    print(f"Reference foil:    {args.reference_foil}")
    print(f"Theta table:       {theta_tsv}")
    print(f"Ytar cuts:         {ytar_cut_file}")

    splits = compute_splits(
        rootfile, ytar_cut_file, num_foils, edges, spec, verbose=args.verbose
    )

    for foil in range(num_foils):
        for ndel in range(n_slices):
            delta_min = edges[ndel]
            delta_max = edges[ndel + 1]

            split, event_count = splits[(foil, ndel)]

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
                    f'{spec.nx},1.0,0.12,0.06,0.18,2,'
                    f'0.003,0.025,0.30,'
                    f'true,{foil},-1,-1,-999,false,'
                    f'{fixed_theta},{xpfp_min},{xpfp_max},true,'
                    f'"{campaign}","{rootfile}")'
                )

                run_expression(expression, args.dry_run)


if __name__ == "__main__":
    main()
