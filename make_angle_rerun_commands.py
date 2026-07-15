#!/usr/bin/env python3
import argparse
import csv
import os
import re
import shlex
import sys
from pathlib import Path


DELTA_EDGES = [
    (-10, -8),
    (-8, -5),
    (-5, 0),
    (0, 5),
    (5, 10),
]

XMACRO = "assign_xfp_xpfp_angleScanBands_split.C"
YMACRO = "assign_yfp_ypfp_angleScanBands.C"


def parse_encoded_number(value):
    return float(value.replace("m", "-").replace("p", "."))


def delta_index(delta_min, delta_max):
    pair = (int(delta_min), int(delta_max))
    if pair not in DELTA_EDGES:
        raise ValueError(f"unsupported delta interval {pair}")
    return DELTA_EDGES.index(pair)


def canonical_number(value):
    return f"{float(value):g}"


def load_campaign_rows(path):
    rows = {}

    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")

        required = {"rungroup", "optics_id", "rootfile"}
        missing = required.difference(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                f"{path} is missing required columns: {sorted(missing)}"
            )

        for row in reader:
            tag = row["rungroup"].strip()
            if not tag:
                continue

            rows[tag] = {
                "optics_id": int(row["optics_id"]),
                "rootfile": row["rootfile"].strip(),
            }

    return rows


def load_y_thetas(path):
    values = {}

    if not path.exists():
        return values

    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")

        for row in reader:
            key = (
                row["tag"].strip(),
                int(row["foil"]),
                int(row["ndel"]),
            )
            values[key] = row["thetaDeg"].strip()

    return values


def load_x_thetas(path):
    values = {}

    if not path.exists():
        return values

    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")

        for row in reader:
            key = (
                row["tag"].strip(),
                int(row["foil"]),
                int(row["ndel"]),
                row.get("zone", "").strip(),
                canonical_number(row.get("xpfpMin", "")),
                canonical_number(row.get("xpfpMax", "")),
            )
            values[key] = row["thetaDeg"].strip()

    return values


def shell_assignment(name, value):
    return f"{name}={shlex.quote(str(value))}"


def shell_root_command(expression):
    quoted = shlex.quote(expression)
    quoted = quoted.replace("__THETA__", '\'"$THETA"\'')
    quoted = quoted.replace("__ROOTFILE__", '\'"${ROOTFILE}"\'')
    return f"hcana -b -l -q {quoted}"


def find_open_evince_pdfs():
    """Return PDF paths appearing in command lines of running Evince processes."""
    pdfs = []

    for proc_dir in Path("/proc").iterdir():
        if not proc_dir.name.isdigit():
            continue

        try:
            if proc_dir.stat().st_uid != os.getuid():
                continue
            raw = (proc_dir / "cmdline").read_bytes()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue

        args = [
            item.decode(errors="replace")
            for item in raw.split(b"\0")
            if item
        ]

        if not args:
            continue

        executable = Path(args[0]).name.lower()
        if "evince" not in executable:
            continue

        for arg in args[1:]:
            if arg.lower().endswith(".pdf"):
                pdfs.append(str(Path(arg).expanduser().resolve()))

    return sorted(set(pdfs))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Emit single-line ROOT rerun commands for angle-scan PDFs. "
            "PDF paths may be supplied on standard input; otherwise the "
            "script automatically checks running Evince processes."
        )
    )
    parser.add_argument("campaign", help="campaign directory, e.g. HMS_5p878GeV")
    args = parser.parse_args()

    project_dir = Path.cwd().resolve()
    campaign_dir = project_dir / args.campaign
    config_dir = campaign_dir / "config"

    config_files = sorted(config_dir.glob("rungroups_*_inputs.tsv"))

    if len(config_files) != 1:
        raise SystemExit(
            f"ERROR: expected exactly one rungroups_*_inputs.tsv in {config_dir}"
        )

    rows = load_campaign_rows(config_files[0])

    xtheta_path = (
        campaign_dir
        / "02b_angle_scan_x"
        / "tsv"
        / "xfp_xpfp_selected_thetas.tsv"
    )

    ytheta_path = (
        campaign_dir
        / "02a_angle_scan_y"
        / "tsv"
        / "yfp_ypfp_selected_thetas.tsv"
    )

    xtheta = load_x_thetas(xtheta_path)
    ytheta = load_y_thetas(ytheta_path)

    if sys.stdin.isatty():
        pdfs = find_open_evince_pdfs()
        if not pdfs:
            raise SystemExit(
                "ERROR: no open Evince PDFs found and no PDF paths were supplied "
                "on standard input."
            )
        print(f"# Found {len(pdfs)} open angle-scan PDF(s) in Evince.")
    else:
        pdfs = sorted({
            line.strip()
            for line in sys.stdin
            if line.strip().lower().endswith(".pdf")
        })
        if not pdfs:
            raise SystemExit("ERROR: no PDF paths received on standard input.")

    for pdf in pdfs:
        base = Path(pdf).name
        print(f"# {pdf}")

        if base.startswith("xfp_xpfp_angleScan_"):
            match = re.fullmatch(
                r"xfp_xpfp_angleScan_(.+)_foil(\d+)_delta_"
                r"(m?\d+(?:p\d+)?)_to_(m?\d+(?:p\d+)?)_"
                r"xpfp_(m?\d+(?:p\d+)?)_to_(m?\d+(?:p\d+)?)\.pdf",
                base,
            )

            if not match:
                print("# could not parse X filename\n")
                continue

            tag, foil, dmin_s, dmax_s, xmin_s, xmax_s = match.groups()

            if tag not in rows:
                print(f"# rungroup not found in campaign TSV: {tag}\n")
                continue

            optics_id = rows[tag]["optics_id"]
            rootfile = rows[tag]["rootfile"]

            dmin = parse_encoded_number(dmin_s)
            dmax = parse_encoded_number(dmax_s)
            xmin = parse_encoded_number(xmin_s)
            xmax = parse_encoded_number(xmax_s)
            n_delta = delta_index(dmin, dmax)

            if xmin <= -998:
                zone = "low"
            elif xmax >= 998:
                zone = "high"
            else:
                zone = "mid"

            key = (
                tag,
                int(foil),
                n_delta,
                zone,
                canonical_number(xmin),
                canonical_number(xmax),
            )

            theta = xtheta.get(key, "EDIT_ME")

            expression_1 = (
                f'{XMACRO}({optics_id},{dmin:g},{dmax:g},"{tag}",'
                f'9,1.0,0.12,0.06,0.18,2,0.003,0.025,0.30,true,'
                f'{foil},-1,-1,-999,false,__THETA__,{xmin:g},{xmax:g},'
                f'true,"{args.campaign}","__ROOTFILE__")'
            )

            expression_2 = (
                f'{XMACRO}({optics_id},{dmin:g},{dmax:g},"{tag}",'
                f'9,1.0,0.10,0.045,0.12,2,0.0015,0.012,0.22,true,'
                f'{foil},-1,-1,-999,false,__THETA__,{xmin:g},{xmax:g},'
                f'true,"{args.campaign}","__ROOTFILE__")'
            )

            print(shell_assignment("THETA", theta))
            print(shell_assignment("ROOTFILE", rootfile))
            print(shell_root_command(expression_1))
            print(shell_root_command(expression_2))
            print()

        elif base.startswith("yfp_ypfp_angleScan_"):
            match = re.fullmatch(
                r"yfp_ypfp_angleScan_(.+)_foil(\d+)_delta_"
                r"(m?\d+(?:p\d+)?)_to_(m?\d+(?:p\d+)?)\.pdf",
                base,
            )

            if not match:
                print("# could not parse Y filename\n")
                continue

            tag, foil, dmin_s, dmax_s = match.groups()

            if tag not in rows:
                print(f"# rungroup not found in campaign TSV: {tag}\n")
                continue

            optics_id = rows[tag]["optics_id"]
            rootfile = rows[tag]["rootfile"]

            dmin = parse_encoded_number(dmin_s)
            dmax = parse_encoded_number(dmax_s)
            n_delta = delta_index(dmin, dmax)

            theta = ytheta.get(
                (tag, int(foil), n_delta),
                "EDIT_ME",
            )

            expression_1 = (
                f'{YMACRO}({optics_id},{dmin:g},{dmax:g},"{tag}",'
                f'9,1.0,0.18,0.08,0.25,2,0.005,0.04,0.45,true,'
                f'{foil},-1,-1,-999,false,__THETA__,'
                f'"{args.campaign}","__ROOTFILE__")'
            )

            expression_2 = (
                f'{YMACRO}({optics_id},{dmin:g},{dmax:g},"{tag}",'
                f'9,1.0,0.14,0.05,0.15,2,0.0025,0.02,0.30,true,'
                f'{foil},-1,-1,-999,false,__THETA__,'
                f'"{args.campaign}","__ROOTFILE__")'
            )

            print(shell_assignment("THETA", theta))
            print(shell_assignment("ROOTFILE", rootfile))
            print(shell_root_command(expression_1))
            print(shell_root_command(expression_2))
            print()

        else:
            print("# unrecognized angle-scan PDF\n")


if __name__ == "__main__":
    main()
