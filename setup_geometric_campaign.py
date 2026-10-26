#!/usr/bin/env python3
"""Create a fresh, explicitly provisional campaign for geometric X/Y selection.

Only the rungroup input table and seed matrix become active configuration.
Overrides referring to old provisional component IDs are copied as references.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys

from spectrometer_config import from_campaign


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def create_campaign(source: Path, destination: Path) -> dict:
    source = source.expanduser().resolve()
    destination = destination.expanduser().absolute()
    if not source.is_dir():
        raise ValueError(f"Source campaign does not exist: {source}")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to overwrite existing destination: {destination}")
    if not destination.parent.is_dir():
        raise ValueError(f"Destination parent does not exist: {destination.parent}")
    source_spec = from_campaign(source)
    destination_spec = from_campaign(destination)
    if source_spec.name != destination_spec.name:
        raise ValueError("Source and destination must use the same HMS/SHMS prefix")
    config = source / "config"
    tables = sorted(config.glob("rungroups_*_inputs.tsv"))
    if len(tables) != 1:
        raise ValueError(f"Expected one rungroups_*_inputs.tsv in {config}")
    table = tables[0]
    seed = config / "oldfit.dat"
    if not seed.is_file():
        raise FileNotFoundError(f"Missing seed optics matrix: {seed}")
    with table.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"rungroup", "optics_id", "foils", "nruns", "runs", "rootfile"}
        fields = reader.fieldnames
        if fields is None or not required.issubset(fields):
            raise ValueError(f"Missing required rungroup columns: {sorted(required-set(fields or []))}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty rungroup table: {table}")
    seen = set()
    for row in rows:
        tag = row["rungroup"]
        if not tag or tag in seen:
            raise ValueError(f"Empty or duplicate rungroup: {tag!r}")
        seen.add(tag)
        if not row["optics_id"].isdigit() or not row["rootfile"]:
            raise ValueError(f"Invalid optics_id or rootfile for {tag}")
        # Existing runners interpret relative replay paths from the repository.
        replay = Path(row["rootfile"]).expanduser()
        if not replay.is_absolute():
            replay = source.parent / replay
        row["rootfile"] = str(replay.resolve())
    references = sorted(p for p in config.iterdir() if p.is_file() and
                        ("override" in p.name.lower() or "veto" in p.name.lower()))
    marker = {
        "schema_version": 1,
        "selection_method": "geometric_provisional",
        "experimental": True,
        "spectrometer": source_spec.name,
        "source_campaign": str(source),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "copied_configs": [table.name, seed.name],
        "source_sha256": {table.name: _sha256(table), seed.name: _sha256(seed)},
        "reference_overrides": [p.name for p in references],
        "notes": [
            "No numbered stage outputs were copied; begin with step 1.",
            "Old override/veto component IDs are reference material, not active settings.",
            "The new campaign uses the repository DATfiles/list_of_optics_run.dat metadata.",
        ],
    }
    # mkdir is exclusive: an existing campaign is never reused or overwritten.
    destination.mkdir()
    try:
        new_config = destination / "config"
        new_config.mkdir()
        with (new_config / table.name).open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        shutil.copy2(seed, new_config / seed.name)
        if references:
            reference_dir = new_config / "previous_selection_reference"
            reference_dir.mkdir()
            for path in references:
                shutil.copy2(path, reference_dir / path.name)
        with (new_config / "geometric_selection.json").open("x") as stream:
            json.dump(marker, stream, indent=2)
            stream.write("\n")
    except Exception:
        # Leave a partial fresh folder for inspection; never remove user data.
        raise
    return marker


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_campaign", type=Path)
    parser.add_argument("destination_campaign", type=Path,
                        help="New path retaining the source HMS_ or SHMS_ prefix")
    args = parser.parse_args(argv)
    try:
        marker = create_campaign(args.source_campaign, args.destination_campaign)
    except (OSError, ValueError, csv.Error) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    destination = args.destination_campaign.expanduser().absolute()
    print(f"Created provisional {marker['spectrometer']} campaign: {destination}")
    print("Copied input table and oldfit.dat; no stage outputs were copied.")
    if marker["reference_overrides"]:
        print("Previous overrides/vetoes are in config/previous_selection_reference for review.")
    print("Begin with run_ytar_ridge_all_rungroups.sh, then the geometric X/Y selector.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
