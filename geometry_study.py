#!/usr/bin/env python3
"""Run ROOT core-event geometry diagnostics using campaign paths and metadata."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from spectrometer_config import from_campaign, run_metadata

PROJECT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign')
    parser.add_argument('--config', type=Path)
    parser.add_argument('--foil', type=float, help='Physical foil position in cm')
    parser.add_argument('--rungroup', action='append', help='Exact name; repeat to select several')
    parser.add_argument('--check', action='store_true', help='Check inputs without running ROOT')
    parser.add_argument('--central-dense-half', action='store_true',
                        help='Central hole only; compare full core with highest-score half per delta slice')
    args = parser.parse_args()
    campaign = Path(args.campaign)
    if not campaign.is_dir():
        campaign = PROJECT / campaign
    campaign = campaign.resolve()
    spec = from_campaign(campaign)
    config = args.config or campaign / 'config/geometry_study.json'
    cfg = json.loads(config.read_text())
    foil = cfg.get('foil', 0.) if args.foil is None else args.foil
    if not math.isfinite(foil):
        raise ValueError('Foil must be finite')
    metadata = Path(os.environ.get('OPTICS_METADATA', PROJECT / 'DATfiles/list_of_optics_run.dat'))
    tables = list((campaign / 'config').glob('rungroups_*_inputs.tsv'))
    if len(tables) != 1:
        raise ValueError('Expected one campaign rungroup input table')
    with tables[0].open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if args.rungroup and set(args.rungroup) - {r['rungroup'] for r in rows}:
        raise ValueError('Unknown rungroup requested')
    def resolve(value):
        path = Path(value)
        if path.is_absolute():
            return path
        # Existing campaign tables can use repository-relative paths.
        return PROJECT / path if path.parts[0] == campaign.name else campaign / path
    jobs = []
    for row in rows:
        name = row['rungroup']
        if args.rungroup and name not in args.rungroup:
            continue
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', name):
            raise ValueError('Unsafe rungroup name')
        meta = run_metadata(int(row['optics_id']), metadata)
        if not any(abs(z - foil) < 1e-8 for z in meta['foils']):
            continue
        angle = float(row.get('angle_deg') or row.get('hms_angle_deg'))
        if abs(angle - meta['angle_deg']) > 1e-6:
            raise ValueError(f'{name}: angle disagrees with optics metadata')
        if spec.name == 'SHMS' and meta['sieve_flag'] != 1:
            raise ValueError('Only centered SHMS sieve is supported')
        core = resolve(cfg['core_file'].format(rungroup=name))
        replay = resolve(row['rootfile'])
        jobs.append((row, core, replay, angle))
    if not jobs:
        raise ValueError('No requested run groups contain this foil')
    offset = int(cfg.get('hole_step', 2))
    holes = cfg.get('holes') or [[spec.nx//2, spec.ny//2],
        [spec.nx//2-offset, spec.ny//2], [spec.nx//2+offset, spec.ny//2],
        [spec.nx//2, spec.ny//2-offset], [spec.nx//2, spec.ny//2+offset]]
    if args.central_dense_half:
        holes = [[spec.nx//2, spec.ny//2]]
    if any(len(h) != 2 or any(type(i) is not int for i in h) or
           not (0 <= h[0] < spec.nx and 0 <= h[1] < spec.ny) for h in holes):
        raise ValueError('Invalid [xscol, yscol] hole indices')
    minimum = cfg.get('min_events', 30)
    fraction = cfg.get('fp_fraction', .5)
    if type(minimum) is not int or minimum < 3 or not 0 < fraction < 1:
        raise ValueError('Require min_events >= 3 and 0 < fp_fraction < 1')
    print(f'{spec.name}; foil {foil:g} cm; holes {holes}', flush=True)
    print('Each run group and foil stays separate. All retained cores are used, including holdout.', flush=True)
    missing = []
    for row, core, replay, angle in jobs:
        print(f"{row['rungroup']}: angle {angle:g} deg\n  core: {core}\n  replay: {replay}", flush=True)
        if not core.is_file():
            missing.append(str(core))
        if not replay.is_file():
            print('  Replay unavailable: using saved replay values in CoreSample.', flush=True)
    if missing:
        raise FileNotFoundError('Missing inputs (no analysis run):\n' + '\n'.join(missing))
    if args.check:
        return
    root = shutil.which('root')
    if not root:
        raise RuntimeError('ROOT executable is required')
    out = resolve(cfg.get('output', '07_diagnostics/geometry')) / f'foil_{foil:g}cm'
    if args.central_dense_half:
        out = out / 'central_dense_half'
    out.mkdir(parents=True, exist_ok=True)
    manifest = dict(campaign=str(campaign), config=cfg, foil=foil, holes=holes,
                    central_dense_half=args.central_dense_half,
                    table=str(tables[0]), metadata=str(metadata), jobs=[], status='running')
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    for row, core, replay, angle in jobs:
        dest = out / row['rungroup']
        dest.mkdir(exist_ok=True)
        hole_string = ';'.join(f'{x},{y}' for x, y in holes)
        values = [str(core), str(replay), str(dest), spec.name, int(row['optics_id']),
                  angle, foil, hole_string, minimum, fraction, cfg.get('foil_width', 2.),
                  args.central_dense_half]
        expression = str(PROJECT / 'diagnostics/validation/geometry/core_geometry.C')
        expression += '(' + ','.join(json.dumps(v) for v in values) + ')'
        with (dest/'terminal.txt').open('w') as log:
            process = subprocess.Popen([root, '-l', '-b', '-q', expression],
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                print(line, end='', flush=True)
                log.write(line)
            if process.wait():
                raise RuntimeError(f"ROOT failed for {row['rungroup']}; see {dest/'terminal.txt'}")
        manifest['jobs'].append(dict(rungroup=row['rungroup'], angle=angle,
                                    core=str(core), replay=str(replay), source_replay_available=replay.is_file(),
                                    output=str(dest)))
    manifest['status'] = 'complete'
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, RuntimeError) as error:
        sys.exit(str(error))
