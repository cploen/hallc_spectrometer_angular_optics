#!/usr/bin/env python3
"""Copy campaign inputs and downstream sample settings into a fresh campaign."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

from setup_geometric_campaign import create_campaign
from spectrometer_config import run_metadata

REPO = Path(__file__).resolve().parent


def create_joint_campaign(source, destination, metadata=None):
    source = source.expanduser().resolve()
    destination = destination.expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError('Refusing to overwrite existing destination: ' + str(destination))
    tables = sorted((source / 'config').glob('rungroups_*_inputs.tsv'))
    if len(tables) != 1:
        raise ValueError('Expected one rungroups_*_inputs.tsv in ' + str(source / 'config'))
    snapshot = source / 'config/optics.dat'
    metadata = (metadata if metadata is not None else
                snapshot if snapshot.exists() else REPO / 'DATfiles/list_of_optics_run.dat')
    metadata = metadata.expanduser().resolve()
    # Validate the handoff before creating any destination files. No replay I/O.
    with tables[0].open(newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            run_metadata(int(row['optics_id']), metadata)
    inherited = [source / 'config' / name for name in
                 ('core_sample.json', 'core_balance.json', 'sieve_mask.json')
                 if (source / 'config' / name).is_file()]
    marker = create_campaign(source, destination)
    config = destination / 'config'
    (config / 'geometric_selection.json').unlink()
    for path in inherited:
        shutil.copy2(path, config / path.name)
    shutil.copy2(metadata, config / 'optics.dat')
    marker['copied_configs'] += [path.name for path in inherited] + ['optics.dat']
    marker['source_sha256'].update({path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                   for path in inherited})
    marker['source_sha256']['optics.dat'] = hashlib.sha256(metadata.read_bytes()).hexdigest()
    marker['metadata_source'] = str(metadata)
    marker['selection_method'] = 'joint_surface_provisional'
    marker['notes'] = [
        'No numbered outputs or historical labels were copied.',
        'Fit inputs are measured FP coordinates and continuous delta only.',
        'Physical names are assigned after geometry is frozen, using one global orientation/origin per foil/axis.',
        'Core settings, allocation policy and sieve mask are inherited when present; metadata is snapshotted in config/optics.dat.',
        'The driver stops after conventional candidates; density cores and allocation are separate commands.',
    ]
    with (config / 'joint_surface_campaign.json').open('x') as stream:
        json.dump(marker, stream, indent=2)
        stream.write('\n')
    with (config / 'joint_surface.json').open('x') as stream:
        json.dump({'schema_version': 1,
                   'selection': {'cal_cut': .65, 'cer_cut': 6., 'tree': 'auto', 'hash_replay': False},
                   'coverage_min_events': 10, 'fit': {'common': {}, 'x': {}, 'y': {}},
                   'global_maps': {}}, stream, indent=2)
        stream.write('\n')
    return marker


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_campaign', type=Path)
    parser.add_argument('destination_campaign', type=Path)
    parser.add_argument('--metadata', type=Path,
                        help='Snapshot this run metadata; defaults to source config/optics.dat or repository metadata')
    args = parser.parse_args(argv)
    try:
        marker = create_joint_campaign(args.source_campaign, args.destination_campaign, args.metadata)
        print('Created ' + marker['spectrometer'] + ' joint-surface campaign: ' + str(args.destination_campaign))
        print('Run: python3 run_joint_surface_campaign.py ' + str(args.destination_campaign) + ' --plan')
        return 0
    except (OSError, ValueError) as exc:
        parser.exit(1, 'ERROR: ' + str(exc) + '\n')


if __name__ == '__main__':
    raise SystemExit(main())
