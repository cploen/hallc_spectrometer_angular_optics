#!/usr/bin/env python3
"""Copy only campaign inputs into a fresh joint-surface test campaign."""
import argparse
import json
from pathlib import Path

from setup_geometric_campaign import create_campaign


def create_joint_campaign(source, destination):
    marker = create_campaign(source, destination)
    destination = destination.expanduser().absolute()
    config = destination / 'config'
    (config / 'geometric_selection.json').unlink()
    marker['selection_method'] = 'joint_surface_provisional'
    marker['notes'] = [
        'No numbered outputs or historical labels were copied.',
        'Fit inputs are measured FP coordinates and continuous delta only.',
        'Physical names are assigned after geometry is frozen, using one global orientation/origin per foil/axis.',
        'The driver stops after conventional pre-GMM candidates; cleanup and SVD are separate decisions.',
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
    args = parser.parse_args(argv)
    try:
        marker = create_joint_campaign(args.source_campaign, args.destination_campaign)
        print('Created ' + marker['spectrometer'] + ' joint-surface campaign: ' + str(args.destination_campaign))
        print('Run: python3 run_joint_surface_campaign.py ' + str(args.destination_campaign) + ' --plan')
        return 0
    except (OSError, ValueError) as exc:
        parser.exit(1, 'ERROR: ' + str(exc) + '\n')


if __name__ == '__main__':
    raise SystemExit(main())
