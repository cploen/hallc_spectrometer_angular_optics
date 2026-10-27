#!/usr/bin/env python3
"""Isolated PyROOT I/O for continuous-delta joint-surface campaigns.

No numerical fit is run here. Replay inputs are read-only; canonical outputs
are exclusively published from temporary files. Candidate trees retain the
legacy GMM contract, while membership is evaluated before writing this file.
"""
import argparse
from array import array
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from spectrometer_config import from_campaign


def require_fresh(path):
    path = Path(path)
    if os.path.lexists(path):
        raise FileExistsError("Refusing to overwrite " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def export_replay(payload):
    """Reuse the established PID/foil/delta exporter, without X/Y cuts."""
    import ROOT
    # Load as a file: the isolated ROOT interpreter need not also have SciPy.
    module_spec = importlib.util.spec_from_file_location('joint_preband_export',
                                                        REPO / 'geometric_selection/export_preband.py')
    exporter = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(exporter)
    spec = from_campaign(payload['campaign'])
    required = {'xtar': 'gtr.x', 'xptar': 'gtr.th', 'yptar': 'gtr.ph'}
    optional = {'reactx': 'react.x', 'reacty': 'react.y', 'reactz': 'react.z',
                'xbpm_tar': 'rb.raster.fr_xbpm_tar', 'ybpm_tar': 'rb.raster.fr_ybpm_tar',
                'frx': 'rb.raster.fr_xa', 'fry': 'rb.raster.fr_ya'}
    source = ROOT.TFile.Open(payload['input'], 'READ')
    if not source or source.IsZombie():
        raise OSError('Cannot open replay ' + payload['input'])
    try:
        tree_names = ('T', 'Tout') if payload.get('tree', 'auto') == 'auto' else (payload['tree'],)
        tree = next((source.Get(n) for n in tree_names if source.Get(n)), None)
        if not tree:
            raise ValueError('Cannot find replay T or Tout')
        # The exporter remaps the suffix to the arm selected by folder name.
        additions = {k: 'H.' + v for k, v in required.items()}
        for key, suffix in optional.items():
            leaf = tree.GetLeaf(spec.branch(suffix))
            if leaf and leaf.GetLenStatic() == 1 and not leaf.GetLeafCount():
                additions[key] = 'H.' + suffix
    finally:
        source.Close()
    exporter.BRANCHES = dict(exporter.BRANCHES, **additions)
    args = exporter.parse_args([
        '--input', payload['input'], '--foil-cuts', payload['foil_cuts'],
        '--run-list', payload['metadata'], '--metadata-run', str(payload['optics_id']),
        '--campaign', payload['campaign'], '--rungroup', payload['rungroup'],
        '--cal-cut', str(payload.get('cal_cut', .65)),
        '--cer-cut', str(payload.get('cer_cut', 6.)), '--tree', payload.get('tree', 'auto'),
        '--output', payload['output'],
    ] + (['--hash-replay'] if payload.get('hash_replay') else []))
    exporter.export(args)


def write_candidates(payload):
    """Write independent X/Y pre-GMM populations from event-level labels."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    spec = from_campaign(payload['campaign'])
    with np.load(payload['events'], allow_pickle=False) as source:
        data = {k: source[k].copy() for k in source.files}
    if np.any(data['entry'] < 0) or np.any(data['entry'] > np.iinfo(np.int32).max):
        raise ValueError('Original replay entry exceeds the unchanged GMM mask signed-32-bit entry limit '
                         '(2147483647); candidate export stopped before writing. Split the replay input '
                         'or upgrade the downstream mask contract explicitly.')
    labels = {axis: np.full(len(data['entry']), -1, np.int32) for axis in ('x', 'y')}
    by_entry = {int(entry): i for i, entry in enumerate(data['entry'])}
    if len(by_entry) != len(data['entry']):
        raise ValueError('Duplicate original replay entries in pre-band events')
    seen_entries = set()
    for member_file in payload['memberships']:
        with np.load(member_file, allow_pickle=False) as member:
            entries = [int(e) for e in member['entry']]
            if seen_entries.intersection(entries) or len(set(entries)) != len(entries):
                raise ValueError('Repeated original entry in membership artifacts')
            seen_entries.update(entries)
            indices = np.asarray([by_entry[int(e)] for e in member['entry']], np.int64)
            for axis in ('x', 'y'):
                labels[axis][indices] = member[axis + '_label']
    if seen_entries != set(by_entry):
        raise ValueError('Membership artifacts do not cover every exported pre-band event')
    paths = []
    for axis in ('x', 'y'):
        paths.extend([payload['outputs'][axis]['root'], payload['outputs'][axis]['tsv']])
    for path in paths:
        require_fresh(path)
    double_fields = ('sumnpe', 'etracknorm', 'delta', 'ytar', 'xtar', 'xptar', 'yptar',
                     'xfp', 'xpfp', 'yfp', 'ypfp', 'xsieve', 'ysieve', 'reactx', 'reacty',
                     'reactz', 'xbpm_tar', 'ybpm_tar', 'frx', 'fry')
    reports = {}
    for axis in ('x', 'y'):
        selected = np.flatnonzero(labels[axis] >= 0)
        slots = spec.nx if axis == 'x' else spec.ny
        if np.any(labels[axis][selected] >= slots):
            raise ValueError('Physical ' + axis.upper() + ' label outside spectrometer grid')
        outputs = payload['outputs'][axis]
        destination = Path(outputs['root'])
        fd, temporary = tempfile.mkstemp(prefix='.joint-candidates-', suffix='.root', dir=destination.parent)
        os.close(fd)
        try:
            out = ROOT.TFile.Open(temporary, 'RECREATE')
            tree = ROOT.TTree('TXCand' if axis == 'x' else 'TYCand',
                             'Continuous-delta joint-surface pre-GMM candidates')
            ints = {key: array('i', [0]) for key in
                    ('run', 'foil', 'ndel', axis + 'scol', 'n' + axis.upper() + 'Matches',
                     'base_pass', 'ytar_cut_pass', axis + 'cut_pass')}
            if axis == 'x':
                ints['nXComponentMatches'] = array('i', [1])
            entry = array('q', [0])
            doubles = {key: array('d', [0.]) for key in double_fields + ('zfoil', 'delta_low', 'delta_high')}
            tree.Branch('entry', entry, 'entry/L')
            for name, value in ints.items():
                tree.Branch(name, value, name + '/I')
            for name, value in doubles.items():
                tree.Branch(name, value, name + '/D')
            counts = np.zeros((len(payload['foil_z_cm']), len(payload['delta_edges']) - 1, slots), np.int64)
            for i in selected:
                foil, ndel, label = int(data['foil'][i]), int(data['ndel'][i]), int(labels[axis][i])
                entry[0] = int(data['entry'][i])
                for name, value in ints.items():
                    value[0] = (int(payload['optics_id']) if name == 'run' else
                                foil if name == 'foil' else ndel if name == 'ndel' else
                                label if name == axis + 'scol' else 1)
                for name in double_fields:
                    doubles[name][0] = float(data[name][i]) if name in data else 0.
                doubles['zfoil'][0] = float(payload['foil_z_cm'][foil])
                doubles['delta_low'][0] = float(payload['delta_edges'][ndel])
                doubles['delta_high'][0] = float(payload['delta_edges'][ndel + 1])
                tree.Fill()
                counts[foil, ndel, label] += 1
            out.cd()
            tree.Write()
            ROOT.TNamed('hallc_spectrometer', spec.name).Write()
            ROOT.TNamed('hallc_geometry_profile', spec.name + '_centered_v1').Write()
            ROOT.TNamed('joint_surface_selection_manifest', json.dumps(payload)).Write()
            out.Close()
            # Atomic no-replace publication; never a partial canonical ROOT file.
            os.link(temporary, destination)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        with Path(outputs['tsv']).open('x', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(['rungroup', 'run', 'foil', 'ndel', axis + 'scol', 'N',
                             'delta_low', 'delta_high', 'zfoil'])
            for foil, zfoil in enumerate(payload['foil_z_cm']):
                for ndel, (low, high) in enumerate(zip(payload['delta_edges'], payload['delta_edges'][1:])):
                    for label in range(slots):
                        writer.writerow([payload['rungroup'], payload['optics_id'], foil, ndel, label,
                                         int(counts[foil, ndel, label]), low, high, zfoil])
        reports[axis] = {'rows': len(selected), 'unique_entries': len(np.unique(data['entry'][selected])),
                         'root': str(destination), 'tsv': outputs['tsv']}
    return reports


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['export', 'candidates'])
    parser.add_argument('--payload', type=Path, required=True)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args(argv)
    try:
        payload = json.loads(args.payload.read_text())
        report = export_replay(payload) if args.mode == 'export' else write_candidates(payload)
        if args.report:
            with require_fresh(args.report).open('x') as stream:
                json.dump(report, stream, indent=2)
        return 0
    except Exception as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
