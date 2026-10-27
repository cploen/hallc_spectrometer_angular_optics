#!/usr/bin/env python3
"""Synthetic raw-replay I/O test for both arms and the unchanged GMM reader.

This checks software contracts, source entries, cuts, types and axis-independent
populations; it is not a replay physics validation. Numerical discovery has its
own independent tests. Run with a Python containing ROOT and NumPy.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import ROOT

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
module_spec = importlib.util.spec_from_file_location('geometric_fixture', REPO / 'tests/geometric_root_smoke.py')
fixture = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(fixture)
io_spec = importlib.util.spec_from_file_location('joint_root_io', REPO / 'joint_surface/root_io.py')
root_io = importlib.util.module_from_spec(io_spec)
io_spec.loader.exec_module(root_io)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gmm-python', help='Optional combined PyROOT/scikit-learn interpreter to run actual cleanup')
    args = parser.parse_args(argv)
    ROOT.gROOT.SetBatch(True)
    with tempfile.TemporaryDirectory(prefix='joint-surface-root-') as temporary:
        directory = Path(temporary)
        metadata = directory / 'metadata.dat'
        metadata.write_text('99002,test,12.5,2,1,3,0\n-3,3\n-10,0,22\n')
        cuts = directory / 'foil_cuts.root'
        fixture.foil_cuts(cuts, [(-1., 1.), (0., 2.)])
        for arm in ('HMS', 'SHMS'):
            spec = fixture.from_campaign(arm + '_fixture')
            for tree_name in ('T', 'Tout'):
                campaign = directory / (arm + '_' + tree_name)
                campaign.mkdir()
                conditions = [(-5., .5, 7., .9), (-10., -.5, 7., .9),
                              (0., -.5, 7., .9), (-5., -.5, 6., .9),
                              (-5., -.5, 7., .65), (16., .5, 7., .9),
                              (spec.delta_max, .5, 7., .9)]
                repetitions = 100 if args.gmm_python else 1
                rows = [dict(zip(['gtr.dp', 'gtr.y', spec.cer + '.npeSum', 'cal.etottracknorm'], values))
                        for repetition in range(repetitions) for values in conditions]
                rng = np.random.RandomState(1234)
                for row in rows:
                    row.update({'dc.x_fp': 3., 'dc.xp_fp': .02, 'dc.y_fp': 4., 'dc.yp_fp': .03,
                                'gtr.x': .1, 'gtr.th': .001, 'gtr.ph': .002,
                                'extcor.xsieve': float(2.54 + rng.normal(0, .025)),
                                'extcor.ysieve': float(1.524 + rng.normal(0, .025))})
                replay = directory / (arm + '_' + tree_name + '.root')
                fixture.write_replay(replay, arm, rows, tree_name)
                digest = hashlib.sha256(replay.read_bytes()).hexdigest()
                output = campaign / 'preband'
                payload = {'campaign': str(campaign), 'input': str(replay), 'foil_cuts': str(cuts),
                           'metadata': str(metadata), 'optics_id': 99002, 'rungroup': 'test', 'output': str(output)}
                root_io.export_replay(payload)
                with np.load(output / 'preband_events.npz') as source:
                    data = {k: source[k].copy() for k in source.files}
                expected_base = [0, 2, 5] if arm == 'SHMS' else [0, 2]
                expected = [7 * repetition + value for repetition in range(repetitions) for value in expected_base]
                assert data['entry'].tolist() == expected
                assert data['foil'].tolist() == ([1, 0, 1] if arm == 'SHMS' else [1, 0]) * repetitions
                assert data['ndel'].tolist() == ([0, 1, 1] if arm == 'SHMS' else [0, 1]) * repetitions
                assert data['xptar'].tolist() == [.001] * len(expected)
                assert data['ypfp'].tolist() == [.03] * len(expected)
                # Deliberately reject distinct events in the two planes.
                x_label, y_label = np.full(len(expected), 4, np.int32), np.full(len(expected), 3, np.int32)
                x_label[0] = -1
                y_label[-1] = -1
                membership = campaign / 'membership.npz'
                np.savez_compressed(membership, entry=data['entry'], x_label=x_label, y_label=y_label)
                outputs = {}
                for axis, folder, prefix in (('x', '04b_candidate_trees_x', 'Xscol'), ('y', '04a_candidate_trees_y', 'Yscol')):
                    outputs[axis] = {'root': str(campaign / folder / 'root' / (prefix + 'Candidates_test.root')),
                                     'tsv': str(campaign / folder / 'tsv' / (prefix + 'Candidates_test_summary.tsv'))}
                candidate_payload = {'campaign': str(campaign), 'rungroup': 'test', 'optics_id': 99002,
                                     'events': str(output / 'preband_events.npz'), 'memberships': [str(membership)],
                                     'foil_z_cm': [-3., 3.], 'delta_edges': [-10., 0., 22.], 'outputs': outputs}
                root_io.write_candidates(candidate_payload)
                for axis in ('x', 'y'):
                    file = ROOT.TFile.Open(outputs[axis]['root'], 'READ')
                    tree = file.Get('TXCand' if axis == 'x' else 'TYCand')
                    labels = x_label if axis == 'x' else y_label
                    assert [int(e.entry) for e in tree] == data['entry'][labels >= 0].tolist()
                    assert tree.GetLeaf('entry').GetTypeName() == 'Long64_t'
                    assert tree.GetLeaf(axis + 'scol').GetTypeName() == 'Int_t'
                    assert tree.GetLeaf('xpfp').GetTypeName() == 'Double_t'
                    for event in tree:
                        assert event.xpfp == .02
                        assert getattr(event, 'n' + axis.upper() + 'Matches') == 1
                    file.Close()
                    with Path(outputs[axis]['tsv']).open() as stream:
                        summary = list(csv.DictReader(stream, delimiter='\t'))
                    assert sum(int(row['N']) for row in summary) == int((labels >= 0).sum())
                    assert len(summary) == 2 * 2 * (spec.nx if axis == 'x' else spec.ny)
                    if args.gmm_python:
                        script = REPO / ('gmm_cleanup_' + axis + 'scol_candidates_v2.py')
                        selected = np.flatnonzero(labels >= 0)[0]
                        command = [args.gmm_python, str(script), '--campaign', str(campaign), '--rungroup', 'test',
                                   '--foil', str(data['foil'][selected]), '--ndel', str(data['ndel'][selected]),
                                   '--min-events', '2', '--keep-frac', '1', '--n-init', '1']
                        result = subprocess.run(command, cwd=REPO, text=True, capture_output=True)
                        assert result.returncode == 0, result.stdout + result.stderr
                        mask_folder = campaign / ('05b_gmm_cleanup_x' if axis == 'x' else '05a_gmm_cleanup_y') / 'root'
                        mask_path = next(mask_folder.glob('gmm_clean_*.root'))
                        mask_file = ROOT.TFile.Open(str(mask_path), 'READ')
                        mask = mask_file.Get('GMMClean')
                        expected_mask = data['entry'][(labels >= 0) & (data['foil'] == data['foil'][selected]) &
                                                      (data['ndel'] == data['ndel'][selected])].tolist()
                        assert sorted(int(e.entry) for e in mask) == sorted(expected_mask)
                        assert all(int(e.gmm_keep) == 1 for e in mask)
                        mask_file.Close()
                        print('PASS unchanged GMM reader, fit and exact keep-all mask:', arm, tree_name, axis)
                try:
                    root_io.write_candidates(candidate_payload)
                except FileExistsError:
                    pass
                else:
                    raise AssertionError('Candidate outputs were overwritten')
                assert hashlib.sha256(replay.read_bytes()).hexdigest() == digest
                # The legacy GMM mask stores entry/I. Never silently truncate.
                overflow_events = campaign / 'overflow_events.npz'
                overflow_members = campaign / 'overflow_members.npz'
                altered = dict(data)
                altered['entry'] = data['entry'].copy()
                altered['entry'][0] = 2147483648
                np.savez_compressed(overflow_events, **altered)
                np.savez_compressed(overflow_members, entry=altered['entry'], x_label=x_label, y_label=y_label)
                overflow_payload = dict(candidate_payload, events=str(overflow_events), memberships=[str(overflow_members)])
                try:
                    root_io.write_candidates(overflow_payload)
                except ValueError as exc:
                    assert '2147483647' in str(exc)
                else:
                    raise AssertionError('A downstream entry overflow was permitted')
                print('PASS', arm, tree_name, 'PID/foil/delta, original entries, candidate schema and counts')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
