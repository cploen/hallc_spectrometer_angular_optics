#!/usr/bin/env python3
"""Bounded ROOT check of cumulative geometric X/Y relabeling.

Two foils and two delta slices each contain two manually specified focal-plane
components with 40 events apiece. A first-slice manual merge must survive all
later relabel calls. No numerical discovery, GMM cleanup or optics fit is run.
"""
import argparse
import json
from pathlib import Path
import re
import sys
import tempfile

import numpy as np
import ROOT

from geometric_root_smoke import foil_cuts, run, write_replay

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from spectrometer_config import from_campaign

RUN = 99011
TAG = 'cumulative'
EDGES = [-8., -5., 0.]
PER_COMPONENT = 40


def points_for(component):
    xp, x = (-.01, -4.) if component == 0 else (.01, 4.)
    return [(xp-.001, x-.5), (xp+.001, x-.5),
            (xp+.001, x+.5), (xp-.001, x+.5), (xp-.001, x-.5)]


def write_components(campaign, axis):
    prefix, stage = ('XpFpXFp', '02b_angle_scan_x') if axis == 'x' else ('YpFpYFp', '02a_angle_scan_y')
    path = campaign/stage/'cuts'/(prefix+'_'+TAG+'_auto_band_cut.root')
    path.parent.mkdir(parents=True)
    output = ROOT.TFile(str(path), 'RECREATE')
    expected = {}
    for foil in range(2):
        for ndel in range(2):
            for component in range(2):
                if axis == 'x':
                    name = f'hXpFpXFp_cut_candidate_{component:02d}_zone_global_nfoil_{foil}_ndel_{ndel}'
                else:
                    name = f'hYpFpYFp_cut_yscol_{component}_nfoil_{foil}_ndel_{ndel}'
                points = points_for(component)
                cut = ROOT.TCutG(name, len(points))
                for i, point in enumerate(points):
                    cut.SetPoint(i, *point)
                cut.SetVarX(axis+'pfp'); cut.SetVarY(axis+'fp'); cut.Write()
                expected[(foil, ndel, component)] = np.asarray(points)
    output.Close()
    return expected


def inspect_components(campaign, axis, expected):
    if axis == 'x':
        path = campaign/'03b_relabel_x/cuts'/('XpFpXFp_'+TAG+'_auto_band_cut_xscolRelabeled.root')
        pattern = re.compile(r'^hXpFpXFp_cut_xscol_(\d+)_nfoil_(\d+)_ndel_(\d+)_part(\d+)_zone_global_src_(\d+)$')
    else:
        path = campaign/'03a_relabel_y/cuts'/('YpFpYFp_'+TAG+'_auto_band_cut_yscolRelabeled.root')
        pattern = re.compile(r'^hYpFpYFp_cut_yscol_(\d+)_nfoil_(\d+)_ndel_(\d+)_part(\d+)_src_(\d+)$')
    output = ROOT.TFile(str(path), 'READ')
    assert output and not output.IsZombie(), path
    names = [key.GetName() for key in output.GetListOfKeys()]
    component_names = [name for name in names if pattern.fullmatch(name)]
    assert len(component_names) == len(set(component_names)) == 8, component_names
    seen = {}
    for name in component_names:
        label, foil, ndel, part, source = map(int, pattern.fullmatch(name).groups())
        identity = (foil, ndel, source)
        assert identity not in seen, ('Duplicate source component', identity)
        expected_label = 3 if source == 0 or (foil, ndel) == (0, 0) else 5
        assert label == expected_label, (axis, identity, label, expected_label)
        cut = output.Get(name)
        actual = np.array([(cut.GetPointX(i), cut.GetPointY(i)) for i in range(cut.GetN())])
        assert np.array_equal(actual, expected[identity]), ('Changed polygon', axis, identity)
        seen[identity] = (label, part)
    assert set(seen) == set(expected), ('Lost foil/delta/source components', axis, seen)
    # Both first-slice pieces remain distinct after their common label is set,
    # including after the other three foil/slice calls rewrite the cut file.
    assert seen[(0, 0, 0)] == (3, 0)
    assert seen[(0, 0, 1)] == (3, 1)
    output.Close()
    return dict(components=len(seen), foil_delta_groups=4,
                first_slice_merged_parts=2, unchanged_polygon_vertices=True)


def check_arm(temp, arm, metadata, root_command):
    campaign = temp/(arm+'_geometric_cumulative')
    config = campaign/'config'; config.mkdir(parents=True)
    spec = from_campaign(campaign)
    marker = dict(selection_method='geometric_provisional', experimental=True,
                  spectrometer=arm)
    (config/'geometric_selection.json').write_text(json.dumps(marker))
    replay = temp/(arm+'_cumulative_replay.root')
    (config/'rungroups_test_inputs.tsv').write_text(
        'rungroup\toptics_id\tangle_deg\tfoils\tnruns\truns\trootfile\n' +
        f'{TAG}\t{RUN}\t12.5\t-3,3\t1\t{RUN}\t{replay}\n')
    rng = np.random.RandomState(1704)
    rows = []
    for foil in range(2):
        for ndel, (low, high) in enumerate(zip(EDGES, EDGES[1:])):
            for component, label in enumerate((3, 5)):
                xp, x = (-.01, -4.) if component == 0 else (.01, 4.)
                for _ in range(PER_COMPONENT):
                    rows.append({'cal.etottracknorm': .9, spec.cer+'.npeSum': 8.,
                                 'gtr.y': (-1. if foil == 0 else 1.)+float(rng.normal(0., .02)),
                                 'gtr.dp': float(rng.uniform(low+.1, high-.1)),
                                 'dc.xp_fp': xp+float(rng.normal(0., .0001)),
                                 'dc.x_fp': x+float(rng.normal(0., .05)),
                                 'dc.yp_fp': xp+float(rng.normal(0., .0001)),
                                 'dc.y_fp': x+float(rng.normal(0., .05)),
                                 'extcor.xsieve': spec.xs(label)+float(rng.normal(0., .005)),
                                 'extcor.ysieve': spec.ys(label)+float(rng.normal(0., .005))})
    write_replay(replay, arm, rows, tree_name='T' if arm == 'HMS' else 'Tout')
    foil_cuts(campaign/'01_ytar_cuts/cuts'/('ytar_ridge_cut_'+TAG+'.root'),
              [(-1.5, -.5), (.5, 1.5)])
    summary = dict(spectrometer=arm, replay_events=len(rows), per_component=PER_COMPONENT)
    for axis in ('x', 'y'):
        expected = write_components(campaign, axis)
        override = config/(axis+'_label_override.tsv')
        # The second first-slice component automatically maps to physical5;
        # forcing physical3 exercises the manual merge interface explicitly.
        entry = f'{RUN} 0 0 global 1 3 force cumulative_merge\n' if axis == 'x' else f'{RUN} 0 0 1 3 force cumulative_merge\n'
        override.write_text(entry)
        run([sys.executable, REPO/('run_relabel_'+axis+'_geometric_all.py'), campaign,
             '--metadata', metadata, '--root-command', root_command,
             '--override-file', override, '--min-events', 30], temp)
        summary[axis] = inspect_components(campaign, axis, expected)
    print('PASS:', arm, 'both axes preserve 8 components across 2 foils x 2 delta slices;',
          'first-slice manual merge retained', flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root-command', default='root')
    parser.add_argument('--keep-output', type=Path, help='Fresh directory for inspection; default is temporary')
    args = parser.parse_args()
    context = tempfile.TemporaryDirectory(prefix='hallc-geometric-cumulative-') if args.keep_output is None else None
    temp = Path(context.name) if context else args.keep_output.absolute()
    if not context:
        temp.mkdir(parents=True, exist_ok=False)
    try:
        metadata = temp/'metadata.dat'
        metadata.write_text(f'{RUN},cumulative,12.5,2,1,3,0\n-3,3\n-8,-5,0\n')
        summaries = [check_arm(temp, arm, metadata, args.root_command) for arm in ('HMS', 'SHMS')]
        (temp/'cumulative_relabel_summary.json').write_text(json.dumps(summaries, indent=2)+'\n')
    finally:
        if context:
            context.cleanup()


if __name__ == '__main__':
    main()
