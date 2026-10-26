#!/usr/bin/env python3
"""Synthetic ROOT campaign integration; never a calibration-quality assertion.

Run with PyROOT, passing --numerical-python when SciPy is in another interpreter.
The legacy GMM command requires PyROOT and scikit-learn in the *same* Python.
Without --gmm-python, compatible synthetic masks explicitly replace that stage.
"""
from array import array
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import ROOT

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from spectrometer_config import from_campaign

ROOT.gROOT.SetBatch(True)


def run(command, cwd, timeout=120):
    result = subprocess.run([str(x) for x in command], cwd=cwd, text=True,
                            capture_output=True, timeout=timeout,
                            env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1',
                                     MPLBACKEND='Agg',
                                     MPLCONFIGDIR='/private/tmp/hms-clustering-mpl',
                                     OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'))
    with (Path(cwd)/'software_smoke_commands.log').open('a') as log:
        log.write('\nCOMMAND: '+str([str(x) for x in command])+'\n')
        log.write(result.stdout+result.stderr+'\nEXIT: '+str(result.returncode)+'\n')
    if result.returncode:
        raise AssertionError(' '.join(map(str, command))+'\n'+result.stdout+result.stderr)
    return result


def macro(name, args, cwd, allow_diagnostic_error=False):
    arguments = ','.join('kFALSE' if x is False else 'kTRUE' if x is True
                         else json.dumps(x) for x in args)
    result = run(['root', '-l', '-b', '-q', '-e',
                  '#include '+json.dumps(str(REPO/name)), '-e',
                  Path(name).stem+'('+arguments+');'], cwd)
    if not allow_diagnostic_error:
        assert 'error:' not in result.stderr.lower(), result.stdout+result.stderr
    return result


def foil_cuts(path, intervals):
    path.parent.mkdir(parents=True, exist_ok=True)
    f = ROOT.TFile(str(path), 'RECREATE')
    for foil, (lo, hi) in enumerate(intervals):
        points = [(lo, -30.), (hi, -30.), (hi, 30.), (lo, 30.), (lo, -30.)]
        cut = ROOT.TCutG('delta_vs_ytar_cut_foil%d' % foil, len(points))
        for i, point in enumerate(points):
            cut.SetPoint(i, *point)
        cut.SetVarX('ytar'); cut.SetVarY('delta'); cut.Write()
    f.Close()


def write_replay(path, arm, rows, tree_name='T'):
    spec = from_campaign(arm+'_fixture')
    f = ROOT.TFile(str(path), 'RECREATE')
    tree = ROOT.TTree(tree_name, 'Synthetic preband replay')
    fields = ['cal.etottracknorm', spec.cer+'.npeSum', 'gtr.dp', 'gtr.x',
              'gtr.y', 'gtr.th', 'gtr.ph', 'dc.x_fp', 'dc.xp_fp', 'dc.y_fp',
              'dc.yp_fp', 'react.x', 'react.y', 'react.z', 'extcor.xsieve',
              'extcor.ysieve', 'rb.raster.fr_xbpm_tar', 'rb.raster.fr_ybpm_tar']
    values = {key: array('d', [0.]) for key in fields}
    for key, value in values.items():
        name = spec.branch(key); tree.Branch(name, value, name+'/D')
    for row in rows:
        for key, value in values.items():
            value[0] = row.get(key, 0.)
        tree.Fill()
    tree.Write(); f.Close()


def check_export(numerical_python, root_python, temp):
    """Replay thresholds, original entries, half-open slices, and foil precedence."""
    module_spec = importlib.util.spec_from_file_location('preband_export_fixture',
                                                        REPO/'geometric_selection/export_preband.py')
    exporter = importlib.util.module_from_spec(module_spec); module_spec.loader.exec_module(exporter)
    parse_args, export = exporter.parse_args, exporter.export
    metadata = temp/'export_metadata.dat'
    metadata.write_text('99002,test,12.5,2,1,3,0\n-3,3\n-10,0,22\n')
    cuts = temp/'export_foils.root'; foil_cuts(cuts, [(-1., 1.), (0., 2.)])
    for arm in ('HMS', 'SHMS'):
        spec = from_campaign(arm+'_fixture')
        # good, exact lower acceptance boundary, exact slice boundary,
        # Cherenkov threshold, calorimeter threshold, high-acceptance point.
        conditions = [(-5., .5, 7., .9), (-10., -.5, 7., .9),
                      (0., -.5, 7., .9), (-5., -.5, 6., .9),
                      (-5., -.5, 7., .65), (16., .5, 7., .9),
                      (spec.delta_max, .5, 7., .9)]
        rows = [dict(zip(['gtr.dp', 'gtr.y', spec.cer+'.npeSum',
                          'cal.etottracknorm'], values)) for values in conditions]
        for tree in ('T', 'Tout'):
            source = temp/(arm+'_'+tree+'.root'); write_replay(source, arm, rows, tree)
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            output = temp/(arm+'_'+tree+'_export')
            args = parse_args(['--input', str(source), '--foil-cuts', str(cuts),
                               '--run-list', str(metadata), '--metadata-run', '99002',
                               '--campaign', arm+'_fixture', '--output', str(output)])
            export(args)
            candidate = next(output.glob('*.npz'))
            with np.load(candidate) as data:
                expected = [0, 2, 5] if arm == 'SHMS' else [0, 2]
                assert data['entry'].tolist() == expected, (arm, tree, data['entry'])
                assert data['foil'].tolist() == ([1, 0, 1] if arm == 'SHMS' else [1, 0])
                assert data['ndel'].tolist() == ([0, 1, 1] if arm == 'SHMS' else [0, 1])
            try:
                export(args)
            except FileExistsError:
                pass
            else:
                raise AssertionError('Exporter reused existing output')
            assert hashlib.sha256(source.read_bytes()).hexdigest() == digest


def synthetic_masks(campaign, axis, rows):
    """Deliberate mask fixture when the combined legacy GMM environment is absent."""
    folder = campaign/('05a_gmm_cleanup_y' if axis == 'y' else '05b_gmm_cleanup_x')/'root'
    folder.mkdir(parents=True)
    path = folder/('gmm_clean_test_foil0_delta_m5_to_0_all'+axis.upper()+'.root')
    f = ROOT.TFile(str(path), 'RECREATE'); tree = ROOT.TTree('GMMClean', 'Software mask fixture')
    values = {name: array('i', [0]) for name in ('entry', 'foil', 'ndel', axis+'scol', 'gmm_keep')}
    for name, value in values.items(): tree.Branch(name, value, name+'/I')
    score = array('d', [1.]); tree.Branch('gmm_score', score, 'gmm_score/D')
    for row in rows:
        for name, value in values.items(): value[0] = 1 if name == 'gmm_keep' else row[name]
        tree.Fill()
    tree.Write(); ROOT.TNamed('hallc_spectrometer', from_campaign(campaign).name).Write(); f.Close()


def check_ridge(temp):
    """Bounded stage-1 check, separate from the fixed downstream cut fixture."""
    directory = temp/'ridge_fixture'; directory.mkdir()
    dat = directory/'DATfiles'; dat.mkdir()
    (dat/'list_of_optics_run.dat').write_text('99003,ridge,12.5,1,1,2,0\n0\n-10,22\n')
    rng = np.random.RandomState(984)
    for arm in ('HMS', 'SHMS'):
        spec = from_campaign(arm+'_ridge')
        campaign = directory/(arm+'_ridge')
        delta = rng.uniform(spec.delta_min+.1, spec.delta_max-.1, 6000)
        ytar = .012*delta+rng.normal(0., .12, len(delta))
        rows = [{'gtr.dp': float(dp), 'gtr.y': float(y), spec.cer+'.npeSum': 8.,
                 'cal.etottracknorm': .9} for dp, y in zip(delta, ytar)]
        replay = directory/(arm+'_ridge.root')
        write_replay(replay, arm, rows, 'T' if arm == 'HMS' else 'Tout')
        macro('ytar_ridge_cut.C', [99003, 'ridge', '-1', '-1', str(replay), str(campaign)], directory)
        output = campaign/'01_ytar_cuts/cuts/ytar_ridge_cut_ridge.root'
        assert output.is_file(), 'Stage 1 did not publish a foil cut'
        f = ROOT.TFile(str(output)); cut = f.Get('delta_vs_ytar_cut_foil0')
        assert cut and cut.GetN() > 3
        assert sum(bool(cut.IsInside(float(y), float(dp))) for y, dp in zip(ytar, delta)) > 3000
        f.Close()
    print('PASS: actual stage-1 ridge macro generated HMS/T and SHMS/Tout foil cuts')


def inspect_candidates(campaign, axis, expected):
    prefix = 'Yscol' if axis == 'y' else 'Xscol'
    stage = '04a_candidate_trees_y' if axis == 'y' else '04b_candidate_trees_x'
    f = ROOT.TFile(str(campaign/stage/'root'/(prefix+'Candidates_test.root')))
    assert f and not f.IsZombie()
    tree = f.Get('TYCand' if axis == 'y' else 'TXCand'); assert tree and tree.GetEntries() > 0
    rows = []
    for event in tree:
        entry = int(event.entry); label = int(getattr(event, axis+'scol'))
        assert label == int(expected[entry]), (axis, entry, label, expected[entry])
        assert int(getattr(event, 'nYMatches' if axis == 'y' else 'nXMatches')) == 1
        rows.append(dict(entry=entry, foil=int(event.foil), ndel=int(event.ndel), **{axis+'scol': label}))
    assert len({row['entry'] for row in rows}) == len(rows)
    f.Close(); return rows


def check_manual_preview(args, temp, campaign, metadata):
    """Round trip reviewed ROOT cuts without replacing canonical selections."""
    canonical = campaign/'02b_angle_scan_x/cuts/XpFpXFp_test_auto_band_cut.root'
    original = hashlib.sha256(canonical.read_bytes()).hexdigest()
    polygons = campaign/'config/reviewed_x.json'
    run([sys.executable, REPO/'geometric_selection/root_io.py', '--extract-cuts', canonical,
         '--axis', 'x', '--foil', 0, '--ndel', 0, '--output', polygons], temp)
    assert len(json.loads(polygons.read_text())) == from_campaign(campaign).nx
    override = campaign/'config/manual_preview.json'
    override.write_text(json.dumps(dict(slices={'test/foil0/delta0':
                                               dict(x=dict(polygons_file=str(polygons)))})))
    run([args.numerical_python, REPO/'run_geometric_selection.py', campaign,
         '--rungroup', 'test', '--axis', 'x', '--metadata', metadata,
         '--root-python', sys.executable, '--overrides', override,
         '--preview', '--foil', 0, '--ndel', 0], temp)
    assert hashlib.sha256(canonical.read_bytes()).hexdigest() == original
    print('PASS: ROOT cut extraction/manual preview retained canonical cut bytes')


def check_campaign(args, temp):
    summaries = []
    metadata = temp/'campaign_metadata.dat'
    metadata.write_text('99001,test,12.5,1,1,2,0\n0\n-5,0\n')
    dat = temp/'DATfiles'; dat.mkdir()
    (dat/'list_of_optics_run.dat').write_text(metadata.read_text())
    for arm in ('HMS', 'SHMS'):
        spec = from_campaign(arm+'_fixture')
        fixture = temp/(arm+'_arrays.npz')
        generator = '''
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(sys.argv[1])/'tests'))
from test_geometric_selection import x_fixture, y_fixture
x, xl, _ = x_fixture(slots=int(sys.argv[3]), per_band=550)
y, yl = y_fixture(per_band=len(x)//5)
order = np.random.RandomState(613).permutation(len(y)); y, yl = y[order], yl[order]
np.savez_compressed(sys.argv[2], x=x, y=y, xscol=xl, yscol=yl)
'''
        run([args.numerical_python, '-c', generator, REPO, fixture, spec.nx], temp)
        with np.load(fixture) as source:
            arrays = {key: source[key].copy() for key in source.files}
        rng = np.random.RandomState(39)
        rows = []
        for i in range(len(arrays['x'])):
            x, y = int(arrays['xscol'][i]), int(arrays['yscol'][i])
            rows.append({'cal.etottracknorm': .9, spec.cer+'.npeSum': 8.,
                         'gtr.dp': float(rng.uniform(-4.9, -.1)),
                         'dc.xp_fp': arrays['x'][i, 0], 'dc.x_fp': arrays['x'][i, 1],
                         'dc.yp_fp': arrays['y'][i, 0], 'dc.y_fp': arrays['y'][i, 1],
                         'extcor.xsieve': spec.xs(x)+float(rng.normal(0., .025)),
                         'extcor.ysieve': spec.ys(y)+float(rng.normal(0., .025))})
        replay = temp/(arm+'_campaign_replay.root')
        write_replay(replay, arm, rows, 'T' if arm == 'HMS' else 'Tout')
        source = temp/(arm+'_source'); config = source/'config'; config.mkdir(parents=True)
        table = config/'rungroups_test_inputs.tsv'
        table.write_text('rungroup\toptics_id\tangle_deg\tfoils\tnruns\truns\trootfile\n'
                         'test\t99001\t12.5\t0\t1\t99001\t'+str(replay)+'\n')
        seed = config/'oldfit.dat'
        seed.write_text('! software fixture constant term\n ---------------------------------------------\n'
                        ' 0 0 0 0 00000\n ---------------------------------------------\n')
        campaign = temp/(arm+'_geometric')
        run([args.numerical_python, REPO/'setup_geometric_campaign.py', source, campaign], temp)
        foil_cuts(campaign/'01_ytar_cuts/cuts/ytar_ridge_cut_test.root', [(-1., 1.)])
        command = [args.numerical_python, REPO/'run_geometric_selection.py', campaign,
                   '--rungroup', 'test', '--axis', 'both', '--metadata', metadata,
                   '--root-python', sys.executable]
        run(command, temp)
        if arm == 'HMS':
            check_manual_preview(args, temp, campaign, metadata)
        # Reusing a campaign must require an explicit rerun, never overwrite by default.
        duplicate = subprocess.run([str(x) for x in command], cwd=temp, text=True, capture_output=True)
        assert duplicate.returncode != 0, 'Selector overwrote existing output'
        run([args.numerical_python, REPO/'run_relabel_x_geometric_all.py', campaign,
             '--metadata', metadata, '--root-command', 'root'], temp)
        xcut = campaign/'03b_relabel_x/cuts/XpFpXFp_test_auto_band_cut_xscolRelabeled.root'
        digest = hashlib.sha256(xcut.read_bytes()).hexdigest()
        duplicate = subprocess.run([args.numerical_python, str(REPO/'run_relabel_x_geometric_all.py'),
                                    str(campaign), '--metadata', str(metadata), '--root-command', 'root'],
                                   cwd=temp, text=True, capture_output=True)
        assert duplicate.returncode != 0 and hashlib.sha256(xcut.read_bytes()).hexdigest() == digest
        run([args.numerical_python, REPO/'run_relabel_y_geometric_all.py', campaign,
             '--metadata', metadata, '--root-command', 'root'], temp)
        ycut = campaign/'03a_relabel_y/cuts/YpFpYFp_test_auto_band_cut_yscolRelabeled.root'
        assert ycut.is_file()
        digest = hashlib.sha256(ycut.read_bytes()).hexdigest()
        duplicate = subprocess.run([args.numerical_python, str(REPO/'run_relabel_y_geometric_all.py'),
                                    str(campaign), '--metadata', str(metadata), '--root-command', 'root'],
                                   cwd=temp, text=True, capture_output=True)
        assert duplicate.returncode != 0 and hashlib.sha256(ycut.read_bytes()).hexdigest() == digest
        candidates = {}
        for axis in ('y', 'x'):
            macro('make_'+axis+'scol_candidate_tree.C',
                  [99001, 'test', str(campaign), str(replay)], temp)
            candidates[axis] = inspect_candidates(campaign, axis, arrays[axis+'scol'])
            assert len(candidates[axis]) >= .995*len(rows), (axis, len(candidates[axis]), len(rows))
            if args.gmm_python:
                run([args.gmm_python, REPO/('gmm_cleanup_'+axis+'scol_candidates_v2.py'),
                     '--campaign', campaign, '--rungroup', 'test', '--foil', 0, '--ndel', 0,
                     '--n-init', 1, '--keep-frac', 1., '--min-events', 30], temp)
            else:
                synthetic_masks(campaign, axis, candidates[axis])
        mask_y = campaign/'05a_gmm_cleanup_y/root'; mask_x = campaign/'05b_gmm_cleanup_x/root'
        macro('make_fit_ntuple_from_gmm.C',
              [99001, -1, 'test', '', str(mask_y), str(mask_x), 6., .65, False,
               '', str(campaign/'06a_fit_ntuple'), str(replay), 'test'], temp)
        output = campaign/'06a_fit_ntuple/root/Optics_99001_-1_fit_tree_gmm.root'
        f = ROOT.TFile(str(output)); tree = f.Get('TFit'); assert tree
        entries = int(tree.GetEntries()); assert entries >= .99*len(rows)
        assert str(f.Get('hallc_spectrometer').GetTitle()) == arm
        fit_counts = {}
        for event in tree:
            i = int(event.entry)
            assert int(event.xscol) == int(arrays['xscol'][i])
            assert int(event.yscol) == int(arrays['yscol'][i])
            label = int(event.yscol); fit_counts[label] = fit_counts.get(label, 0)+1
        f.Close()
        result = macro('fit_opt_matrix_gmm.C',
                       ['toy', -1, 100000, -1, str(campaign/'06a_fit_ntuple/root'),
                        str(campaign/'06b_svd_fit'), str(campaign/'config/oldfit.dat'),
                        str(campaign/'config/rungroups_test_inputs.tsv'), str(metadata)],
                       temp, allow_diagnostic_error=True)
        # The stable SVD intentionally caps each (foil, delta, yscol) at 1000.
        svd_entries = sum(min(count, 1000) for count in fit_counts.values())
        assert 'number to fit = '+str(svd_entries)+' ' in result.stdout, result.stdout+result.stderr
        assert (campaign/('06b_svd_fit/matrices/nps_'+arm.lower()+'_newfit_toy.dat')).is_file()
        checks = {}
        mask_counts = {}
        for axis in ('x', 'y'):
            cut = campaign/('02b_angle_scan_x' if axis == 'x' else '02a_angle_scan_y')/'cuts'/(
                ('XpFpXFp' if axis == 'x' else 'YpFpYFp')+'_test_auto_band_cut.root')
            root = ROOT.TFile(str(cut))
            checks[axis] = json.loads(str(root.Get('geometric_cut_membership_checks').GetTitle()))
            root.Close()
            directory = mask_x if axis == 'x' else mask_y
            count = 0
            for path in directory.glob('*.root'):
                root = ROOT.TFile(str(path)); mask = root.Get('GMMClean')
                count += sum(int(event.gmm_keep) == 1 for event in mask); root.Close()
            mask_counts[axis] = count
        summaries.append(dict(spectrometer=arm, replay_events=len(rows),
                              x_candidates=len(candidates['x']), y_candidates=len(candidates['y']),
                              gmm_kept=mask_counts, fit_ntuple_events=entries,
                              svd_events=svd_entries, root_membership_checks=checks,
                              used_legacy_gmm=bool(args.gmm_python)))
        print('PASS:', arm, 'selector/cuts -> relabel -> candidates -> GMM -> TFit -> toy SVD;',
              entries, 'TFit events;', svd_entries, 'SVD events')
    (temp/'software_smoke_summary.json').write_text(json.dumps(summaries, indent=2))
    print('Legacy GMM:', 'executed' if args.gmm_python else
          'not executed; compatible synthetic masks used (requires combined PyROOT/sklearn environment)')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--numerical-python', default=sys.executable)
    parser.add_argument('--gmm-python', help='Python importing ROOT, sklearn, NumPy, and Matplotlib')
    parser.add_argument('--keep-output', type=Path,
                        help='Fresh directory for inspecting fixtures; default is temporary')
    args = parser.parse_args()
    context = tempfile.TemporaryDirectory(prefix='hallc-geometric-root-') if args.keep_output is None else None
    temp = Path(context.name) if context else args.keep_output.absolute()
    if not context: temp.mkdir(parents=True, exist_ok=False)
    try:
        check_export(args.numerical_python, sys.executable, temp)
        check_ridge(temp)
        check_campaign(args, temp)
    finally:
        if context: context.cleanup()


if __name__ == '__main__':
    main()
