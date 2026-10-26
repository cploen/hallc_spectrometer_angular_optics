#!/usr/bin/env python3
"""Alternative step 2: label-free X/Y discovery with conventional ROOT cuts.

Run from the repository root in a fresh campaign made by
setup_geometric_campaign.py. Stable angle-scan scripts remain available.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import uuid
from importlib.metadata import version
import numpy as np
from spectrometer_config import run_metadata
from run_relabel_x_geometric_all import load_groups, require_provisional_campaign
from geometric_selection import fit_x, detect_x, predict_x, fit_y, predict_y, global_coordinates
from geometric_selection.cuts import build_x_polygons, build_y_polygons, polygon_predictions
from geometric_selection.diagnostics import render

REPO = Path(__file__).resolve().parent


def plain(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [plain(v) for v in value]
    return value


def save(path, value):
    with Path(path).open('x') as stream:
        json.dump(plain(value), stream, indent=2, allow_nan=False); stream.write('\n')


def delta_tag(low, high):
    return ('%g_to_%g' % (low, high)).replace('-', 'm').replace('.', 'p')


def cut_file(campaign, axis, tag):
    stage = '02b_angle_scan_x' if axis == 'x' else '02a_angle_scan_y'
    prefix = 'XpFpXFp' if axis == 'x' else 'YpFpYFp'
    return campaign/stage/'cuts'/(prefix+'_'+tag+'_auto_band_cut.root')


def run(command, log=None):
    print(shlex.join([str(x) for x in command]), flush=True)
    result = subprocess.run([str(x) for x in command], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if log:
        with Path(log).open('x') as stream: stream.write(result.stdout)
    if result.returncode:
        raise RuntimeError('Command failed (%d); %s\n%s' % (result.returncode, log or '', result.stdout[-6000:]))


def read_overrides(path, groups, campaign):
    explicit = path is not None
    if path is None:
        path = campaign/'config/geometric_overrides.json'
    if not path.exists():
        if explicit: raise FileNotFoundError('Missing override file '+str(path))
        return {}, None
    value = json.loads(path.read_text())
    if set(value) != {'slices'} or not isinstance(value['slices'], dict):
        raise ValueError('Override JSON must contain only a slices object')
    known = {row['rungroup'] for row in groups}
    for key, pair in value['slices'].items():
        parts = key.split('/')
        if len(parts) != 3 or parts[0] not in known or not re.fullmatch(r'foil\d+', parts[1]) or not re.fullmatch(r'delta\d+', parts[2]):
            raise ValueError('Unknown or malformed override key '+key)
        for axis, item in pair.items():
            allowed = {'angle', 'frozen_model', 'polygons_file'} | ({'learning_delta', 'accepted_indices', 'boundaries_u'} if axis == 'x' else set())
            if axis not in ('x', 'y') or not isinstance(item, dict) or set(item)-allowed:
                raise ValueError('Unknown override fields for '+key+'/'+axis)
            if 'frozen_model' in item and ('angle' in item or 'learning_delta' in item):
                raise ValueError('A frozen model already fixes angle/window; remove the conflicting override for '+key+'/'+axis)
    return value['slices'], path.resolve()


def fit_slice(axis, data, selection, foil_selection, ndel, edges, slots, override, base):
    raw = np.c_[data[axis+'pfp'][selection], data[axis+'fp'][selection]]
    if not len(raw): raise ValueError('No finite target events in this foil/delta slice')
    if len(raw) < 100 and not ('frozen_model' in override or 'polygons_file' in override):
        raise ValueError('Fewer than 100 target events; provide a reviewed model/cut strategy or retain stable/manual selection')
    def resolve(name):
        path = Path(name).expanduser()
        return path if path.is_absolute() else base/path
    if 'polygons_file' in override and 'frozen_model' not in override:
        # Fully manual fallback remains usable even when discovery cannot start.
        polygons = json.loads(resolve(override['polygons_file']).read_text())
        if not isinstance(polygons, list): raise ValueError('polygons_file must contain a list of index/points objects')
        mean = raw.mean(0); scale = raw.std(0, ddof=1) if len(raw) > 1 else np.ones(2)
        scale = np.where(np.isfinite(scale) & (scale > 0), scale, 1.)
        angle = float(override.get('angle', 135. if axis == 'x' else 0.))
        if not np.isfinite(angle): raise ValueError('Manual display angle must be finite')
        theta = np.deg2rad(angle)
        z = (raw-mean)/scale
        labels, matches = polygon_predictions(polygons, raw)
        if (matches > 1).any(): raise ValueError('Manual polygons assign events to several bands')
        transform = dict(mean=mean.tolist(), scale=scale.tolist(), theta_deg=angle)
        model = dict(axis=axis, manual_partition=True, manual_polygons=str(resolve(override['polygons_file'])),
                     overrides=override, target_events=len(raw), **({'geometry': transform} if axis == 'x' else {'model': transform}))
        prediction = dict(labels=labels, q=z @ [np.cos(theta), np.sin(theta)], p=z @ [-np.sin(theta), np.cos(theta)])
        return raw, model, prediction, polygons, None, np.full(len(raw), -1, int)
    frozen = json.loads(resolve(override['frozen_model']).read_text()) if 'frozen_model' in override else None
    if frozen and frozen.get('axis') != axis:
        raise ValueError('Frozen model axis does not match requested slice')
    learning_u = None
    if axis == 'x':
        interval = override.get('learning_delta', [edges[0], edges[2]] if ndel == 0 and len(edges) >= 3 else edges[ndel:ndel+2])
        if len(interval) != 2 or not np.isfinite(interval).all() or interval[0] >= interval[1]:
            raise ValueError('learning_delta needs two increasing finite limits')
        learn = foil_selection & (data['delta'] >= interval[0]) & (data['delta'] < interval[1])
        learning = np.c_[data['xpfp'][learn], data['xfp'][learn]]
        geometry = frozen['geometry'] if frozen else fit_x(learning, ndel=ndel, angle=override.get('angle'), max_slots=slots)
        if len(geometry['possible_indices']) != slots:
            raise ValueError('Frozen X model has incompatible spectrometer slot count')
        occupancy = frozen['occupancy'] if frozen else detect_x(geometry, raw)
        occupancy = {k: occupancy[k] for k in ('accepted_indices', 'mixture_model', 'candidates')}
        if 'accepted_indices' in override:
            ids = override['accepted_indices']
            if not ids or any(not isinstance(i, int) or i < 0 or i >= slots for i in ids) or ids != sorted(set(ids)):
                raise ValueError('accepted_indices must be increasing distinct physical slot hypotheses')
            occupancy['accepted_indices'] = ids
        prediction = predict_x(geometry, occupancy, raw)
        model = dict(axis=axis, geometry=geometry, occupancy=occupancy, learning_delta=interval,
                     learning_events=int(learn.sum()), target_events=len(raw))
        if 'boundaries_u' in override: model['boundaries_u'] = override['boundaries_u']
        polygons = build_x_polygons(geometry, occupancy, raw, boundaries_u=override.get('boundaries_u'))
        if 'boundaries_u' in override:
            valid = prediction['finite_geometry']
            prediction['labels'][valid] = np.asarray(occupancy['accepted_indices'])[np.searchsorted(override['boundaries_u'], prediction['u'][valid])]
        learning_u = global_coordinates(geometry, learning)['u']
    else:
        fitted = frozen['model'] if frozen else fit_y(raw, angle=override.get('angle'), max_slots=slots)
        if fitted.get('config', {}).get('max_slots', slots) != slots or fitted['band_count'] > slots:
            raise ValueError('Frozen Y model has incompatible spectrometer slot count')
        details = predict_y(fitted, raw, details=True)
        prediction = dict(labels=details['labels'], q=details['coordinates'][:, 0], p=details['coordinates'][:, 1],
                          posterior=details['posterior'], mahalanobis2=details['mahalanobis2'])
        model = dict(axis=axis, model=fitted, target_events=len(raw))
        polygons = build_y_polygons(fitted, raw)
    automatic_labels = prediction['labels'].copy()
    if 'polygons_file' in override:
        polygons = json.loads(resolve(override['polygons_file']).read_text())
        if not isinstance(polygons, list): raise ValueError('polygons_file must contain a list of index/points objects')
        prediction['labels'], matches = polygon_predictions(polygons, raw)
        if (matches > 1).any(): raise ValueError('Manual polygons assign events to several bands')
        model['manual_polygons'] = str(resolve(override['polygons_file']))
    elif axis == 'y':
        # Contours approximate curved decision regions. Refine only when actual
        # event membership requires it; no retuning of the fitted statistics.
        for grid in (401, 801, 1601):
            if grid != 401: polygons = build_y_polygons(fitted, raw, grid=grid)
            exported, matches = polygon_predictions(polygons, raw)
            if not (matches > 1).any() and np.mean(exported != prediction['labels']) <= .002: break
        model['polygon_grid'] = grid
    exported, matches = polygon_predictions(polygons, raw)
    model['polygon_check'] = dict(multiple_matches=int((matches > 1).sum()),
                                  disagreement_fraction=float(np.mean(exported != prediction['labels'])))
    model['overrides'] = override
    if model['polygon_check']['multiple_matches'] or model['polygon_check']['disagreement_fraction'] > .002:
        raise ValueError('Polygon approximation needs review: '+json.dumps(model['polygon_check']))
    return raw, model, prediction, polygons, learning_u, automatic_labels


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--axis', choices=['x', 'y', 'both'], default='both')
    parser.add_argument('--rungroup', action='append')
    parser.add_argument('--metadata', type=Path, default=REPO/'DATfiles/list_of_optics_run.dat')
    parser.add_argument('--root-python', default=os.environ.get('HALLC_ROOT_PYTHON', sys.executable))
    parser.add_argument('--overrides', type=Path)
    parser.add_argument('--cal-cut', type=float, default=.65)
    parser.add_argument('--cer-cut', type=float, default=6.)
    parser.add_argument('--preview', action='store_true', help='Write editable preview cuts and diagnostics outside canonical cuts')
    parser.add_argument('--foil', type=int, help='Only with --preview')
    parser.add_argument('--ndel', type=int, help='Only with --preview')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args(argv)
    try:
        campaign = args.campaign.resolve()
        spec = require_provisional_campaign(campaign)
        groups = load_groups(campaign, args.rungroup)
        all_groups = load_groups(campaign)
        overrides, override_path = read_overrides(args.overrides, all_groups, campaign)
        valid_keys = set()
        for row in all_groups:
            meta = run_metadata(int(row['optics_id']), args.metadata)
            valid_keys.update('%s/foil%d/delta%d' % (row['rungroup'], foil, ndel)
                              for foil in range(len(meta['foils'])) for ndel in range(len(meta['edges'])-1))
        if set(overrides)-valid_keys:
            raise ValueError('Override addresses nonexistent foil/delta: '+', '.join(sorted(set(overrides)-valid_keys)))
        if not args.preview and (args.foil is not None or args.ndel is not None):
            raise ValueError('Subset foil/delta selection requires --preview; canonical files must cover the rungroup')
        axes = ['x', 'y'] if args.axis == 'both' else [args.axis]
        for row in groups:
            meta = run_metadata(int(row['optics_id']), args.metadata)
            if spec.name == 'SHMS' and meta['sieve_flag'] != 1:
                raise ValueError('Only centered SHMS sieve is supported')
            if not args.preview:
                for axis in axes:
                    if cut_file(campaign, axis, row['rungroup']).exists():
                        raise FileExistsError('Cuts already exist; use --preview for review or create a fresh campaign: '+str(cut_file(campaign, axis, row['rungroup'])))
            print('%s: %d foils x %d delta slices; %s' % (row['rungroup'], len(meta['foils']), len(meta['edges'])-1, ','.join(axes)))
        if args.dry_run:
            return 0
        attempt = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
        for row in groups:
            tag = row['rungroup']
            meta = run_metadata(int(row['optics_id']), args.metadata)
            job = campaign/'00_geometric_work'/tag/attempt
            job.mkdir(parents=True, exist_ok=False)
            foil_cuts = campaign/'01_ytar_cuts/cuts'/('ytar_ridge_cut_'+tag+'.root')
            command = [args.root_python, REPO/'geometric_selection/export_preband.py', '--input', row['rootfile'],
                       '--foil-cuts', foil_cuts, '--run-list', args.metadata, '--metadata-run', row['optics_id'],
                       '--campaign', campaign, '--rungroup', tag, '--output', job/'preband',
                       '--cal-cut', str(args.cal_cut), '--cer-cut', str(args.cer_cut)]
            run(command, job/'export.log')
            with np.load(job/'preband/preband_events.npz') as source:
                data = {k: source[k].copy() for k in source.files}
            for axis in axes:
                stage = campaign/('02b_angle_scan_x' if axis == 'x' else '02a_angle_scan_y')
                work = job/axis; work.mkdir()
                payload = dict(schema_version=1, axis=axis, spectrometer=spec.name, campaign=str(campaign),
                               rungroup=tag, provisional=True, no_x_region_split=True, maximum_cut_disagreement=.002,
                               export_provenance=str(job/'preband/provenance.json'), metadata=str(args.metadata.resolve()),
                               override_file=str(override_path) if override_path else None,
                               cut_output=str(work/'preview_cuts.root' if args.preview else cut_file(campaign, axis, tag)), slices=[])
                payload['numerical_versions'] = {name: version(name) for name in ('numpy', 'scipy', 'matplotlib')}
                for foil in range(len(meta['foils'])):
                    for ndel, (low, high) in enumerate(zip(meta['edges'][:-1], meta['edges'][1:])):
                        if args.foil is not None and foil != args.foil: continue
                        if args.ndel is not None and ndel != args.ndel: continue
                        key = '%s/foil%d/delta%d' % (tag, foil, ndel)
                        override = overrides.get(key, {}).get(axis, {})
                        finite = np.isfinite(data[axis+'pfp']) & np.isfinite(data[axis+'fp']) & np.isfinite(data['delta'])
                        fs = (data['foil'] == foil) & finite
                        select = fs & (data['ndel'] == ndel)
                        base = axis+'fp_'+axis+'pfp_geometric_'+tag+'_foil%d_delta_%s' % (foil, delta_tag(low, high))
                        try:
                            raw, model, prediction, polygons, learning, original = fit_slice(axis, data, select, fs, ndel,
                                meta['edges'], spec.nx if axis == 'x' else spec.ny, override,
                                override_path.parent if override_path else campaign/'config')
                        except Exception as exc:
                            save(work/('unresolved_foil%d_ndel%d.json' % (foil, ndel)), dict(slice=key, error=str(exc), overrides=override))
                            raise RuntimeError('%s %s unresolved: %s; diagnostics retained in %s' % (key, axis, exc, work)) from exc
                        model_path = work/(base+'_model.json'); save(model_path, model)
                        arrays = {k: v[select] for k, v in data.items()}
                        arrays.update(model_label=prediction['labels'], automatic_label=original,
                                      **{k: v for k, v in prediction.items() if k != 'labels'})
                        assignments = work/(base+'_assignments.npz'); np.savez_compressed(assignments, **arrays)
                        plot_path = stage/'plots'/attempt/(base+'.pdf')
                        render(plot_path, axis, raw, prediction, polygons, model,
                               '%s foil %d (z=%g cm), delta [%g,%g)' % (tag, foil, meta['foils'][foil], low, high), arrays, learning)
                        payload['slices'].append(dict(foil=foil, ndel=ndel, delta_low=low, delta_high=high,
                            model=str(model_path), assignments=str(assignments), polygons=polygons,
                            diagnostic_root=str(stage/'root'/attempt/(base+'.root')), diagnostic_pdf=str(plot_path)))
                        print('%s %s: %d provisional bands, N=%d' % (key, axis, len(np.unique(prediction['labels'][prediction['labels'] >= 0])), len(raw)), flush=True)
                if not payload['slices']:
                    raise ValueError('No requested foil/delta slices exist')
                payload['source_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                            for p in (REPO/'geometric_selection').glob('*.py')}
                payload_path = work/'root_payload.json'; save(payload_path, payload)
                run([args.root_python, REPO/'geometric_selection/root_io.py', '--payload', payload_path,
                     '--report', work/'root_membership_report.json'], work/'root_writer.log')
                checks = json.loads((work/'root_membership_report.json').read_text())
                summary = stage/'tsv'/attempt/(axis+'fp_'+axis+'pfp_geometric_'+tag+'_summary.tsv')
                summary.parent.mkdir(parents=True, exist_ok=True)
                fields = ['tag', 'run', 'foil', 'ndel', 'deltaMin', 'deltaMax', 'thetaDeg',
                          'Nevents', 'nBands', 'learningDeltaMin', 'learningDeltaMax',
                          'cutComponents', 'cutUnassigned', 'cutClassifierDisagreements', 'manualPartition']
                with summary.open('x', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t'); writer.writeheader()
                    for item, check in zip(payload['slices'], checks):
                        stored = json.loads(Path(item['model']).read_text())
                        transform = stored['geometry'] if axis == 'x' else stored['model']
                        window = stored.get('learning_delta', [item['delta_low'], item['delta_high']])
                        with np.load(item['assignments']) as assignments:
                            labels = assignments['model_label']; bands = len(np.unique(labels[labels >= 0]))
                        writer.writerow(dict(tag=tag, run=row['optics_id'], foil=item['foil'], ndel=item['ndel'],
                            deltaMin=item['delta_low'], deltaMax=item['delta_high'], thetaDeg=transform['theta_deg'],
                            Nevents=check['events'], nBands=bands, learningDeltaMin=window[0], learningDeltaMax=window[1],
                            cutComponents=check['cut_components'], cutUnassigned=check['unassigned'],
                            cutClassifierDisagreements=check['classifier_disagreements'], manualPartition=stored.get('manual_partition', False)))
                print('Wrote '+payload['cut_output'], flush=True)
        return 0
    except Exception as exc:
        print('ERROR: '+str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
