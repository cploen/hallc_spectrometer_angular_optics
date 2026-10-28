"""Campaign adapter: replay -> feature-only surfaces -> named candidates.

The fit stage receives measured focal-plane coordinates and delta only. The
separate label stage can use reconstructed sieve positions solely to choose
one orientation and one integer origin for each already frozen band family.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from run_relabel_x_geometric_all import load_groups
from spectrometer_config import from_campaign, run_metadata

STAGES = ('ridge', 'export', 'fit', 'label', 'candidates')


def plain(value):
    import numpy as np
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    return value


def save(path, value):
    path = Path(path)
    if os.path.lexists(path):
        raise FileExistsError('Refusing to overwrite ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(plain(value), stream, indent=2, allow_nan=False)
        stream.write('\n')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def work(campaign, tag):
    return campaign / '00_joint_surface_work' / tag


def load_config(campaign):
    spec = from_campaign(campaign)
    marker = campaign / 'config/joint_surface_campaign.json'
    if not marker.is_file():
        raise ValueError('Use setup_joint_surface_campaign.py to create a fresh joint-surface campaign')
    value = json.loads(marker.read_text())
    if (value.get('selection_method') != 'joint_surface_provisional' or
            value.get('spectrometer') != spec.name or value.get('experimental') is not True):
        raise ValueError('Invalid joint-surface campaign marker: ' + str(marker))
    legacy_veto = campaign / 'config/pre_gmm_veto.tsv'
    if legacy_veto.is_file() and any(line.strip() and not line.lstrip().startswith('#')
                                     for line in legacy_veto.read_text().splitlines()):
        raise ValueError('Active config/pre_gmm_veto.tsv contains legacy TCutG vetoes, which do not apply '
                         'to continuous surface assignments. Move it to previous_selection_reference; '
                         'no legacy cut veto has been applied. Review global_maps or geometry parameters instead.')
    path = campaign / 'config/joint_surface.json'
    config = json.loads(path.read_text())
    allowed = {'schema_version', 'selection', 'fit', 'global_maps', 'coverage_min_events'}
    if set(config) - allowed:
        raise ValueError('Unknown joint-surface configuration keys: ' + str(sorted(set(config) - allowed)))
    if config.get('schema_version') != 1:
        raise ValueError('Expected joint-surface configuration schema_version 1')
    selection = config.get('selection', {})
    if set(selection) - {'cal_cut', 'cer_cut', 'tree', 'hash_replay'}:
        raise ValueError('Unknown selection option')
    if config.get('coverage_min_events', 10) < 1:
        raise ValueError('coverage_min_events must be positive')
    if set(config.get('fit', {})) - {'common', 'x', 'y'}:
        raise ValueError('fit configuration contains only common, x and y objects')
    return spec, config


def run_io(mode, payload, base, root_python):
    file = base / (mode + '_payload.json')
    save(file, payload)
    command = [root_python, str(REPO / 'joint_surface/root_io.py'), mode, '--payload', str(file)]
    if mode == 'candidates':
        command += ['--report', str(base / 'candidate_report.json')]
    print(shlex.join(command), flush=True)
    result = subprocess.run(command, cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log = base / (mode + '.log')
    with log.open('x') as stream:
        stream.write(result.stdout)
    if result.returncode:
        raise RuntimeError('ROOT I/O failed; inspect ' + str(log) + '\n' + result.stdout[-4000:])
    print(result.stdout[-1500:].rstrip(), flush=True)


def verify_geometry(base):
    record = json.loads((base / 'geometry_frozen.json').read_text())
    for relative, digest in record['sha256'].items():
        if sha256(base / relative) != digest:
            raise ValueError('Frozen geometry changed: ' + relative + '; use a fresh campaign for refitting')
    return record


def fit_group(campaign, row, meta, config, spec):
    import numpy as np
    import scipy
    from joint_surface.engine import fit_axis, predict_axis
    base = work(campaign, row['rungroup'])
    if (base / 'fit_started.json').exists():
        raise FileExistsError('Fit already started; inspect this attempt or use a fresh campaign: ' + str(base))
    with np.load(base / 'preband/preband_events.npz', allow_pickle=False) as source:
        data = {k: source[k].copy() for k in source.files}
    save(base / 'fit_started.json', {'config': config, 'feature_inputs': ['xfp', 'xpfp', 'yfp', 'ypfp', 'delta'],
                                    'historical_labels_used': False, 'reconstructed_sieve_used': False,
                                    'python_version': sys.version, 'numpy_version': np.__version__,
                                    'scipy_version': scipy.__version__,
                                    'source_sha256': {name: sha256(REPO / name) for name in
                                                      ('joint_surface/engine.py', 'joint_surface/campaign.py',
                                                       'spectrometer_config.py')}})
    frozen = {'fit_started.json': sha256(base / 'fit_started.json'),
              'preband/preband_events.npz': sha256(base / 'preband/preband_events.npz')}
    if (base / 'preband/provenance.json').is_file():
        frozen['preband/provenance.json'] = sha256(base / 'preband/provenance.json')
    for foil, zfoil in enumerate(meta['foils']):
        take = data['foil'] == foil
        folder = base / ('foil' + str(foil))
        for axis in ('x', 'y'):
            selected_config = dict(config.get('fit', {}).get('common', {}),
                                   **config.get('fit', {}).get(axis, {}))
            selected_config.setdefault('max_slots', spec.nx if axis == 'x' else spec.ny)
            selected_config.setdefault('coverage_min_events', config.get('coverage_min_events', 10))
            raw = np.column_stack((data[axis + 'pfp'][take], data[axis + 'fp'][take]))
            delta = data['delta'][take]
            print(f"Fitting {row['rungroup']} foil{foil} {axis.upper()}: {len(delta)} pre-band events", flush=True)
            finite = np.isfinite(np.column_stack((raw, delta))).all(axis=1)
            finite_indices = np.flatnonzero(finite)
            model, metadata = fit_axis(raw[finite], delta[finite], meta['edges'], config=selected_config)
            # Expose the actual training membership in the complete foil table.
            for name in ('training_indices', 'model_choice_indices'):
                metadata[name] = finite_indices[np.asarray(metadata[name], int)].tolist()
            metadata['preband_events'] = len(delta)
            metadata['nonfinite_feature_events'] = int((~finite).sum())
            model_file = folder / (axis + '_model.json')
            save(model_file, model)
            save(folder / (axis + '_fit_metadata.json'), metadata)
            prediction = predict_axis(model, raw, delta, coverage_min_events=config.get('coverage_min_events', 10))
            output = folder / (axis + '_prediction.npz')
            if output.exists():
                raise FileExistsError(str(output))
            np.savez_compressed(output, **prediction)
            for artifact in (model_file, folder / (axis + '_fit_metadata.json'), output):
                frozen[str(artifact.relative_to(base))] = sha256(artifact)
    save(base / 'geometry_frozen.json', {'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(),
                                       'sha256': frozen, 'historical_labels_used': False,
                                       'reconstructed_sieve_used_in_fit': False})


def choose_global_map(latent, reconstructed, slots, pitch, coordinate, manual=None):
    """One whole-family affine integer map; never event-by-event relabeling."""
    import numpy as np
    latent = np.asarray(latent, int)
    reconstructed = np.asarray(reconstructed, float)
    valid = (latent >= 0) & np.isfinite(reconstructed)
    ids = np.unique(latent[latent >= 0])
    if not len(ids):
        raise ValueError('No accepted latent bands available for physical naming')
    if manual is not None:
        if set(manual) != {'orientation', 'origin'}:
            raise ValueError('Manual global map needs only orientation and origin')
        sign, origin = manual['orientation'], manual['origin']
        if sign not in (-1, 1) or not isinstance(origin, int):
            raise ValueError('Map orientation must be -1 or 1; origin must be an integer')
        if np.any(sign * ids + origin < 0) or np.any(sign * ids + origin >= slots):
            raise ValueError('Manual global map puts a discovered band outside the physical grid')
        return {'orientation': sign, 'origin': origin, 'source': 'manual_global_map',
                'reconstructed_sieve_used': False, 'latent_ids': ids.tolist()}
    medians = {int(index): float(np.median(reconstructed[valid & (latent == index)]))
               for index in ids if np.any(valid & (latent == index))}
    if not medians:
        raise ValueError('No finite reconstructed sieve positions for naming; provide a manual global map')
    options = []
    for sign in (-1, 1):
        for origin in range(-slots, 2 * slots):
            mapped = sign * ids + origin
            if np.any(mapped < 0) or np.any(mapped >= slots):
                continue
            errors = [abs(median - coordinate(sign * index + origin)) / pitch for index, median in medians.items()]
            options.append({'orientation': sign, 'origin': origin, 'median_grid_error_pitch': float(np.median(errors)),
                            'maximum_band_median_error_pitch': float(np.max(errors))})
    if not options:
        raise ValueError('Discovered latent index span exceeds the physical grid; review fit, do not patch individual labels')
    options.sort(key=lambda value: (value['median_grid_error_pitch'], value['maximum_band_median_error_pitch']))
    return dict(options[0], source='global_reconstructed_sieve_alignment', reconstructed_sieve_used=True,
                latent_ids=ids.tolist(), reconstructed_band_medians_cm=medians,
                alternative_maps=options[1:4])


def apply_map(labels, mapping):
    import numpy as np
    labels = np.asarray(labels, int)
    return np.where(labels >= 0, mapping['orientation'] * labels + mapping['origin'], -1).astype(np.int32)


def label_group(campaign, row, meta, config, spec):
    import numpy as np
    base = work(campaign, row['rungroup'])
    verify_geometry(base)
    fitted_config = json.loads((base / 'fit_started.json').read_text())['config']
    if ({k: v for k, v in fitted_config.items() if k != 'global_maps'} !=
            {k: v for k, v in config.items() if k != 'global_maps'}):
        raise ValueError('Fit/selection/count configuration changed after fitting; use a fresh campaign. '
                         'Only global_maps may be edited between fit and label stages.')
    if (base / 'physical_maps.json').exists():
        raise FileExistsError('Physical maps already exist; use a fresh campaign for different maps')
    with np.load(base / 'preband/preband_events.npz', allow_pickle=False) as source:
        data = {k: source[k].copy() for k in source.files}
    maps, cases = {}, []
    for foil, zfoil in enumerate(meta['foils']):
        take = data['foil'] == foil
        folder = base / ('foil' + str(foil))
        member = {key: data[key][take] for key in ('entry', 'xfp', 'xpfp', 'yfp', 'ypfp', 'delta')}
        key = row['rungroup'] + '/foil' + str(foil)
        maps[key] = {}
        for axis in ('x', 'y'):
            with np.load(folder / (axis + '_prediction.npz'), allow_pickle=False) as source:
                prediction = {k: source[k].copy() for k in source.files}
            slots, pitch = (spec.nx, spec.dx) if axis == 'x' else (spec.ny, spec.dy)
            coordinate = spec.xs if axis == 'x' else spec.ys
            mapping = choose_global_map(prediction['labels'], data[axis + 'sieve'][take], slots, pitch, coordinate,
                                        config.get('global_maps', {}).get(key, {}).get(axis))
            maps[key][axis] = mapping
            member[axis + '_label'] = apply_map(prediction['labels'], mapping)
            member['physical_' + axis] = member[axis + '_label']
            member['latent_' + axis] = prediction['labels'].astype(np.int32)
            member['proposed_' + axis] = prediction.get('proposed_label', prediction['labels']).astype(np.int32)
            member['accepted_' + axis] = np.asarray(prediction['accepted'], bool)
            member['reason_' + axis] = np.asarray(prediction['reason'], str)
            fit_metadata = json.loads((folder / (axis + '_fit_metadata.json')).read_text())
            for source_name, mask_name in (('training_indices', 'train_' + axis),
                                           ('model_choice_indices', 'model_choice_' + axis)):
                member[mask_name] = np.zeros(len(member['entry']), bool)
                member[mask_name][np.asarray(fit_metadata[source_name], int)] = True
            for name in ('q', 'r', 'u', 'posterior', 'coverage', 'occupancy', 'interpolation'):
                if name in prediction:
                    member[axis + '_' + name] = prediction[name]
            print(f"{key} {axis.upper()}: physical = {mapping['orientation']} * latent + {mapping['origin']} ({mapping['source']})", flush=True)
        output = folder / 'membership.npz'
        if output.exists():
            raise FileExistsError(str(output))
        np.savez_compressed(output, **member)
        relative = str(output.relative_to(campaign / '00_joint_surface_work'))
        cases.append({'key': row['rungroup'] + '_foil' + str(foil),
                      'label': f"{row['rungroup']} foil {zfoil:g} cm", 'runs': row.get('runs', ''),
                      'foil_z_cm': zfoil, 'delta_edges': meta['edges'], 'events': relative,
                      'predictions': {'main': relative}, 'label_kind': 'physical'})
    save(base / 'physical_maps.json', {'geometry_sha256': sha256(base / 'geometry_frozen.json'), 'maps': maps})
    artifacts = [base / 'physical_maps.json'] + [base / ('foil' + str(foil)) / 'membership.npz'
                                               for foil in range(len(meta['foils']))]
    save(base / 'labels_frozen.json', {'sha256': {str(path.relative_to(base)): sha256(path) for path in artifacts}})
    save(base / 'viewer_cases.json', cases)
    # Publish/replace the small aggregate index; never change event artifacts.
    aggregate = []
    for item in sorted((campaign / '00_joint_surface_work').glob('*/viewer_cases.json')):
        aggregate.extend(json.loads(item.read_text()))
    manifest = campaign / '00_joint_surface_work/viewer_manifest.json'
    temporary = manifest.with_name('.viewer_manifest.new')
    with temporary.open('w') as stream:
        json.dump({'schema_version': 1, 'cases': aggregate}, stream, indent=2)
    os.replace(temporary, manifest)


def candidate_payload(campaign, row, meta):
    base = work(campaign, row['rungroup'])
    verify_geometry(base)
    frozen = json.loads((base / 'physical_maps.json').read_text())
    if frozen['geometry_sha256'] != sha256(base / 'geometry_frozen.json'):
        raise ValueError('Geometry changed since physical naming')
    labels_frozen = json.loads((base / 'labels_frozen.json').read_text())
    for relative, digest in labels_frozen['sha256'].items():
        if sha256(base / relative) != digest:
            raise ValueError('Frozen membership changed: ' + relative)
    outputs = {}
    for axis, stage, prefix in (('x', '04b_candidate_trees_x', 'Xscol'), ('y', '04a_candidate_trees_y', 'Yscol')):
        stem = prefix + 'Candidates_' + row['rungroup']
        outputs[axis] = {'root': str(campaign / stage / 'root' / (stem + '.root')),
                         'tsv': str(campaign / stage / 'tsv' / (stem + '_summary.tsv'))}
    return {'campaign': str(campaign), 'rungroup': row['rungroup'], 'optics_id': int(row['optics_id']),
            'events': str(base / 'preband/preband_events.npz'), 'foil_z_cm': meta['foils'], 'delta_edges': meta['edges'],
            'memberships': [str(base / ('foil' + str(foil)) / 'membership.npz') for foil in range(len(meta['foils']))],
            'outputs': outputs, 'geometry_frozen_sha256': sha256(base / 'geometry_frozen.json'),
            'physical_maps_sha256': sha256(base / 'physical_maps.json')}


def parser():
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument('campaign', type=Path)
    value.add_argument('stage', nargs='?', choices=STAGES)
    value.add_argument('--plan', action='store_true', help='Print stages through pre-GMM candidates; no farm access needed')
    value.add_argument('--dry-run', action='store_true')
    value.add_argument('--rungroup', action='append', help='Run only this group; repeatable')
    value.add_argument('--root-python', default=os.environ.get('HALLC_ROOT_PYTHON', sys.executable))
    value.add_argument('--root-command', default='hcana', help='ROOT/HCANA executable used only for foil ridge construction')
    value.add_argument('--metadata', type=Path,
                       help='Defaults to campaign config/optics.dat, then repository metadata')
    return value


def main(argv=None):
    args = parser().parse_args(argv)
    if not args.plan and args.stage is None:
        parser().error('Choose one explicit stage or --plan')
    if args.plan and args.stage:
        parser().error('--plan cannot be combined with an execution stage')
    try:
        campaign = args.campaign.expanduser().resolve()
        spec, config = load_config(campaign)
        groups = load_groups(campaign, args.rungroup)
        snapshot = campaign / 'config/optics.dat'
        metadata = (args.metadata if args.metadata is not None else
                    snapshot if snapshot.exists() else REPO / 'DATfiles/list_of_optics_run.dat')
        metadata = metadata.expanduser().resolve()
        if snapshot.exists() and sha256(metadata) != sha256(snapshot):
            raise ValueError('--metadata differs from the campaign metadata snapshot config/optics.dat; '
                             'use a fresh campaign for different metadata')
        stages = STAGES if args.plan else (args.stage,)
        planning = args.plan or args.dry_run
        print(f'Campaign: {campaign} ({spec.name}); continuous-delta joint surfaces', flush=True)
        for row in groups:
            meta = run_metadata(int(row['optics_id']), metadata)
            if spec.name == 'SHMS' and meta['sieve_flag'] != 1:
                raise ValueError('Only centered SHMS sieve (SieveFlag=1) is currently supported')
            base = work(campaign, row['rungroup'])
            replay = Path(row['rootfile']).expanduser()
            replay = replay if replay.is_absolute() else REPO / replay
            foil_cuts = campaign / '01_ytar_cuts/cuts' / ('ytar_ridge_cut_' + row['rungroup'] + '.root')
            for stage in stages:
                print(f"[{stage}] {row['rungroup']}: {len(meta['foils'])} foil(s), delta {meta['edges'][0]:g} to {meta['edges'][-1]:g}%", flush=True)
                if planning:
                    command = [sys.executable, str(REPO / 'run_joint_surface_campaign.py'), str(campaign), stage,
                               '--rungroup', row['rungroup'], '--root-python', args.root_python,
                               '--metadata', str(metadata), '--root-command', args.root_command]
                    print('  ' + shlex.join(command), flush=True)
                    continue
                if stage in ('ridge', 'export') and not replay.is_file():
                    raise FileNotFoundError('Replay not mounted: ' + str(replay))
                if stage == 'ridge':
                    if foil_cuts.exists():
                        raise FileExistsError('Refusing to replace foil cuts: ' + str(foil_cuts))
                    expression = 'ytar_ridge_cut.C(' + ','.join([str(int(row['optics_id']))] +
                        [json.dumps(v) for v in (row['rungroup'], '-1', '-1', str(replay), str(campaign))]) + ')'
                    subprocess.run(shlex.split(args.root_command) + ['-b', '-l', '-q', expression], cwd=REPO, check=True)
                    if not foil_cuts.is_file():
                        raise RuntimeError('Foil ridge macro did not publish expected output: ' + str(foil_cuts))
                elif stage == 'export':
                    if not foil_cuts.is_file():
                        raise FileNotFoundError('Missing foil cuts: ' + str(foil_cuts))
                    selection = dict(config.get('selection', {}))
                    payload = dict(selection, campaign=str(campaign), input=str(replay.resolve()), foil_cuts=str(foil_cuts),
                                   metadata=str(metadata), optics_id=int(row['optics_id']), rungroup=row['rungroup'],
                                   output=str(base / 'preband'))
                    run_io('export', payload, base, args.root_python)
                elif stage == 'fit':
                    fit_group(campaign, row, meta, config, spec)
                elif stage == 'label':
                    label_group(campaign, row, meta, config, spec)
                elif stage == 'candidates':
                    run_io('candidates', candidate_payload(campaign, row, meta), base, args.root_python)
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 1
