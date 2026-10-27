#!/usr/bin/env python3
"""Reproduce the bounded 14,000-event numerical campaign adapter smoke.

This generates controlled pre-band inputs, fits both axes, freezes geometry,
assigns physical names, and prepares the conventional candidate payload. It
needs NumPy/SciPy but no ROOT and never reads historical event labels. The
generated truth is checked only after fitting and naming. ROOT/GMM contracts
are tested separately by joint_surface_root_smoke.py.
"""
import argparse
from contextlib import nullcontext
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'tests'))
from test_joint_surface_engine import sample
from joint_surface.campaign import (candidate_payload, fit_group, label_group,
                                    save, sha256, verify_geometry)
from spectrometer_config import from_campaign


def exercise(directory):
    campaign = directory / 'HMS_numeric'
    base = campaign / '00_joint_surface_work/test'
    (base / 'preband').mkdir(parents=True, exist_ok=False)
    raw, delta, truth, r, frame = sample(14000, count=7)
    entries = 2 * np.arange(len(delta), dtype=np.int64)
    edges = [-10., -8., -5., 0., 5., 10.]
    # Identical controlled patterns exercise both independent adapter paths.
    # Neither the generated truth nor reconstructed sieve positions enter fit.
    data = dict(entry=entries, foil=np.zeros(len(delta), np.int32),
                ndel=(np.searchsorted(edges, delta, side='right') - 1).astype(np.int32),
                xfp=raw[:, 1], xpfp=raw[:, 0], yfp=raw[:, 1], ypfp=raw[:, 0], delta=delta,
                xsieve=(truth - 4) * 2.54, ysieve=(truth - 4) * 1.524,
                sumnpe=np.full(len(delta), 8.), etracknorm=np.full(len(delta), .9),
                xtar=np.zeros(len(delta)), ytar=np.zeros(len(delta)),
                xptar=np.zeros(len(delta)), yptar=np.zeros(len(delta)))
    np.savez_compressed(base / 'preband/preband_events.npz', **data)
    config = {'schema_version': 1, 'coverage_min_events': 10,
              'fit': {'common': {'phase_starts': 4, 'phase_iterations': 60, 'iterations': 3,
                                 'irls': 4, 'max_fit_events': 7000,
                                 'families': ['bilinear', 'curved', 'shift', 'spacing']}},
              'global_maps': {}}
    row = {'rungroup': 'test', 'optics_id': '99002', 'runs': 'synthetic'}
    meta = {'foils': [0.], 'edges': edges}
    spec = from_campaign(campaign)
    save(directory / 'numeric_smoke_configuration.json', {'config': config, 'row': row, 'metadata': meta,
                                                        'generated_events': len(delta), 'generated_bands': 7})
    fit_group(campaign, row, meta, config, spec)
    verify_geometry(base)
    frozen_digest = sha256(base / 'geometry_frozen.json')
    label_group(campaign, row, meta, config, spec)
    payload = candidate_payload(campaign, row, meta)
    save(base / 'numeric_smoke_candidate_payload.json', payload)
    assert sha256(base / 'geometry_frozen.json') == frozen_digest
    verify_geometry(base)
    summary = {'events': len(delta), 'historical_inputs_used': False, 'axes': {}}
    with np.load(base / 'foil0/membership.npz', allow_pickle=False) as member:
        np.testing.assert_array_equal(member['entry'], entries)
        for axis in ('x', 'y'):
            accepted = member['accepted_' + axis]
            assert np.array_equal(accepted, member[axis + '_label'] >= 0)
            np.testing.assert_array_equal(member['entry'][accepted], entries[accepted])
            training, choice = member['train_' + axis], member['model_choice_' + axis]
            assert not np.any(training & choice)
            assert np.all(training | choice)
            assert int(training.sum()) == 11200
            assert int(choice.sum()) == 2800
            acceptance = float(accepted.mean())
            agreement = float(np.mean(member[axis + '_label'][accepted] == truth[accepted]))
            assert acceptance > .95, (axis, acceptance)
            assert agreement > .995, (axis, agreement)
            model = json.loads((base / 'foil0' / (axis + '_model.json')).read_text())
            metadata = json.loads((base / 'foil0' / (axis + '_fit_metadata.json')).read_text())
            assert not metadata['reference_labels_used']
            assert not metadata['reconstructed_sieve_inputs_used']
            assert not metadata['physical_names_assigned']
            assert model['count'] == 7
            summary['axes'][axis] = {'accepted': int(accepted.sum()), 'acceptance': acceptance,
                                      'truth_agreement_among_accepted': agreement,
                                      'geometry_training_events': int(training.sum()),
                                      'model_choice_events': int(choice.sum()), 'family': model['family']}
    maps = json.loads((base / 'physical_maps.json').read_text())['maps']['test/foil0']
    for axis in ('x', 'y'):
        assert maps[axis]['orientation'] == 1
        assert maps[axis]['origin'] == 0
        assert maps[axis]['source'] == 'global_reconstructed_sieve_alignment'
    summary['viewer_manifest'] = str(campaign / '00_joint_surface_work/viewer_manifest.json')
    # Keep entries identical but change a measured coordinate. Publication must
    # refuse it even though every downstream entry join would still succeed.
    preband = base / 'preband/preband_events.npz'
    original_bytes = preband.read_bytes()
    tampered = dict(data)
    tampered['xfp'] = data['xfp'].copy()
    tampered['xfp'][0] += .01
    try:
        np.savez_compressed(preband, **tampered)
        try:
            candidate_payload(campaign, row, meta)
        except ValueError as exc:
            assert 'preband/preband_events.npz' in str(exc)
        else:
            raise AssertionError('Changed pre-band coordinates passed the frozen-data check')
    finally:
        preband.write_bytes(original_bytes)
    verify_geometry(base)
    summary['changed_preband_coordinates_refused'] = True
    save(directory / 'numeric_smoke_summary.json', summary)
    print(json.dumps(summary, indent=2))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        help='Keep artifacts in this new directory, including HMS_numeric and its viewer manifest')
    args = parser.parse_args(argv)
    if args.output is None:
        context = tempfile.TemporaryDirectory(prefix='joint-surface-numeric-')
    else:
        directory = args.output.expanduser().absolute()
        if directory.exists() or directory.is_symlink():
            raise FileExistsError('Refusing to overwrite ' + str(directory))
        directory.mkdir(parents=True)
        context = nullcontext(str(directory))
    with context as directory:
        exercise(Path(directory))
        if args.output is None:
            print('PASS; temporary artifacts removed. Use --output NEW_DIRECTORY to retain a viewer fixture.')
        else:
            print('PASS; artifacts retained in ' + str(args.output.expanduser().absolute()))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
