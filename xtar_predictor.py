"""Separate single-response x-target model. No transport-matrix writes or replay hooks.

The caller must supply an audited measured-input design, frozen training IDs,
folds and scoring cells. This module does not select events or derive features
from foil/hole/truth labels. Native target and exported predictions are cm;
beam scores are pooled cell MSE in mm^2. Production feature choice is pending
the installed HCANA/input audit.
"""
import json
import numpy as np
from elastic_net import DEFAULTS as EN_DEFAULTS, tune, validate
from beam_search import DEFAULTS as BEAM_DEFAULTS, Folds, search, fit_all

TARGETS = ('xtar',)
SCALES = np.array([10.])  # cm -> mm; never apply the angular mrad scale


def fit_training_design(x, xtar_truth_cm, folds, cells, *, en_config=None,
                        beam_config=None, known_offset_cm=None):
    """Fit on already-verified training rows only; return a separate model.

    x must not include a constant: centered scaling handles a free predictor
    intercept, which does not change any transport-function constant.
    Existing training folds are required, never generated here. The returned
    OOF predictions are adaptive development evidence, not an independent test.
    """
    x = np.asarray(x, dtype=float)
    target = np.asarray(xtar_truth_cm, dtype=float)
    if target.ndim != 1:
        raise ValueError('Expected one native xtarT response in cm')
    anchor = np.zeros(len(target)) if known_offset_cm is None else np.asarray(known_offset_cm, dtype=float)
    if anchor.shape != target.shape or not np.isfinite(anchor).all():
        raise ValueError('Expected one finite measured anchor in cm per training row')
    y = (target-anchor)[:, None]
    folds, cells = np.asarray(folds), np.asarray(cells)
    en = dict(EN_DEFAULTS, folds=len(set(folds)))
    if en_config is not None:
        en.update(en_config)
    validate(en)
    beam = dict(BEAM_DEFAULTS, beam=30)
    if beam_config is not None:
        beam.update(beam_config)
    if beam['beam'] != 30:
        raise ValueError('The separate xtar experiment requires beam width 30')
    selected = tune(x, y, folds, cells, en, targets=TARGETS, scales=SCALES)
    choice = selected['choices'][0]
    seeds = {(), tuple(map(int, np.flatnonzero(selected['beta'][:, 0])))}
    for f in range(en['folds']):
        support = selected['fold_supports'][0, f, choice['l1'], choice['alpha']]
        seeds.add(tuple(map(int, np.flatnonzero(support))))
    prepared = Folds(x, y, folds, cells, en['rcond'], beam['threads'],
                     targets=TARGETS, scales=SCALES)
    candidates = []
    result = search(lambda mask: prepared.evaluate(0, mask), sorted(seeds),
                    x.shape[1], beam, candidates.append)
    mask = result['chosen']['mask']
    coefficients, info = fit_all(x, y, [mask], en['rcond'], targets=TARGETS, scales=SCALES)
    oof, fold_info = prepared.predict(0, mask)
    return dict(coefficients_cm=coefficients[:, 0], selected_columns=list(mask),
                final_fit=info[0], fold_fits=fold_info, search=result,
                candidates=candidates, enet_choices=selected['choices'],
                enet_cv=selected['cv'], enet_summary=selected['summary'],
                oof_prediction_cm=oof+anchor, enet_config=en, beam_config=beam,
                uses_known_offset=known_offset_cm is not None,
                full_design_svd=selected['svd_info'], scaling=selected['state'])


def export_model(path, fitted, *, features, provenance, validity, anchor=None):
    """Export a new JSON file only. Feature expressions and units are explicit.

    features is a list of {name, expression, unit, source} dictionaries, in
    design-column order. Provenance and validity must be supplied by the
    audited data adapter. This is not a HCANA matrix format.
    """
    required = {'name', 'expression', 'unit', 'source'}
    if (not features or any(set(f) != required or not all(isinstance(v, str) and v for v in f.values())
                            for f in features) or
            len({f['name'] for f in features}) != len(features)):
        raise ValueError('Explicit distinct measured feature definitions are required')
    if not provenance or not validity:
        raise ValueError('Source hashes and validity information are required')
    if bool(anchor) != fitted.get('uses_known_offset', False):
        raise ValueError('Export anchor must match the training target offset convention')
    c = np.asarray(fitted['coefficients_cm'])
    if c.shape != (len(features)+1,) or not np.isfinite(c).all():
        raise ValueError('Export coefficient/feature mismatch')
    model = dict(schema='separate_xtar_predictor_v1', target='xtarT',
                 output_unit='cm', reporting_unit='mm', reporting_scale=10.,
                 features=features, anchor=anchor, intercept_cm=float(c[0]),
                 coefficients_cm=c[1:].tolist(), selected_columns=fitted['selected_columns'],
                 provenance=provenance, validity=validity,
                 enet_config=fitted['enet_config'], beam_config=fitted['beam_config'],
                 final_fit=fitted['final_fit'],
                 cv_status='adaptive development evidence; not an independent test',
                 replay_status='offline predictor; not installed')
    with open(path, 'x') as stream:
        json.dump(model, stream, indent=2, allow_nan=False)
        stream.write('\n')
    with open(path) as stream:
        exported = json.load(stream)
    original = c[0] + np.eye(len(features)) @ c[1:]
    reread = predict(exported, np.eye(len(features)), anchor_cm=np.zeros(len(features)) if anchor else None)
    np.testing.assert_allclose(reread, original, rtol=1e-14, atol=1e-14)
    return exported


def predict(model, measured_design, *, anchor_cm=None):
    """Consume the documented measured design; produce one position in cm."""
    if model.get('schema') != 'separate_xtar_predictor_v1' or model.get('output_unit') != 'cm':
        raise ValueError('Not a separate cm-valued xtar model')
    x = np.asarray(measured_design, dtype=float)
    c = np.asarray(model['coefficients_cm'], dtype=float)
    if (x.ndim != 2 or x.shape[1] != len(model['features']) or c.shape != (x.shape[1],)
            or not np.isfinite(x).all() or not np.isfinite(c).all()
            or not np.isfinite(model['intercept_cm'])):
        raise ValueError('Prediction design/coefficients invalid')
    if model.get('anchor'):
        if anchor_cm is None:
            raise ValueError('This model requires its measured beam/pointing anchor')
        anchor = np.asarray(anchor_cm, dtype=float)
        if anchor.shape != (len(x),) or not np.isfinite(anchor).all():
            raise ValueError('Invalid measured anchor')
    else:
        if anchor_cm is not None:
            raise ValueError('Unanchored model must not receive an extra anchor')
        anchor = 0.
    return anchor + model['intercept_cm'] + x @ c
