#!/usr/bin/env python3
"""Campaign-independent sieve slit scattering diagnostic of saved replay ytar."""
import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import uproot

PROJECT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT))
from spectrometer_config import from_campaign


def populations(data):
    """Rank within foil/delta/hole; ties use original entry, never ytar."""
    quality = data['quality']
    if not np.isin(quality, [0, 1, 2]).all():
        raise ValueError('Unknown quality code')
    core = quality == 2
    if not np.isfinite(data['core_score'][core]).all():
        raise ValueError('Nonfinite core score')
    dense = np.zeros(len(core), dtype=bool)
    keys = np.column_stack([data[k] for k in ('foil', 'ndel', 'xscol', 'yscol')])
    _, inverse = np.unique(keys, axis=0, return_inverse=True)
    for group in np.unique(inverse[core]):
        idx = np.flatnonzero(core & (inverse == group))
        order = np.lexsort((data['entry'][idx], -data['core_score'][idx]))
        dense[idx[order[:(len(idx)+1)//2]]] = True
    return dict(dense_half=dense, outer_core=core & ~dense,
                shoulders=quality == 1, unsupported=quality == 0,
                out_of_core=quality != 2)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_table(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--config', type=Path)
    parser.add_argument('--rungroup', action='append')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    campaign = args.campaign.resolve() if args.campaign.is_dir() else (PROJECT / args.campaign).resolve()
    spec = from_campaign(campaign)
    config = args.config or campaign / 'config/sieve_slit.json'
    cfg = json.loads(config.read_text())
    def resolve(value):
        path = Path(value)
        return path if path.is_absolute() else campaign / path
    table = resolve(cfg['rungroup_table'])
    with table.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if args.rungroup:
        if set(args.rungroup) - {r['rungroup'] for r in rows}:
            raise ValueError('Unknown rungroup')
        rows = [r for r in rows if r['rungroup'] in args.rungroup]
    bins = cfg.get('bins', 100)
    if type(bins) is not int or bins < 2:
        raise ValueError('bins must be an integer >= 2')
    branches = ['entry', 'run', 'foil', 'zfoil', 'ndel', 'xscol', 'yscol', 'quality', 'core_score', 'ytar']
    jobs = []
    for row in rows:
        path = resolve(cfg['core_file'].format(rungroup=row['rungroup']))
        with uproot.open(path) as root:
            data = root['CoreSample'].arrays(branches, library='np')
        if not len(data['entry']) or len(np.unique(data['entry'])) != len(data['entry']):
            raise ValueError(f'{path}: empty input or duplicate entries')
        if not np.all(data['run'] == int(row['optics_id'])):
            raise ValueError(f'{path}: optics ID mismatch')
        if not np.isfinite(data['ytar']).all() or not np.isfinite(data['zfoil']).all():
            raise ValueError(f'{path}: nonfinite coordinates')
        masks = populations(data)
        jobs.append((row, path, data, masks))
        print(f"{spec.name} {row['rungroup']}: {len(data['entry'])} events; {masks['out_of_core'].sum()} out of core", flush=True)
    if args.check:
        return
    out = resolve(cfg.get('output', '07_diagnostics/sieve_slit'))
    out.mkdir(parents=True, exist_ok=True)
    summary, histogram = [], []
    manifest = dict(status='running', spectrometer=spec.name, config=cfg,
                    config_sha256=digest(config), table_sha256=digest(table),
                    script_sha256=digest(Path(__file__)), inputs=[],
                    versions={p: importlib.metadata.version(p) for p in ['numpy', 'matplotlib', 'uproot']})
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    labels = dict(dense_half='Densest 50% of cores', outer_core='Remaining core events',
                  shoulders='Shoulders / between cores', unsupported='Unsupported holes',
                  out_of_core='All out-of-core events')
    colors = dict(dense_half='#2166ac', shoulders='#d6604d', out_of_core='#762a83')
    for row, path, data, masks in jobs:
        name = row['rungroup']
        manifest['inputs'].append(dict(rungroup=name, path=cfg['core_file'].format(rungroup=name), sha256=digest(path)))
        for foil in np.unique(data['foil']):
            selected = data['foil'] == foil
            z = np.unique(data['zfoil'][selected])
            if len(z) != 1:
                raise ValueError('Foil label has multiple positions')
            values = data['ytar'][selected]
            low, high = float(values.min()), float(values.max())
            margin = max((high-low)*.02, .01)
            edges = np.linspace(low-margin, high+margin, bins+1)
            stats = {}
            for key, mask in masks.items():
                y = data['ytar'][selected & mask]
                counts, _ = np.histogram(y, edges)
                stats[key] = dict(rungroup=name, zfoil_cm=float(z[0]), population=key, n=len(y),
                                  mean_cm=float(np.mean(y)) if len(y) else '',
                                  std_cm=float(np.std(y)) if len(y) else '',
                                  median_cm=float(np.median(y)) if len(y) else '',
                                  p05_cm=float(np.quantile(y, .05)) if len(y) else '',
                                  p95_cm=float(np.quantile(y, .95)) if len(y) else '')
                summary.append(stats[key])
                for a, b, count in zip(edges[:-1], edges[1:], counts):
                    histogram.append(dict(rungroup=name, zfoil_cm=float(z[0]), population=key,
                                          low_cm=a, high_cm=b, count=int(count)))
            for key, color in colors.items():
                y = data['ytar'][selected & masks[key]]
                fig, ax = plt.subplots(figsize=(9, 5.5), layout='constrained')
                ax.hist(y, bins=edges, histtype='stepfilled', alpha=.65, color=color)
                ax.set(xlabel=r'$y_{\mathrm{tar}}$ (cm)', ylabel='Events / bin', xlim=(edges[0], edges[-1]))
                fig.suptitle(f'Sieve slit scattering study · {spec.name}', fontsize=16)
                ax.set_title(f'{labels[key]} · foil {z[0]:g} cm\n{name} · all holes and δ slices', fontsize=11)
                text = f'N = {len(y):,}'
                if len(y):
                    text += f"\nMean = {np.mean(y):.3f} cm\nSD = {np.std(y):.3f} cm"
                ax.text(.98, .96, text, transform=ax.transAxes, ha='right', va='top', fontsize=10)
                ax.grid(alpha=.2)
                fig.savefig(out/f'{name}_foil{int(foil)}_{key}.png', dpi=160)
                plt.close(fig)
    write_table(out/'summary.tsv', summary)
    write_table(out/'histograms.tsv', histogram)
    manifest['status'] = 'complete'
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    main()
