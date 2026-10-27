#!/usr/bin/env python3
"""Render frozen pilot comparisons; no fitting or membership changes.

The original pilot data are not distributed with this branch. Supply the
project directory explicitly to reproduce these documentation figures.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path(__file__).parents[1] / 'joint_surface_figures')
    args = parser.parse_args()
    study = args.project / 'outputs/provisional_xy_joint_surfaces_20261027'
    dest = args.output
    dest.mkdir(parents=True, exist_ok=True)
    hashes = {}

    def load(path):
        hashes[str(path.relative_to(args.project))] = hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path, allow_pickle=False) as data:
            return {k: data[k].copy() for k in data.files}

    def events(case):
        return load(study / 'postfit/evaluation/joined_predictions' / (case + '_A_min10_test.npz'))

    def save(fig, name):
        fig.savefig(dest / (name + '.png'), dpi=190, bbox_inches='tight', facecolor='white')
        plt.close(fig)

    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9})
    cmap = plt.get_cmap('tab10')
    norm = BoundaryNorm(np.arange(-.5, 9.5), cmap.N)
    rng = np.random.default_rng(20261027)
    central = events('central1544')
    other = events('rg07_minus8')
    # Include all rejections, and randomly sample assigned events. This display
    # emphasizes rejections and is explicitly not a density estimate.
    fig = plt.figure(figsize=(10.6, 4.4))
    for i, axis in enumerate(['x', 'y']):
        ax = fig.add_subplot(1, 2, i + 1, projection='3d')
        accepted = central['selected_' + axis]
        keep = np.flatnonzero(accepted)
        if len(keep) > 4800:
            keep = rng.choice(keep, 4800, replace=False)
        ax.scatter(central[axis+'fp'][keep], central[axis+'pfp'][keep]*1000, central['delta'][keep],
                   c=central['physical_'+axis][keep], cmap=cmap, norm=norm, s=2.5, alpha=.55, rasterized=True)
        reject = ~accepted
        ax.scatter(central[axis+'fp'][reject], central[axis+'pfp'][reject]*1000, central['delta'][reject],
                   color='#666666', marker='x', s=8, alpha=.5, rasterized=True)
        ax.set(xlabel=axis+'fp (cm)', ylabel=axis+"'fp (mrad)", zlabel='delta (%)',
               title=f'{axis.upper()} plane: {accepted.sum():,}/{len(accepted):,} assigned')
        ax.view_init(elev=21, azim=-60 if axis=='x' else -66)
    fig.tight_layout(w_pad=2)
    save(fig, 'central_3d')

    fig, axs = plt.subplots(2, 2, figsize=(10.5, 6.6))
    for row, (label, data) in enumerate([('Run 1544, foil 0 cm', central), ('Runs 6983/6984, foil -8 cm', other)]):
        use = (data['delta'] >= -10) & (data['delta'] < -8)
        for col, axis in enumerate(['x', 'y']):
            ax = axs[row, col]
            accepted = use & data['selected_'+axis]
            reject = use & ~data['selected_'+axis]
            ax.scatter(data[axis+'fp'][accepted], data[axis+'pfp'][accepted]*1000,
                       c=data['physical_'+axis][accepted], cmap=cmap, norm=norm, s=3, alpha=.65, rasterized=True)
            ax.scatter(data[axis+'fp'][reject], data[axis+'pfp'][reject]*1000, c='#666666', marker='x', s=10, alpha=.6)
            ax.set(xlabel=axis+'fp (cm)', ylabel=axis+"'fp (mrad)", title=f'{label}: {axis.upper()}, {accepted.sum():,}/{use.sum():,}')
            ax.grid(alpha=.2)
    fig.tight_layout()
    save(fig, 'low_delta_planes')

    co = load(study / 'postfit/coordinates/central1544.npz')
    assert np.array_equal(co['entry'][central['indices']], central['entry'])
    xs, ys = co['xsieve'][central['indices']], co['ysieve'][central['indices']]
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for ax, bounds in zip(axs, [(-10, -8), (0, 5)]):
        use = (central['delta'] >= bounds[0]) & (central['delta'] < bounds[1])
        accepted = use & central['selected_xy']
        reject = use & ~central['selected_xy']
        ax.scatter(xs[accepted], ys[accepted], c=central['physical_x'][accepted], cmap=cmap, norm=norm, s=2, alpha=.45, rasterized=True)
        ax.scatter(xs[reject], ys[reject], color='#666666', marker='x', s=10)
        for x in (np.arange(9)-4)*2.54:
            ax.axvline(x, color='#aaaaaa', lw=.45, zorder=0)
        for y in (np.arange(9)-4)*1.524:
            ax.axhline(y, color='#aaaaaa', lw=.45, zorder=0)
        ax.set(xlabel='Reconstructed x_sieve (cm)', ylabel='Reconstructed y_sieve (cm)',
               title=f'{bounds[0]}% <= delta < {bounds[1]}%: {accepted.sum():,} pairs', xlim=(-12.5,12.5), ylim=(-8,8))
        ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    save(fig, 'sieve_plane')

    registry = json.loads((study/'case_registry.json').read_text())
    model_path = args.project / registry['central1544']['models']['x']['A']['10']['model']
    hashes[str(model_path.relative_to(args.project))] = hashlib.sha256(model_path.read_bytes()).hexdigest()
    model = json.loads(model_path.read_text())['model']
    raw = np.c_[central['xpfp'], central['xfp']]
    z = (raw-model['mean'])/model['scale']
    angle = np.deg2rad(model['theta_deg'])
    q = z @ np.array([np.cos(angle),np.sin(angle)])
    r = z @ np.array([-np.sin(angle),np.cos(angle)])
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.6))
    for ax, coord, name in zip(axs, [q,r], ['q: standardized separation','r: perpendicular standardized coordinate']):
        for low, high, color in [(-10,-8,'#0072b2'),(0,5,'#d55e00')]:
            use=(central['delta']>=low)&(central['delta']<high)
            ax.hist(coord[use], bins=100, histtype='step', density=True, color=color,
                    label=f'[{low},{high})%: {use.sum():,} events')
        ax.set(xlabel=name, ylabel='Observed projection density')
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    save(fig, 'actual_projections')

    # This physical-coordinate plot shows the narrow correlated acceptance
    # domain; it is not a proof of separately identifiable delta coefficients.
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.6))
    for ax, data, label in zip(axs,[central,other],['Run 1544, foil 0 cm','Runs 6983/6984, foil -8 cm']):
        ax.scatter(data['delta'], data['xfp'], c=data['physical_x'], cmap=cmap, norm=norm,
                   s=1.8, alpha=.25, rasterized=True)
        ax.set(xlabel='Measured delta (%)', ylabel='xfp (cm)', title=label)
        ax.grid(alpha=.2)
    fig.tight_layout()
    save(fig, 'acceptance_domain')

    # Mathematics is rendered separately from data, with symbol definitions in
    # the document. These are equations, not additional fitted curves.
    equations = {
        'rotation': r'$z_1=(s-\mu_s)/\sigma_s,\quad z_2=(x-\mu_x)/\sigma_x,\qquad q=z_1\cos\theta+z_2\sin\theta,\quad r=-z_1\sin\theta+z_2\cos\theta$',
        'basis': r'$F(t,R,D)=a_0+a_1R+a_2t+a_3tR+a_4R^2+a_5tR^2+a_6t^2+a_7t^2R+a_8t^2R^2+D(b_0+b_1t+b_2R+b_3tR)$',
        'objective': r'$\widehat{c}=\arg\min_c\left[\sum_i w_i\rho\left(q_i-\phi(j_i,r_i,\delta_i)^Tc\right)+c^T\Lambda c\right],\qquad \partial F/\partial j\geq m>0$',
        'phase': r'$C(b)=\left|\sum_i w_i\exp\{2\pi\,\mathrm{i}\,\psi_b(q_i,r_i,\delta_i)\}\right|^2,\qquad j_i\simeq\operatorname{round}(\psi_b-\alpha)$',
        'classify': r'$\widehat{j}_i=\arg\min_{j\in\mathcal{O}(\delta_i)}|q_i-F(j,r_i,\delta_i)|,\qquad (\mathrm{X\ row},\mathrm{Y\ column})=(s_Xj_X+o_X,\ s_Yj_Y+o_Y)$',
        'boundary': r'$q_{j+1/2}(r,\delta)=\frac{1}{2}\{F(j,r,\delta)+F(j+1,r,\delta)\}$',
    }
    for name, text in equations.items():
        fig=plt.figure(figsize=(11, .65));fig.text(.01,.5,text,va='center',fontsize=12);save(fig,'equation_'+name)
    comparison=study/'plots/worst_low_delta_x_3d.png'
    hashes[str(comparison.relative_to(args.project))]=hashlib.sha256(comparison.read_bytes()).hexdigest()
    shutil.copyfile(comparison,dest/'worst_low_delta_x_3d.png')
    (dest/'figure_manifest.json').write_text(json.dumps({'source_sha256':hashes,
        'fit':'Frozen pilot A, minimum count 10', 'coordinates':'Measured focal plane; separately identified reconstructed sieve plane',
        'fitting_performed':False, 'central_3d_sampling':'4800 random accepted events per axis plus all rejections; not a density estimate',
        'other_plots':'All common evaluation events in indicated intervals; no thinning'},indent=2)+'\n')


if __name__=='__main__':
    main()
