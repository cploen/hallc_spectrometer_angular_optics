"""Observed-density and component diagnostics; no fitted PDF used as data."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from .x_geometry import global_centers


def render(path, axis, raw, prediction, polygons, model, title, data, learning=None):
    path = Path(path)
    if path.exists() or path.with_suffix('.png').exists():
        raise FileExistsError(str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    labels = prediction['labels']
    q, p = prediction['q'], prediction['p']
    phi = 0.
    unfolded = axis == 'x' and not model.get('manual_partition')
    if unfolded:
        fan = model['geometry']
        middle = float(np.median(p))
        step = max(float(np.ptp(p))*1e-4, 1e-6)
        centers = global_centers(fan, [middle-step, middle+step])
        slope = np.median((centers[:, 1]-centers[:, 0])/(2*step))
        phi = -float(np.arctan(slope))
    c, s = np.cos(phi), np.sin(phi)
    view_q, view_p = c*q+s*p, -s*q+c*p
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
        _density(ax[0, 0], view_q, view_p)
        ax[0, 0].set(xlabel='standardized q (display rotation)', ylabel='standardized p (display rotation)',
                     title='Observed density and exported dividers')
        transform = model['geometry'] if axis == 'x' else model['model']
        theta = np.deg2rad(transform['theta_deg'])
        for poly in polygons:
            z = (np.asarray(poly['points'])-transform['mean'])/transform['scale']
            pq = z @ [np.cos(theta), np.sin(theta)]
            pp = z @ [-np.sin(theta), np.cos(theta)]
            ax[0, 0].plot(c*pq+s*pp, -s*pq+c*pp, color='crimson', lw=.6)
        ax[1, 0].hist(view_q, bins=240, histtype='step', color='black')
        ax[1, 0].set(xlabel='display q', ylabel='observed events', title='Actual 1D projection; curvature remains visible')
        if unfolded:
            u = prediction['u']; valid = np.isfinite(u)
            _density(ax[0, 1], u[valid], p[valid])
            bins = np.linspace(-.5, len(transform['possible_indices'])-.5, 241)
            ax[1, 1].hist(u[valid], bins=bins, histtype='step', color='black', label='target slice')
            if learning is not None:
                h, _ = np.histogram(learning[np.isfinite(learning)], bins=bins)
                if h.sum():
                    ax[1, 1].step((bins[:-1]+bins[1:])/2, h*valid.sum()/h.sum(), color='steelblue', alpha=.7, label='geometry window (normalized)')
            ids = np.asarray(model['occupancy']['accepted_indices'])
            dividers = model.get('boundaries_u', ((ids[1:]+ids[:-1])/2).tolist())
            for divider in dividers:
                for a in ax[:, 1]: a.axvline(divider, color='crimson', lw=.7)
            ax[0, 1].set(xlabel='unfolded band coordinate u', ylabel='model p', title='Nonlinear unfolding (model-dependent)')
            ax[1, 1].set(xlabel='unfolded band coordinate u', ylabel='observed events', title='Actual unfolded histogram')
            ax[1, 1].legend(fontsize=8)
        else:
            for label in np.unique(labels[labels >= 0]):
                take = labels == label
                ax[0, 1].scatter(q[take], p[take], s=.25, label=str(label), rasterized=True)
            ax[0, 1].set(xlabel='model q', ylabel='model p', title=('Manual cut assignments' if model.get('manual_partition') else 'Gaussian decision assignments')+' (provisional IDs)')
            ax[1, 1].hist(q, bins=240, histtype='step', color='black')
            fitted = model.get('model', {})
            for mean in fitted.get('means', []):
                ax[1, 1].axvline(mean[0], lw=.7, color='steelblue')
            ax[1, 1].set(xlabel='model q', ylabel='observed events', title='Observed projection'+(' and fitted component means' if fitted.get('means') else ''))
        fig.suptitle('PROVISIONAL '+title+'\nRotation = %.3f deg; display = %.3f deg' %
                     (transform['theta_deg'], transform['theta_deg']+np.rad2deg(phi)))
        fig.savefig(path.with_suffix('.png'), dpi=130)
        pdf.savefig(fig); plt.close(fig)
        fig, ax = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        _density(ax[0], raw[:, 0], raw[:, 1])
        for polygon in polygons:
            points = np.asarray(polygon['points']); ax[0].plot(points[:, 0], points[:, 1], lw=.65)
        ax[0].set(xlabel=axis+'pfp (rad)', ylabel=axis+'fp (cm)', title='Raw coordinates and editable TCutG boundaries')
        if 'xsieve' in data and 'ysieve' in data:
            for label in np.unique(labels):
                take = labels == label
                ax[1].scatter(data['ysieve'][take], data['xsieve'][take], s=.3, rasterized=True, label=str(label))
            ax[1].set(xlabel='ysieve (cm)', ylabel='xsieve (cm)', title='Preliminary components at sieve; physical labels in step 3')
        fig.suptitle(title); pdf.savefig(fig); plt.close(fig)
        if axis == 'y' and model['model'].get('scan_summary'):
            scan = model['model']['scan_summary']
            fig, ax = plt.subplots(2, 1, figsize=(10, 6), constrained_layout=True)
            ax[0].plot([r['angle'] for r in scan], [r['score'] for r in scan], '.-')
            ax[1].plot([r['angle'] for r in scan], [r['bands'] for r in scan], '.-')
            ax[0].set(ylabel='spacing score'); ax[1].set(xlabel='rotation angle (deg)', ylabel='supported seed bands')
            fig.suptitle(title+' — Y angle scan'); pdf.savefig(fig); plt.close(fig)
        for label in np.unique(labels[labels >= 0]):
            take = labels == label
            fig, ax = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
            _density(ax[0], raw[take, 0], raw[take, 1]); ax[0].set(xlabel=axis+'pfp', ylabel=axis+'fp')
            ax[1].hist(q[take], bins=100, histtype='step'); ax[1].set(xlabel='model q', ylabel='observed events')
            if 'xsieve' in data and 'ysieve' in data:
                _density(ax[2], data['ysieve'][take], data['xsieve'][take]); ax[2].set(xlabel='ysieve (cm)', ylabel='xsieve (cm)')
            fig.suptitle(title+' — provisional component %d, N=%d' % (label, take.sum()))
            pdf.savefig(fig); plt.close(fig)


def _density(ax, x, y):
    if len(x):
        ax.hist2d(x, y, bins=150, cmin=1, cmap='viridis')
        if len(x) > 200:
            ax.set_xlim(np.quantile(x, [.001, .999])); ax.set_ylim(np.quantile(y, [.001, .999]))
