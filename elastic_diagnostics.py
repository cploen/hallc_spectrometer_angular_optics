"""Compact angular-fit comparison: residuals, coverage, stability and event clouds."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from elastic_net import TARGETS, SCALES

LABELS = ("x′tar (mrad)", 'ytar (cm)', "y′tar (mrad)", 'ztar (cm)')


def residuals(a, truth, predicted, offset):
    r = (predicted-truth)*SCALES
    total = predicted+offset
    s, c = np.sin(a['angle']), np.cos(a['angle'])
    deg = np.rad2deg(a['angle'])
    ymis = .1*(.52-.012*np.abs(deg)+.002*deg**2)
    beam = (a['ytarT']+ymis-a['ztarT']*(s-a['yptarT']*c))/(c+a['yptarT']*s)
    denom = s-total[:, 2]*c
    if np.any(np.abs(denom) < 1e-10):
        raise ValueError('Near-singular reconstructed ztar; inspect candidate before proceeding')
    z = (100*total[:, 1]+ymis-beam*(c+total[:, 2]*s))/denom
    if not np.isfinite(z).all(): raise ValueError('Nonfinite ztar reconstruction')
    return np.column_stack((r, z-a['ztarT']))


def groups(a, pool):
    for p in sorted(set(pool)-{'blocked'}):
        pm = pool == p
        yield dict(pool=p, level='overall', zfoil='', ndel='', rungroup='', xscol='', yscol=''), pm
        for z in np.unique(a['ztarT']):
            for d in np.unique(a['ndel']):
                m = pm & (a['ztarT'] == z) & (a['ndel'] == d)
                if not m.any(): continue
                base = dict(pool=p, level='foil_delta', zfoil=float(z), ndel=int(d), rungroup='', xscol='', yscol='')
                yield base, m
                for rg in np.unique(a['rungroup'][m]):
                    sm = m & (a['rungroup'] == rg)
                    yield dict(base, level='setting', rungroup=rg), sm
                    for x, y in sorted(set(zip(a['xscol'][sm], a['yscol'][sm]))):
                        yield dict(base, level='hole', rungroup=rg, xscol=int(x), yscol=int(y)), sm & (a['xscol'] == x) & (a['yscol'] == y)


def metrics(a, truth, predictions, pool, qa_min, offset):
    rr = {k: residuals(a, truth, pred, offset) for k, pred in predictions.items()}
    rows = []
    for g, mask in groups(a, pool):
        for model, values in rr.items():
            for j, target in enumerate((*TARGETS, 'ztar')):
                v = values[mask, j]
                rows.append(dict(g, model=model, target=target, n=len(v), qa_low=len(v) < qa_min,
                                 bias=float(v.mean()), rms=float(np.sqrt(np.mean(v*v))),
                                 p95=float(np.percentile(np.abs(v), 95))))
    return rows


def save(fig, path):
    fig.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(fig)


def training_plots(out, result):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for j, (name, ax) in enumerate(zip(TARGETS, axes)):
        rows = [r for r in result['summary'] if r['target'] == name]
        for l1 in sorted({r['l1'] for r in rows}):
            r = sorted([r for r in rows if r['l1'] == l1 and r['converged']], key=lambda v: v['alpha'])
            if r:
                ax.errorbar([v['alpha'] for v in r], [v['mse']*SCALES[j]**2 for v in r],
                            yerr=[v['se']*SCALES[j]**2 for v in r], label=f'L1 {l1:g}', marker='.', capsize=2)
        c = result['choices'][j]
        ax.scatter(c['alpha'], c['mse']*SCALES[j]**2, marker='*', s=100, c='black', label='Chosen')
        ax.axhline(c['svd_cv_mse']*SCALES[j]**2, color='gray', ls='--', label='Scaled SVD')
        ax.set(xscale='log', yscale='log', title=name, xlabel='Alpha', ylabel='Mean cell MSE (physical units²)')
        ax.legend(fontsize=8)
    save(fig, out/'plots/cv.png')
    fig, axes = plt.subplots(2, 1, figsize=(13, 5), constrained_layout=True)
    for j, name in enumerate(TARGETS):
        axes[0].plot(np.arange(len(result['beta']))+1, np.abs(result['beta'][:, j]), '.', label=name)
        axes[1].plot(np.arange(len(result['frequency']))+1, result['frequency'][:, j], '.', label=name)
    axes[0].set_yscale('symlog', linthresh=1e-5)
    axes[0].set_ylabel('|Standardized coefficient|')
    axes[1].set(xlabel='Nonconstant term index (seed order)', ylabel='Fraction of folds selected', ylim=(-.05, 1.05))
    for ax in axes: ax.legend()
    save(fig, out/'plots/terms.png')


def evaluation_plots(out, a, truth, predictions, offset, pool, spec, qa_min):
    rr = {k: residuals(a, truth, pred, offset) for k, pred in predictions.items()}
    pools = ('core', 'noncore', 'unsupported')
    fig, axes = plt.subplots(3, 4, figsize=(15, 10), constrained_layout=True)
    for i, p in enumerate(pools):
        m = pool == p
        for j, ax in enumerate(axes[i]):
            if not m.any():
                ax.text(.5, .5, 'No events', transform=ax.transAxes, ha='center'); continue
            # Full-range histograms: no hidden tail or silent plotting truncation.
            values = np.concatenate([r[m, j] for r in rr.values()])
            bins = np.histogram_bin_edges(values, bins=100)
            for model, r in rr.items(): ax.hist(r[m, j], bins=bins, histtype='step', label=model, log=True)
            ax.set(title=f'{p}, N={sum(m):,}', xlabel='Residual '+LABELS[j], ylabel='Events / bin'); ax.legend()
    save(fig, out/'plots/residuals.png')
    foils, deltas = np.unique(a['ztarT']), np.unique(a['ndel'])
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), constrained_layout=True)
    for i, p in enumerate(pools[:2]):
        for j, ax in enumerate(axes[i]):
            values = np.full((len(foils), len(deltas)), np.nan)
            for zi, z in enumerate(foils):
                for di, d in enumerate(deltas):
                    m = (pool == p) & (a['ztarT'] == z) & (a['ndel'] == d)
                    if m.any():
                        ref = np.sqrt(np.mean(rr['svd'][m, j]**2))
                        val = np.sqrt(np.mean(rr['enet'][m, j]**2))
                        values[zi, di] = val/max(ref, 1e-30)
                        ax.text(di, zi, f'{values[zi,di]:.2f}\nN={sum(m)}', ha='center', va='center', fontsize=7)
            im = ax.imshow(values, cmap='coolwarm', norm=TwoSlopeNorm(vmin=.5, vcenter=1, vmax=1.5), aspect='auto')
            ax.set(xticks=range(len(deltas)), xticklabels=deltas, yticks=range(len(foils)), yticklabels=foils,
                   xlabel='Delta slice', ylabel='Foil (cm)', title=f'{p}: {LABELS[j]}')
    fig.colorbar(im, ax=axes.ravel().tolist(), label='Elastic / SVD RMS (below 1 is lower)', extend='both', shrink=.7)
    save(fig, out/'plots/foil_delta.png')
    for d in deltas:
        fig, axes = plt.subplots(len(foils), 4, figsize=(16, 3.3*len(foils)), squeeze=False, constrained_layout=True)
        maxcount = 1
        cells = []
        for zi, z in enumerate(foils):
            for pi, p in enumerate(pools[:2]):
                m = (pool == p) & (a['ztarT'] == z) & (a['ndel'] == d)
                for x, y in sorted(set(zip(a['xscol'][m], a['yscol'][m]))):
                    h = m & (a['xscol'] == x) & (a['yscol'] == y)
                    n = int(sum(h)); maxcount = max(maxcount, n)
                    ref = np.sqrt(np.mean(rr['svd'][h][:, [0, 2]]**2))
                    val = np.sqrt(np.mean(rr['enet'][h][:, [0, 2]]**2))
                    cells.append((zi, pi, spec.ys(int(y)), spec.xs(int(x)), n, val/max(ref, 1e-30)))
        for zi, pi, y, x, n, ratio in cells:
            axes[zi, pi].scatter(y, x, c=[n], norm=LogNorm(1, max(2, maxcount)), cmap='viridis', s=260)
            axes[zi, pi].text(y, x, str(n), ha='center', va='center', fontsize=7, color='black' if LogNorm(1, max(2, maxcount))(n) > .65 else 'white')
            axes[zi, pi+2].scatter(y, x, c=[ratio], norm=TwoSlopeNorm(vmin=.5, vcenter=1, vmax=1.5), cmap='coolwarm', s=260,
                                  edgecolors='black' if n < qa_min else 'none')
            axes[zi, pi+2].text(y, x, f'{ratio:.2f}', ha='center', va='center', fontsize=7, color='white' if abs(ratio-1) > .3 else 'black')
        for zi, z in enumerate(foils):
            for j, ax in enumerate(axes[zi]):
                ax.title.set_fontsize(9)
                ax.set(xlim=(spec.ys(0)-1, spec.ys(spec.ny-1)+1), ylim=(spec.xs(0)-1, spec.xs(spec.nx-1)+1),
                       xlabel='Y sieve (cm)', ylabel='X sieve (cm)', aspect='equal',
                       title=f'Foil {z:g}: '+('core N', 'noncore N', 'core RMS ratio', 'noncore RMS ratio')[j])
        fig.suptitle(f'Delta slice {d}: counts use log color; RMS ratio = elastic / SVD\nBlack circle outlines: fewer than {qa_min} events; blank = no events in that pool', fontsize=13)
        fig.colorbar(plt.cm.ScalarMappable(norm=LogNorm(1, max(2, maxcount)), cmap='viridis'), ax=axes[:, :2].ravel().tolist(), shrink=.5, label='Events')
        fig.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(vmin=.5, vcenter=1, vmax=1.5), cmap='coolwarm'), ax=axes[:, 2:].ravel().tolist(), shrink=.5, label='Angular RMS ratio', extend='both')
        save(fig, out/f'plots/slice_{d}_holes.png')
        for p in pools[:2]:
            fig, axes = plt.subplots(len(foils), 3, figsize=(11, 3*len(foils)), squeeze=False, constrained_layout=True)
            xy = [(a['ys'], a['xs'])]
            for model in ('svd', 'enet'):
                total = predictions[model]+offset
                xy.append((100*total[:, 1]+spec.sieve_distance*total[:, 2], a['xtar']+spec.sieve_distance*total[:, 0]))
            m = (pool == p) & (a['ndel'] == d)
            if m.any():
                allx = np.concatenate([x[m] for y, x in xy]); ally = np.concatenate([y[m] for y, x in xy])
                limits = ((min(spec.ys(0)-1, ally.min()), max(spec.ys(spec.ny-1)+1, ally.max())),
                          (min(spec.xs(0)-1, allx.min()), max(spec.xs(spec.nx-1)+1, allx.max())))
            else: limits = ((spec.ys(0)-1, spec.ys(spec.ny-1)+1), (spec.xs(0)-1, spec.xs(spec.nx-1)+1))
            for zi, z in enumerate(foils):
                h = m & (a['ztarT'] == z)
                for j, ((y, x), ax) in enumerate(zip(xy, axes[zi])):
                    ax.scatter(y[h], x[h], s=.4, alpha=.45, rasterized=True)
                    ax.set(xlim=limits[0], ylim=limits[1], aspect='equal', xlabel='Y sieve (cm)', ylabel='X sieve (cm)',
                           title=f'Foil {z:g}, N={sum(h):,}: '+('saved replay', 'scaled SVD', 'elastic net')[j])
            fig.suptitle(f'Protected {p}, delta slice {d}: identical events in every comparison')
            save(fig, out/f'plots/slice_{d}_{p}_clouds.png')
