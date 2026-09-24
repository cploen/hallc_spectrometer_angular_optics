#!/usr/bin/env python3
"""Sieve clouds before/after the frozen Huber trial; reuses saved residuals."""
import argparse
import html
import json
import shutil
import zipfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from fit_huber import load_frozen
from reallocate_core import load_inputs
from core_sample import digest, write_tsv
from spectrometer_config import from_campaign

MODELS=('balanced_squared_full','all_supported_huber_full')
TITLES=('Before: balanced SVD','After: expanded smooth Huber')


def foil_name(z):
    return f'foil_{"m" if z<0 else "p" if z>0 else ""}{abs(z):g}cm'


def draw(axes, histograms, z, d, yedges, xedges, norm, counts, title_prefix=''):
    for j,ax in enumerate(axes):
        h=histograms[(z,d,j)]
        mesh=ax.pcolormesh(yedges,xedges,np.ma.masked_equal(h.T,0),cmap='viridis',norm=norm,rasterized=True)
        n=counts[(z,d)]
        ax.set(title=f'{title_prefix}{TITLES[j]}\nN = {n:,}; outside frame = {n-int(h.sum())}',
               xlabel='Y sieve (cm)',ylabel='X sieve (cm)',xlim=(-7,7),ylim=(-13,13),aspect='equal')
        ax.tick_params(labelsize=9)
    return mesh


def overviews(output,foils,slices,edges,hists,yedges,xedges,norm,counts,footer):
    for d in slices:
        fig,axes=plt.subplots(len(foils),2,figsize=(10.5,24))
        fig.subplots_adjust(left=.08,right=.84,top=.945,bottom=.06,hspace=.47,wspace=.28)
        for i,z in enumerate(foils):
            mesh=draw(axes[i],hists,z,d,yedges,xedges,norm,counts,title_prefix=f'Foil {z:+g} cm\n')
            for ax in axes[i]: ax.title.set_fontsize(10)
        colorax=fig.add_axes([.88,.35,.018,.3])
        fig.colorbar(mesh,cax=colorax,label='Events per bin (shared log scale)')
        label=f'Delta slice {d}: [{edges[d]:g}, {edges[d+1]:g})%' if d>=0 else f'All delta slices: [{edges[0]:g}, {edges[-1]:g})%'
        fig.suptitle(f'HMS 6.667 GeV · every foil\n{label}',fontsize=15,y=.985)
        fig.text(.5,.018,footer,ha='center',fontsize=9)
        tag=f'delta_{d}' if d>=0 else 'all_delta'
        fig.savefig(output/'plots'/f'all_foils_{tag}.png',dpi=140,facecolor='white');plt.close(fig)


def run(campaign,study,output):
    if output.exists(): raise FileExistsError(output)
    manifest=json.loads((study/'manifest.json').read_text())
    if manifest['status']!='complete': raise ValueError('Study is not complete')
    saved=study/'residuals.npz'
    if digest(saved)!=manifest['outputs']['residuals.npz']: raise ValueError('Saved residual checksum changed')
    a,inputs,_=load_frozen(campaign)
    keep=(a['sample']==2)&~a['excluded']&(a['quality']>=1)
    a={k:v[keep] for k,v in a.items()}
    length=from_campaign(campaign).sieve_distance
    with np.load(saved) as data:
        for k in ('entry','rungroup','quality','ndel'):
            if not np.array_equal(a[k],data[k]): raise ValueError(f'Reserved identity mismatch: {k}')
        if not np.array_equal(a['zfoil'],data['zfoil']): raise ValueError('Foil mismatch')
        xy=[]
        for model in MODELS:
            residual=data[model]
            # Saved residual = reconstructed target - existing truth target.
            # Match elastic_diagnostics.evaluation_plots exactly, with saved xtar.
            xs=a['xtar']+length*(a['xptarT']+residual[:,0]/1000.)
            ys=a['ytarT']+residual[:,1]+length*(a['yptarT']+residual[:,2]/1000.)
            if not np.isfinite(xs).all() or not np.isfinite(ys).all(): raise ValueError('Nonfinite projection')
            xy.append((ys,xs))
    edges=load_inputs(campaign,'min10')['edges']
    foils=sorted(np.unique(a['zfoil']).tolist())
    slices=list(range(len(edges)-1))+[-1]
    yedges=np.linspace(-7,7,281);xedges=np.linspace(-13,13,326)
    hists,counts,records={},{},[]
    for z in foils:
        for d in slices:
            sel=(a['zfoil']==z)&((a['ndel']==d) if d>=0 else True)
            counts[(z,d)]=int(sel.sum())
            for j,(ys,xs) in enumerate(xy):
                hists[(z,d,j)]=np.histogram2d(ys[sel],xs[sel],bins=(yedges,xedges))[0]
                records.append(dict(zfoil=z,ndel=d,model=MODELS[j],n=int(sel.sum()),
                    core=int(np.sum(sel&(a['quality']==2))),shoulder=int(np.sum(sel&(a['quality']==1))),
                    outside_frame=int(sel.sum()-hists[(z,d,j)].sum())))
    # Sum consistency: the combined figure contains each event exactly once.
    for z in foils:
        for j in range(2):
            np.testing.assert_array_equal(sum(hists[(z,d,j)] for d in slices[:-1]),hists[(z,-1,j)])
    output.mkdir(parents=True);(output/'plots').mkdir()
    norm=LogNorm(1,max(2,max(h.max() for h in hists.values())))
    footer=('Identical reserved cores + shoulders in both panels; run groups pooled by physical foil.\n'
            'Saved xtar; existing sieve projection. Same axes, bins and log color scale throughout.')
    figures=[]
    for z in foils:
        for d in slices:
            tag=f'delta_{d}' if d>=0 else 'all_delta'
            label=f'Delta slice {d}: [{edges[d]:g}, {edges[d+1]:g})%' if d>=0 else f'All delta slices: [{edges[0]:g}, {edges[-1]:g})%'
            fig,axes=plt.subplots(1,2,figsize=(9.4,8),layout='constrained')
            mesh=draw(axes,hists,z,d,yedges,xedges,norm,counts)
            fig.colorbar(mesh,ax=axes,label='Events per bin (shared log scale)',shrink=.7,pad=.025)
            fig.suptitle(f'HMS 6.667 GeV · foil {z:+g} cm\n{label}',fontsize=15)
            fig.supxlabel(footer,fontsize=9)
            filename=f'{foil_name(z)}_{tag}.png'
            fig.savefig(output/'plots'/filename,dpi=160,facecolor='white');plt.close(fig)
            figures.append(dict(foil=z,delta=d,file=filename,label=label))
        print(f'Plotted foil {z:+g} cm: five slices and all slices combined',flush=True)
    overviews(output,foils,slices,edges,hists,yedges,xedges,norm,counts,footer)
    write_tsv(output/'counts.tsv',records)
    np.savez_compressed(output/'histograms.npz',yedges=yedges,xedges=xedges,
        **{f'{foil_name(z)}_delta{d}_model{j}':h for (z,d,j),h in hists.items()})
    text='# Before/after Huber sieve plots\n\n'
    text+='Before: balanced full-basis SVD (68,535 training events). After: smooth Huber on all 396,320 supported development events. '
    text+='Every panel uses identical reserved cores plus shoulders; no fitting is repeated.\n\n'
    text+='[All foils, all delta slices combined](plots/all_foils_all_delta.png) · [Browser gallery](gallery.html)\n\n'
    text+='| Foil (cm) | Slice 0 | Slice 1 | Slice 2 | Slice 3 | Slice 4 | All slices |\n|---|---|---|---|---|---|---|\n'
    for z in foils:
        text+=f'| {z:+g} | '+' | '.join(f'[{"combined" if d<0 else str(d)}](plots/{foil_name(z)}_{"all_delta" if d<0 else "delta_"+str(d)}.png)' for d in slices)+' |\n'
    text+='\nAll-foil overview for each slice: '+', '.join(f'[{d}](plots/all_foils_delta_{d}.png)' for d in slices[:-1])+'.\n\n'
    text+='Delta boundaries (%): '+', '.join(map(str,edges))+'.\n\n'
    text+=footer+'\n\nBin sizes are 0.05 cm in Y and 0.08 cm in X. Outside-frame counts are printed on each panel and saved in `counts.tsv`. '
    text+='Coordinates follow the existing `elastic_diagnostics.py` sieve-cloud projection; this is not a fresh iterative replay. '
    text+='This before/after comparison includes both sample expansion and the loss change.\n'
    (output/'README.md').write_text(text)
    gallery='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Before/after Huber sieve plots</title>
<style>body{font:16px system-ui;max-width:1100px;margin:32px auto;padding:0 20px;color:#18232f;background:#f6f8fa}nav{position:sticky;top:0;background:#f6f8faf5;padding:14px 0}a{color:#15527b;margin-right:18px}img{max-width:100%;height:auto;background:white;border:1px solid #ddd}section{scroll-margin-top:80px}h1{font-size:28px}h2{margin-top:40px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:22px}figure{margin:0}figcaption{margin:7px 0 16px}</style>
<h1>Sieve reconstruction before and after the Huber trial</h1><p>Left: balanced SVD. Right: expanded smooth Huber. Identical reserved cores and shoulders in every pair.</p><nav>'''
    gallery+=''.join(f'<a href="#f{foil_name(z)}">Foil {z:+g} cm</a>' for z in foils)
    gallery+='<a href="plots/all_foils_all_delta.png">All foils combined view</a></nav>'
    for z in foils:
        gallery+=f'<section id="f{foil_name(z)}"><h2>Foil {z:+g} cm</h2><div class="grid">'
        for f in [f for f in figures if f['foil']==z]:
            gallery+=f'<figure><a href="plots/{f["file"]}"><img loading="lazy" src="plots/{f["file"]}" alt="{html.escape(f["label"])} before and after Huber"></a><figcaption>{html.escape(f["label"])}</figcaption></figure>'
        gallery+='</div></section>'
    gallery+='<p>Click a plot for full resolution. Same axes, bins, and log color scale throughout. Saved xtar and existing projection; no refit or iterative replay.</p></html>'
    (output/'gallery.html').write_text(gallery)
    shutil.copy2(__file__,output/'plot_huber_sieve.py')
    provenance=dict(source_study=str(study.resolve()),source_manifest=digest(study/'manifest.json'),
        residual_sha256=digest(saved),models=list(MODELS),input_hashes=inputs,
        events=len(a['entry']),core=int(sum(a['quality']==2)),shoulder=int(sum(a['quality']==1)),
        color_range=[1,float(norm.vmax)],delta_edges=edges,individual_plots=len(figures),overview_plots=len(slices),
        outputs={str(p.relative_to(output)):digest(p) for p in output.rglob('*') if p.is_file()})
    (output/'manifest.json').write_text(json.dumps(provenance,indent=2)+'\n')
    archive=output.with_suffix('.zip')
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as zipout:
        for p in output.rglob('*'):
            if p.is_file(): zipout.write(p,str(p.relative_to(output.parent)))
    print(f'Complete: {output}\nArchive: {archive}',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign',type=Path)
    parser.add_argument('--study',default='expansion_20260924')
    parser.add_argument('--name',default='huber_sieve_20260924')
    args=parser.parse_args()
    run(args.campaign,args.campaign/'06f_huber'/args.study,args.campaign/'07_diagnostics'/args.name)
