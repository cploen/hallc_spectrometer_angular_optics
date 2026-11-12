"""Readable scientific figures for the additive correction experiment."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

EDGES=np.array([-15,-12,-10,-8,-5,0,5,8,10,12,14,15.])
TARGETS=('xptar','yptar','ytar','xtar')
UNITS=('mrad','mrad','mm','mm')
LABELS=dict(oldfit='Original 2018',refit='Completed 2018 refit',EN='Refit + EN correction',
            Beam='Refit + beam correction',final='Accepted result')


def report(out,a,pred,truth,records,en,beam,baseline_score,model_info,verdicts,minimum):
    nd=np.searchsorted(EDGES,a['delta'],side='right')-1
    foils=sorted(np.unique(a['ztarT']))
    populations=dict(core=(a['quality']==2)&~a['excluded'],
                     shoulder=(a['quality']==1)&~a['excluded'])
    stats={}
    for population,pop in populations.items():
        for model,p in pred.items():
            for j,target in enumerate(TARGETS):
                sd=np.full((len(foils),len(EDGES)-1),np.nan)
                bias=sd.copy();count=np.zeros(sd.shape,int)
                for i,z in enumerate(foils):
                    for d in range(len(EDGES)-1):
                        m=pop&(a['ztarT']==z)&(nd==d)
                        count[i,d]=m.sum()
                        if m.any():
                            r=p[m,j]-truth[m,j];sd[i,d]=r.std();bias[i,d]=r.mean()
                stats[population,model,target]=(sd,bias,count)

    pairs=[('refit','oldfit'),('EN','refit'),('Beam','refit'),('Beam','EN'),('final','oldfit')]
    for population in populations:
        for metric in ('sd_change','bias_change','absolute_sd'):
            fig,axes=plt.subplots(len(pairs),4,figsize=(32,19),layout='constrained')
            for row,(model,reference) in enumerate(pairs):
                for j,target in enumerate(TARGETS):
                    ax=axes[row,j]
                    sd,bias,n=stats[population,model,target]
                    bs,bb,_=stats[population,reference,target]
                    if metric=='sd_change':
                        value=100*(sd/np.maximum(bs,1e-12)-1);title='Standard deviation change (%)'
                        im=ax.imshow(np.ma.masked_invalid(value),cmap='RdBu_r',vmin=-20,vmax=20,aspect='auto')
                    elif metric=='bias_change':
                        value=bias-bb;title=f'Signed bias change ({UNITS[j]})'
                        finite=np.abs(value[np.isfinite(value)])
                        vmax=max(float(finite.max()) if finite.size else 0.,1e-5)
                        im=ax.imshow(np.ma.masked_invalid(value),cmap='RdBu_r',vmin=-vmax,vmax=vmax,aspect='auto')
                    else:
                        value=sd;title=f'Standard deviation ({UNITS[j]})'
                        im=ax.imshow(np.ma.masked_invalid(value),cmap='viridis',aspect='auto')
                    ax.set_title(f'{target}: {LABELS[model]}\n'+(f'minus {LABELS[reference]}\n' if metric!='absolute_sd' else '')+title,fontsize=10,pad=10)
                    ax.set(xticks=range(len(EDGES)-1),yticks=range(len(foils)),yticklabels=[f'{z:g}' for z in foils],
                           xticklabels=[f'[{l:g},{h:g})' for l,h in zip(EDGES[:-1],EDGES[1:])],ylabel='Foil (cm)',xlabel='Saved delta (%)')
                    ax.tick_params(axis='x',labelrotation=40,labelsize=8)
                    for i,d in np.ndindex(value.shape):
                        if n[i,d]:
                            color='white' if metric=='sd_change' and abs(value[i,d])>10 else 'black'
                            sparse='*' if n[i,d]<minimum else ''
                            text=(f'{value[i,d]:+.1f}%' if metric=='sd_change' else f'{value[i,d]:.3g}')
                            ax.text(d,i,text+f'\nN={n[i,d]}{sparse}',ha='center',va='center',fontsize=6.4,color=color)
                    fig.colorbar(im,ax=ax,shrink=.7,pad=.015)
            caption='Blue = narrower; red = wider. Colors saturate at ±20%; text is uncapped.' if metric=='sd_change' else 'Bias is shown separately from residual width; no fitted alignment shifts.'
            fig.suptitle(f'{population.title()} — {caption}\nIdentical reserved events; common iterative reconstruction; stored constants retained, zero external offsets. * N<{minimum}: uncertified sparse cell.',fontsize=14)
            fig.savefig(out/'plots'/f'{population}_{metric}.png',dpi=160);plt.close(fig)

    # Four panels parallel to the user's delta-refit summary, with explicit
    # validation status and constraints instead of an unqualified best-fit star.
    fig,ax=plt.subplots(2,2,figsize=(15,10),layout='constrained')
    for name,color in [('EN','C1'),('Beam','C0')]:
        rr=[r for r in records if r['name']==name and r['valid']]
        k=sorted(set(r['terms'] for r in rr))
        best=[min((r for r in rr if r['terms']==n),key=lambda r:r['score']) for n in k]
        ax[0,0].plot(k,[np.sqrt(r['score']) for r in best],'.-',label=name,color=color)
        ax[0,1].plot(k,[100*(1-r['score']/baseline_score) for r in best],'.-',label=name,color=color)
        ax[1,0].plot(k,[r['condition'] for r in best],'.-',label=name,color=color)
    for name,spec in [('EN',en),('Beam',beam)]:
        ax[0,0].plot(spec['terms'],np.sqrt(spec['score']),'*',ms=14,label=f'Selected {name}')
    ax[0,0].axhline(np.sqrt(baseline_score),color='.5',ls='--',label='Unchanged')
    ax[0,0].set(title='A. Grouped validation residual width',ylabel='Within-cell xtar standard deviation (mm)')
    ax[0,1].set(title='B. Reduction in within-cell variance',ylabel='Variance improvement (%)')
    ax[1,0].set(title='C. Conditioning of scaled correction design',ylabel='Maximum fold condition number',yscale='log')
    for axis in (ax[0,0],ax[0,1],ax[1,0]):
        axis.set_xlabel('Nonzero correction terms; constant fixed at zero');axis.grid(alpha=.2);axis.legend(fontsize=8)
    pop=populations['core']|populations['shoulder']
    # Pooled within-cell width, using event counts rather than cell weights.
    values=[]
    for model in ('oldfit','refit','EN','Beam'):
        ss=n=0
        for z in foils:
            for d in range(len(EDGES)-1):
                m=pop&(a['ztarT']==z)&(nd==d)
                if m.any():ss+=m.sum()*np.var(pred[model][m,3]-truth[m,3]);n+=m.sum()
        values.append(np.sqrt(ss/n))
    bars=ax[1,1].bar(range(4),values,color=['.6','.4','C1','C0'])
    ax[1,1].set(xticks=range(4),xticklabels=['Original 2018','2018 refit','EN correction','Beam correction'],
                 ylabel='Within-cell xtar standard deviation (mm)',title='D. Reserved development sample')
    for b,v in zip(bars,values):ax[1,1].text(b.get_x()+b.get_width()/2,v,f'{v:.4f}',ha='center',va='bottom',fontsize=9)
    status='; '.join(f'{name}: {"PASS" if v["passed"] else "REJECTED"} downstream guards' for name,v in verdicts.items())
    fig.suptitle('Conservative xtar correction — HMS 2018\n'+status,fontsize=16)
    fig.supxlabel('Run-group validation; ≤6 correction terms, beam width 30, regularization retained.\n2% downstream width limit overall and per populated foil/delta cell; mean biases checked separately. Reserved sample previously inspected.',fontsize=10)
    fig.savefig(out/'plots/xtar_refit_summary.png',dpi=170);plt.close(fig)
    (out/'REPORT.md').write_text('# Conservative xtar correction\n\n'+status+'\n\n'
        'Read xtar_refit_summary.png, core_sd_change.png and shoulder_sd_change.png. '
        'Signed bias changes are separate pages. Red width changes mean worse. '
        'Final transport matrix is unchanged; the additive correction is a separate JSON model. '
        'The manifest records acceptance or fallback to zero correction. '
        'These are fixed-vertex offline development results, not matched HCANA replays.\n')
