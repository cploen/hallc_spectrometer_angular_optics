"""Complete physical-foil and setting diagnostics; preview never claims event membership."""
import json
import os
import shutil
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Circle
from core_sample import read_tsv, write_tsv, digest
from spectrometer_config import from_campaign


def table(rows, fields):
    def fmt(v):
        if isinstance(v,float): return f'{v:.4g}'
        return str(v)
    return '| '+' | '.join(fields)+' |\n|'+'|'.join('---' for _ in fields)+'|\n'+''.join('| '+' | '.join(fmt(r.get(k,'')) for k in fields)+' |\n' for r in rows)


def diagnostics(out,campaign,inp,result,data):
    plots=out/'plots'; plots.mkdir()
    pages=out/'pages'; pages.mkdir()
    preview=data is None
    selected_title='Planned training events (quota)' if preview else 'Selected training events'
    stages=['Available development core events',selected_title,'Planned surplus core events' if preview else 'Surplus core events']
    fields=['available','quota','surplus']
    spec=from_campaign(campaign)
    blocked=inp['blocked']
    leaves=result['leaves']
    aggregates={}
    settings={}
    for f in inp['required']:
        for nd,(lo,hi) in enumerate(zip(inp['edges'],inp['edges'][1:])):
            aggregates[(f,nd)]=[r for r in leaves if r['zfoil']==f and r['ndel']==nd]
    for name,m in inp['meta'].items():
        for nf,f in enumerate(m['foils']):
            for nd in range(len(inp['edges'])-1):
                settings[(name,nf,nd)]=[r for r in leaves if r['rungroup']==name and r['foil']==nf and r['ndel']==nd]
    def counts(rows,field):
        values=defaultdict(int)
        for r in rows:
            values[(r['xscol'],r['yscol'])]+=r[field]
        return values
    vmax=max([2]+[max(counts(rows,'available').values(),default=0) for rows in aggregates.values()])
    norm=LogNorm(1,vmax)
    cmap=plt.get_cmap('viridis')
    def circle(ax,rows,field,title,fraction=False):
        available=counts(rows,'available'); values=counts(rows,field)
        for ix in range(spec.nx):
            for iy in range(spec.ny):
                pair=ix,iy; x,y=spec.ys(iy),spec.xs(ix)
                if pair in blocked:
                    ax.text(x,y,'×',color='#a11b1b',ha='center',va='center',fontsize=20)
                elif pair not in values:
                    ax.plot(x,y,'.',color='#adb4bb',ms=4)
                else:
                    v=values[pair]
                    color_value=v/available[pair] if fraction and available[pair] else v
                    rgb=cmap(Normalize(0,1)(color_value) if fraction else norm(v))[:3] if v else (1,1,1)
                    ax.add_patch(Circle((x,y),.64,facecolor=rgb,edgecolor='#8b9299',lw=.6))
                    label=f'{color_value:.2f}' if fraction and available[pair] else str(v)
                    ax.text(x,y,label,ha='center',va='center',fontsize=8.5,color='white' if sum(a*b for a,b in zip(rgb,[.2126,.7152,.0722]))<.48 else '#111')
        pos=[v for k,v in values.items() if v>0 and k not in blocked]
        avpos=sum(v>0 for k,v in available.items() if k not in blocked)
        total=sum(values.values()); av=sum(available.values())
        fraction_text=f'{total/av:.1%}' if av else 'n/a'
        span=f'{min(pos):,}–{max(pos):,}' if pos else 'none'
        ax.set_title(title,fontsize=12,pad=10)
        ax.set(xlim=(-7.4,7.4),ylim=(-12,12),xlabel='Y sieve (cm)',ylabel='X sieve (cm)')
        ax.set_aspect('equal')
        ax.text(.5,-.13,f'N = {total:,} | fraction = {fraction_text}\nPositive holes {len(pos)}/{avpos} | range {span}',transform=ax.transAxes,ha='center',va='top',fontsize=10)
    footer='× blocked  ·  hollow 0: observed label, zero events  ·  gray dot: absent label (illumination unknown)'
    yedges=np.linspace(-7,7,281); xedges=np.linspace(-13,13,326)
    hist={}
    def cloud_data(rows, identity):
        if data is None: return None
        pieces=[[],[],[]]
        # Empty panels still retain each configured setting/foil/delta identity.
        for name,nf,nd in identity:
            a=data[name]
            base=(a['foil']==nf)&(a['ndel']==nd)&(a['quality']==2)&(a['sample']!=2)&(a['balance_excluded']==0)
            for j,keep in enumerate((base,base&(a['sample']==1),base&(a['sample']==3))):
                pieces[j].append(np.column_stack((a['ysieve'][keep],a['xsieve'][keep])))
        ans=[]
        for p in pieces:
            xy=np.concatenate(p) if p else np.empty((0,2))
            h=np.histogram2d(xy[:,0],xy[:,1],bins=[yedges,xedges])[0]
            ans.append((h,len(xy)))
        return ans
    for key,rows in aggregates.items():
        f,nd=key
        ids=[(name,nf,nd) for name,m in inp['meta'].items() for nf,z in enumerate(m['foils']) if z==f]
        hist[('a',key)]=cloud_data(rows,ids)
    for key,rows in settings.items(): hist[('s',key)]=cloud_data(rows,[key])
    cloudmax=max([2]+[int(h.max()) for hs in hist.values() if hs is not None for h,n in hs])
    cloudnorm=LogNorm(1,cloudmax)
    def cloud(ax,hs,j,title):
        if hs is None:
            ax.text(.5,.5,'Event clouds unavailable\nFrozen event masks missing\nCounts preview only',transform=ax.transAxes,ha='center',va='center',fontsize=12)
            subtitle='No event membership inferred'
        else:
            h,n=hs[j]
            ax.pcolormesh(yedges,xedges,np.ma.masked_equal(h.T,0),cmap=cmap,norm=cloudnorm,rasterized=True)
            subtitle=f'N = {n:,}; outside frame = {n-int(h.sum()):,}'
        ax.set_title(title+'\n'+subtitle,fontsize=11,pad=10)
        ax.set(xlim=(-7,7),ylim=(-13,13),xlabel='Y sieve (cm)',ylabel='X sieve (cm)')
        ax.set_aspect('equal')
    def save(fig,stem):
        fig.savefig(plots/f'{stem}.png',dpi=160,facecolor='white'); plt.close(fig)
    links={}
    for kind,collection in [('a',aggregates),('s',settings)]:
        for key,rows in collection.items():
            if kind=='a':
                f,nd=key; stem=f'foil_{f:g}_delta{nd}'
                label=f'Physical foil {f:+g} cm · combined settings'
            else:
                name,nf,nd=key; f=inp['meta'][name]['foils'][nf]
                stem=f'{name}_foil{nf}_delta{nd}'; label=f'{name} · local foil {nf} · {f:+g} cm'
            lo,hi=inp['edges'][nd:nd+2]
            title=f'HMS core balance · {label}\nDelta {nd}: [{lo:g}, {hi:g})% · width {hi-lo:g} percentage points'
            fig,axes=plt.subplots(1,3,figsize=(17,9))
            fig.subplots_adjust(left=.05,right=.91,bottom=.23,top=.81,wspace=.26)
            for ax,field,t in zip(axes,fields,stages): circle(ax,rows,field,t)
            fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.94,.3,.013,.43]),label='Events per hole · shared campaign log scale')
            fig.suptitle(title,fontsize=15)
            fig.text(.5,.04,footer+'\n'+('Same-label counts summed across settings; see setting views for focal-plane coverage.' if kind=='a' else 'Distinct setting leaves; protected holdout excluded.'),ha='center',fontsize=10)
            save(fig,stem+'_comparison')
            fig,ax=plt.subplots(figsize=(7,9)); fig.subplots_adjust(bottom=.2,top=.82,right=.83)
            circle(ax,rows,'quota',selected_title)
            fig.suptitle(f'{f:+g} cm · delta {nd}\n'+('Combined settings' if kind=='a' else name),fontsize=14)
            fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.86,.28,.022,.43]),label='Events per hole')
            save(fig,stem+'_training')
            fig,ax=plt.subplots(figsize=(7,9)); fig.subplots_adjust(bottom=.2,top=.84,right=.83)
            circle(ax,rows,'quota','Per-hole training / available',True)
            fig.suptitle(f'{f:+g} cm · delta {nd} · '+('combined settings' if kind=='a' else name),fontsize=12)
            fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=cmap),cax=fig.add_axes([.86,.28,.022,.43]),label='Training / available')
            save(fig,stem+'_fraction')
            fig,axes=plt.subplots(1,3,figsize=(17,9)); fig.subplots_adjust(left=.05,right=.91,bottom=.13,top=.8,wspace=.26)
            for j,ax in enumerate(axes): cloud(ax,hist[(kind,key)],j,stages[j])
            if not preview: fig.colorbar(plt.cm.ScalarMappable(norm=cloudnorm,cmap=cmap),cax=fig.add_axes([.94,.24,.013,.43]),label='Events per bin · shared campaign log scale')
            fig.suptitle(title,fontsize=15)
            fig.text(.5,.035,'Fixed axes and bins: 0.05 cm Y × 0.08 cm X. Protected holdout excluded.\n'+('Combined settings may broaden the cloud; inspect individual settings before interpreting.' if kind=='a' else 'Frozen accepted quality only; no density boundaries.'),ha='center',fontsize=10)
            save(fig,stem+'_clouds')
            fig,ax=plt.subplots(figsize=(7,9)); fig.subplots_adjust(bottom=.13,top=.8,right=.83)
            cloud(ax,hist[(kind,key)],1,selected_title)
            fig.suptitle(f'{f:+g} cm · delta {nd}\n'+('Combined settings' if kind=='a' else name),fontsize=13)
            if not preview:
                fig.colorbar(plt.cm.ScalarMappable(norm=cloudnorm,cmap=cmap),cax=fig.add_axes([.86,.24,.02,.43]),label='Events per bin')
            save(fig,stem+'_training_cloud')
            page=pages/f'{stem}.md'; links[(kind,key)]=page.name
            text=f'# {label} · delta {nd}\n\n'+('**Counts preview: training and surplus are quotas, not selected membership. Event clouds are unavailable.**\n\n' if preview else '')
            text+=f'Interval [{lo:g}, {hi:g})%, width {hi-lo:g} percentage points.\n\n'
            text+=f'![Available / training / surplus](../plots/{stem}_comparison.png)\n\n![Training fraction](../plots/{stem}_fraction.png)\n\n![Event clouds](../plots/{stem}_clouds.png)\n\n'
            text+=f'[Full-resolution training map](../plots/{stem}_training.png) · [Full-resolution training cloud](../plots/{stem}_training_cloud.png)\n\n'
            if kind=='a':
                text+='Same-label counts are summed across settings; this does not imply identical focal-plane coverage.\n\n'
                for name,m in inp['meta'].items():
                    for nf,z in enumerate(m['foils']):
                        if z==f: text+=f'- [{name}, local foil {nf}]({name}_foil{nf}_delta{nd}.md)\n'
            text+='\n'+table(rows,['rungroup','foil','xscol','yscol','available','q_hole','hole_cap','capacity','quota','actual_selected','surplus','qa_low','blocked'])
            page.write_text(text)
    # Complete foil×delta matrices, including zeros.
    for field,filename in [('available','available'),('quota','training_quotas'),('actual_selected','actual_selected'),('fraction','selected_available_fraction')]:
        matrix=[]
        for f in inp['required']:
            row={'zfoil':f}
            for nd in range(len(inp['edges'])-1):
                d=next(r for r in result['deltas'] if r['zfoil']==f and r['delta_low']==inp['edges'][nd])
                row[f'delta{nd}']=d['quota']/d['available'] if field=='fraction' and d['available'] else ('' if field=='fraction' else d[field])
            matrix.append(row)
        write_tsv(out/'tsv'/f'{filename}_matrix.tsv',matrix)
    setting_rows=[]
    for (name,nf,nd),rows in settings.items():
        row=dict(rungroup=name,foil=nf,zfoil=inp['meta'][name]['foils'][nf],ndel=nd,delta_low=inp['edges'][nd],delta_high=inp['edges'][nd+1])
        for k in ('available','capacity','quota','surplus','excluded_core','total','holdout','holdout_core','noncore_development','geometry_excluded_development'):
            row[k]=sum(r[k] for r in rows)
        row['actual_selected']='' if preview else row['quota']
        row['q_hole']=rows[0]['q_hole'] if rows else 0
        row['hole_cap']=rows[0]['hole_cap'] if rows else 0
        row['vanished']=row['available']>0 and row['quota']==0
        setting_rows.append(row)
    write_tsv(out/'tsv/settings.tsv',setting_rows)
    write_tsv(out/'tsv/qa_low.tsv',[r for r in leaves if r['qa_low']],list(leaves[0]) if leaves else None)
    historical=campaign/'06b_svd_fit/diagnostics/balance'
    history=[]
    if (historical/'counts.tsv').exists():
        saved=read_tsv(historical/'counts.tsv'); prov=json.loads((historical/'provenance.json').read_text())
        if sum(int(r['training']) for r in saved)!=prov['selected']: raise ValueError('Historical totals disagree with verified provenance')
        (out/'historical').mkdir()
        for filename in ('counts.tsv','provenance.json','optics.dat'):
            shutil.copy2(historical/filename,out/'historical'/filename)
        for d in result['deltas']:
            nd=inp['edges'].index(d['delta_low'])
            total=sum(int(r['training']) for r in saved if inp['meta'][r['rungroup']]['foils'][int(r['foil'])]==d['zfoil'] and int(r['ndel'])==nd)
            history.append(dict(zfoil=d['zfoil'],ndel=nd,historical_GMM_SVD=total,new_core_quota=d['quota'],new_actual_selected=d['actual_selected']))
        write_tsv(out/'tsv/historical_comparison.tsv',history)
        write_tsv(out/'tsv/historical_foil_totals.tsv',[dict(zfoil=f,historical_GMM_SVD=sum(r['historical_GMM_SVD'] for r in history if r['zfoil']==f),new_core_quota=result['foil_budget']) for f in inp['required']])
    delta_links=[]
    for nd,(lo,hi) in enumerate(zip(inp['edges'],inp['edges'][1:])):
        for typ in ('circles','clouds'):
            fig,axes=plt.subplots(1,len(inp['required']),figsize=(25,9))
            fig.subplots_adjust(left=.035,right=.94,bottom=.22,top=.81,wspace=.34)
            for ax,f in zip(np.atleast_1d(axes),inp['required']):
                if typ=='circles': circle(ax,aggregates[(f,nd)],'quota',f'{f:+g} cm')
                else: cloud(ax,hist[('a',(f,nd))],1,f'{f:+g} cm')
            if typ=='circles' or not preview:
                fig.colorbar(plt.cm.ScalarMappable(norm=norm if typ=='circles' else cloudnorm,cmap=cmap),cax=fig.add_axes([.962,.27,.01,.43]),label='Events per hole' if typ=='circles' else 'Events per bin')
            fig.suptitle(f'Delta {nd}: [{lo:g}, {hi:g})% · width {hi-lo:g} pp · '+selected_title+'\nPhysical foils in target order · same-label aggregation across settings',fontsize=18)
            fig.text(.5,.055,footer if typ=='circles' else 'Fixed frame and bins across all foils. Combined settings may broaden appearance; inspect setting views.',ha='center',fontsize=12)
            save(fig,f'delta{nd}_{typ}')
        ds=[r for r in result['deltas'] if r['delta_low']==lo]
        text=f'# Delta {nd}: [{lo:g}, {hi:g})%\n\nWidth {hi-lo:g} percentage points. Counts per configured slice, without width weighting.\n\n'
        text+=('**Counts preview only. Actual membership and event clouds require frozen event files.**\n\n' if preview else '')
        for d in ds:d['selected_fraction']=d['quota']/d['available'] if d['available'] else ''
        text+=table(ds,['zfoil','available','q_delta','delta_cap','hole_capacity','capacity','quota','actual_selected','selected_fraction','surplus'])
        text+=f'\n![Across-foil training maps](../plots/delta{nd}_circles.png)\n\n![Across-foil event clouds](../plots/delta{nd}_clouds.png)\n\n'
        for f in inp['required']: text+=f'- [{f:+g} cm: full comparison and settings]({links[("a",(f,nd))]})\n'
        (pages/f'delta{nd}.md').write_text(text); delta_links.append(f'- [Delta {nd}: [{lo:g}, {hi:g})%, width {hi-lo:g} pp](pages/delta{nd}.md)\n')
    index='# HMS 6.667 core balance\n\n'
    index+=('**Counts-only preview. All training counts are planned quotas. Actual selected counts are blank; no event IDs or real allocation clouds have been produced.**\n\n' if preview else '**Frozen core reallocation; no SVD fit run.**\n\n')
    index+=f"Fixed hole/delta multipliers {inp['cfg']['balance']['hole']:g}/{inp['cfg']['balance']['delta']:g}; qa_min={inp['cfg']['qa_min']} is reporting only.\n\nExact campaign budget: **{result['training_total']:,} = {len(inp['required'])} × {result['foil_budget']:,}**.\n\n"
    if result['zero_capacity_foils']: index+=f"**STOP: zero feasible capacity at {result['zero_capacity_foils']}; no export or solver use.**\n\n"
    index+=table(result['foils'],['zfoil','available','capacity','quota','actual_selected','hole_loss','delta_loss','budget_loss','surplus','limiting'])+'\n'
    index+=''.join(delta_links)+'\n'
    index+='[Exact missing inputs](INPUT_AVAILABILITY.md) · [All leaves](tsv/leaves.tsv) · [Setting totals](tsv/settings.tsv) · [QA: populated holes below 10](tsv/qa_low.tsv) · [Delta capacities](tsv/deltas.tsv)\n\n'
    index+='P25 references use original positive eligible counts. Geometry exclusions precede references. No automatic relaxation, reserve, weights, or legacy 80%/400 caps. Surplus includes all unselected eligible core events.\n\n'
    index+='Equal foil totals do not enforce matching delta cells or equal setting contributions. Count equality does not establish equal leverage, conditioning, or matrix stability.\n\n'
    if history:
        index+='Historical comparison uses different populations: historical GMM → SVD versus frozen accepted development cores. Historical membership was reconstructed and verified against the saved log; original solver ID lists were not saved.\n\n'
        index+='[Actual historical foil/delta totals alongside new quotas](tsv/historical_comparison.tsv)\n\n'
        for fn in ('HMS_6p667_SVD_SIEVE_BALANCE.md','HMS_6p667_RG01_EVENT_CLOUDS.md'):
            index+=f'- [Preserved historical context: {fn}]({os.path.relpath(historical/fn,out)})\n'
    (out/'HMS_6p667_CORE_BALANCE.md').write_text(index)
    (out/'diagnostic_provenance.json').write_text(json.dumps(dict(circle_scale=[1,vmax],cloud_scale=[1,cloudmax],
         frame=[-7,7,-13,13],bins=[280,325],aggregate_views=len(aggregates),setting_views=len(settings),
         mode=result['mode'],qa_min=inp['cfg']['qa_min']),indent=2)+'\n')
