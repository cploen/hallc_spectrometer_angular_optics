#!/usr/bin/env python3
"""Conservative additive HMS xtar correction; all transport functions frozen.

Run with --checkout pointing at the established optics checkout. No replay
files are installed or modified. All positions are cm except matrix inputs (m).
"""
import argparse
import itertools
import json
from pathlib import Path
import shutil
import sys
import warnings
import numpy as np

EDGES = np.array([-15,-12,-10,-8,-5,0,5,8,10,12,14,15.])
TARGETS = ('xptar','yptar','ytar','xtar')
UNITS = ('mrad','mrad','mm','mm')
OLD_HASH = 'b2036a935db3f17189fcf489d1017da030a8bb58de7d5a091384ef3d3b60c92c'
REFIT_HASH = '8644a345f39c6a352b0cd0431f4f06f1ac8e976902b167ada1aa4fb097c53fd6'


def subset(a, keep):
    return {k:v[keep] for k,v in a.items()}


def features(a):
    """No truth, hole/foil/run identity, or reconstructed output as predictors."""
    u = np.column_stack((a['xfp']/100,a['xpfp'],a['yfp']/100,
                         a['ypfp'],a['reactycalc']/100))
    names = ['xfp_m','xpfp_rad','yfp_m','ypfp_rad','beam_y_m']
    cols, labels = [], []
    for degree in (1,2):
        for ids in itertools.combinations_with_replacement(range(5),degree):
            cols.append(np.prod(u[:,ids],axis=1)*np.cos(a['angle']))
            labels.append('cos(angle)*'+'*'.join(names[i] for i in ids))
    x = np.column_stack(cols)
    if not np.isfinite(x).all(): raise ValueError('Nonfinite features')
    return x, labels


class Reconstruction:
    def __init__(self,a,rows):
        from sieve_afterburner import coefficients
        self.a = a
        degree = np.minimum(np.abs(np.rad2deg(a['angle'])),50.)
        self.anchor = -a['reactycalc']-.1*(2.37-.086*degree+.0012*degree**2)
        self.z = a['rz']*np.cos(a['angle'])
        self.poly = coefficients(a,rows)
        self.base, self.base_limited = self.predict(np.zeros(len(self.z)))

    def evaluate(self,xtar):
        from sieve_afterburner import evaluate
        return evaluate(self.poly,xtar,[0.,0.,0.])

    def predict(self,correction):
        """Apply g to the vertex relation once per update, never cumulatively.

        At g=0 this is exactly sieve_afterburner.reconstruct: fixed saved
        reaction vertex, maximum five updates, 2 mrad xptar stopping criterion.
        """
        correction = np.asarray(correction)
        if correction.shape != self.anchor.shape or not np.isfinite(correction).all():
            raise ValueError('Invalid correction')
        result = self.evaluate(self.anchor)
        xtar = self.anchor.copy()
        active = np.ones(len(xtar),bool)
        for _ in range(5):
            proposal = self.anchor-result[:,0]*self.z+correction
            candidate = self.evaluate(proposal)
            change = np.abs(candidate[:,0]-result[:,0])*1000
            xtar[active] = proposal[active]
            result[active] = candidate[active]
            active &= change>2
            if not active.any(): break
        out = np.column_stack((result[:,0]*1000,result[:,2]*1000,
                               result[:,1]*1000,xtar*10))
        if not np.isfinite(out).all(): raise ValueError('Nonfinite reconstruction')
        return out,active


def truth(a):
    return np.column_stack((a['xptarT']*1000,a['yptarT']*1000,
                            a['ytarT']*10,a['xtarT']*10))


class Metrics:
    def __init__(self,a,minimum=100):
        self.a, self.minimum = a,minimum
        self.eligible = np.isin(a['quality'],[1,2]) & ~a.get('excluded',np.zeros(len(a['entry']),bool))
        self.groups = [('overall',self.eligible)]
        self.groups += [('core',(a['quality']==2)&self.eligible),('shoulder',(a['quality']==1)&self.eligible)]
        nd = np.searchsorted(EDGES,a['delta'],side='right')-1
        self.cell = np.array([f'{z:g}:{d}' for z,d in zip(a['ztarT'],nd)])
        _,self.inverse=np.unique(self.cell,return_inverse=True)
        for population,mask in list(self.groups):
            for z in np.unique(a['ztarT']):
                for d in range(len(EDGES)-1):
                    selected = mask&(a['ztarT']==z)&(nd==d)
                    if selected.any(): self.groups.append((f'{population}|{z:g}|{d}',selected))
        # Group membership is fixed across candidates. Store integer member
        # lists once so scoring does not scan all event masks repeatedly.
        self.members=[(name,np.flatnonzero(mask)) for name,mask in self.groups if mask.any()]
        self.score_members=[np.flatnonzero((self.inverse==i)&self.eligible) for i in np.unique(self.inverse)]

    def measure(self,residual):
        return {name:(len(mask),residual[mask].mean(axis=0),residual[mask].std(axis=0))
                for name,mask in self.members}

    def score(self,residual):
        # Event-count pooling of within-cell variances; no preferential weights.
        ss = n = 0
        for idx in self.score_members:
            r = residual[idx,3]
            if not len(r):continue
            ss += len(r)*np.var(r)
            n += len(r)
        return ss/n


def guards(metrics,base,candidate,limit=.02):
    """All named populated cells and populations must meet both guardrails."""
    b,c = metrics.measure(base),metrics.measure(candidate)
    worst,width_bad,bias_bad = -np.inf,[],[]
    for name,(n,bias,sd) in b.items():
        if n<metrics.minimum: continue
        cb,cs = c[name][1:]
        changes = cs[:3]/np.maximum(sd[:3],1e-12)-1
        worst = max(worst,float(changes.max()))
        for j in range(3):
            if changes[j]>limit+1e-10: width_bad.append((name,TARGETS[j],float(changes[j])))
            # Separate mean guard: no extra absolute bias beyond 2% of the
            # baseline width, in native reporting units. Not credited as gain.
            if abs(cb[j])>abs(bias[j])+limit*sd[j]+1e-10:
                bias_bad.append((name,TARGETS[j],float(abs(cb[j])-abs(bias[j]))))
    return dict(passed=not width_bad and not bias_bad,worst_width_increase=float(worst),
                width_violations=width_bad,bias_violations=bias_bad)


def scaled_fit(x,y,mask,kind,alpha,l1=.9):
    """No intercept/centering; all correction coefficients are penalized.

    Ridge uses direct augmented X SVD, never X-transpose-X inversion.
    """
    from sklearn.linear_model import ElasticNet
    from sklearn.exceptions import ConvergenceWarning
    beta = np.zeros(x.shape[1]); mask = list(mask)
    if not mask:return beta,dict(rank=0,condition=1.,singular_values=[])
    scale = np.sqrt(np.mean(x[:,mask]**2,axis=0))
    if np.any(scale<=1e-15): raise ValueError('Zero correction feature')
    z = x[:,mask]/scale
    ys = max(float(np.sqrt(np.mean(y*y))),1e-12)
    singular = np.linalg.svd(z,compute_uv=False)
    rank = int(np.sum(singular>singular[0]*1e-12))
    info = dict(rank=rank,condition=float(singular[0]/singular[-1]),
                singular_values=singular.tolist())
    if rank!=len(mask) or info['condition']>1e4:raise ValueError('Ill-conditioned basis')
    if kind=='enet':
        with warnings.catch_warnings():
            warnings.simplefilter('error',ConvergenceWarning)
            model = ElasticNet(alpha=alpha,l1_ratio=l1,fit_intercept=False,
                               tol=1e-8,max_iter=30000).fit(z,y/ys)
        coeff = model.coef_*ys
    else:
        aug = np.vstack((z,np.eye(len(mask))*np.sqrt(len(z)*alpha)))
        coeff = np.linalg.lstsq(aug,np.r_[y,np.zeros(len(mask))],rcond=1e-12)[0]
    beta[mask] = coeff/scale
    return beta,info


def attach_vertex(campaign,tag,sample,a,read_tsv):
    import uproot
    table, = (campaign/'05c_core_sample'/tag).glob('rungroups_*_inputs.tsv')
    settings = {r['rungroup']:r for r in read_tsv(table)}
    a = dict(a); a['rz'] = np.empty(len(a['entry']))
    for run in np.unique(a['rungroup']):
        idx = np.flatnonzero(a['rungroup']==run)
        p = campaign/'06c_core_ntuple'/tag/sample/'root'/f'Optics_{settings[run]["optics_id"]}_-1_fit_tree_gmm.root'
        with uproot.open(p) as f:r=f['TFit'].arrays(['entry','ztar','reactycalc'],library='np')
        order=np.argsort(r['entry'])
        if not np.array_equal(r['entry'][order],a['entry'][idx]):raise ValueError('Vertex ID mismatch')
        np.testing.assert_allclose(r['reactycalc'][order],a['reactycalc'][idx],rtol=0,atol=0)
        a['rz'][idx]=r['ztar'][order]
    if not np.isfinite(a['rz']).all():raise ValueError('Nonfinite saved vertex')
    return a


def validate_geometry(a):
    degree=np.minimum(np.abs(np.rad2deg(a['angle'])),50.)
    vertex=-a['reactycalc']-.1*(2.37-.086*degree+.0012*degree**2)
    z=a['ztarT']*np.cos(a['angle'])
    xp=(a['xsT']-vertex)/(168-z)
    xt=vertex-xp*z
    np.testing.assert_allclose(a['xptarT'],xp,rtol=1e-12,atol=1e-13)
    np.testing.assert_allclose(a['xtarT'],xt,rtol=1e-12,atol=1e-10)
    np.testing.assert_allclose(a['xtarT']+168*a['xptarT'],a['xsT'],rtol=0,atol=1e-10)


def run(args):
    sys.path.insert(0,str(args.checkout))
    from core_sample import digest,read_tsv,write_tsv
    from fit_elastic import load_sample,matrix_rows,checked
    from sklearn.model_selection import GroupKFold
    from threadpoolctl import threadpool_limits
    out=args.out
    out.mkdir(parents=True,exist_ok=False)
    (out/'code').mkdir();(out/'matrices').mkdir();(out/'plots').mkdir()
    for name in ('xtar_correction.py','xtar_correction_plots.py','test_xtar_correction.py','README.md'):
        p=Path(__file__).with_name(name)
        if p.exists():shutil.copy2(p,out/'code'/name)
    campaign=args.campaign.resolve();tag='inner20'
    old=campaign/'config/oldfit.dat'
    refit=campaign/'06e_beam_search/beam30_inner20/matrices/beam.dat'
    hashes={};checked(old,OLD_HASH,hashes);checked(refit,REFIT_HASH,hashes)
    for name,p in [('oldfit',old),('refit',refit)]:shutil.copy2(p,out/'matrices'/f'{name}.dat')
    manifest=dict(schema='conservative_xtar_correction_v1',status='running',
        campaign=str(campaign),settings=vars(args)|{'checkout':str(args.checkout),'out':str(out),'campaign':str(campaign)},
        constants='All stored constants retained for oldfit and refit; zero external additions; no correction intercept.',
        convention='Common fixed-saved-vertex iterative xtar; correction added once to each proposal. Not a new HCANA replay.',
        comparison_metric='Within foil/delta cell standard deviation; signed bias separately.',
        selection='Run-group 3-fold validation, no preferential delta weights; prior reserved events are development evidence.',
        limits='2% downstream standard deviation overall and each populated foil/delta cell, core/shoulder/union separately; mean bias guard = 2% baseline width.',
        sparse='Cells below minimum are displayed, not certified; corrections are not silently clipped.',inputs=hashes)
    def save():(out/'manifest.json').write_text(json.dumps(manifest,indent=2,default=str,allow_nan=False)+'\n')
    save()
    try:
        a,hh=load_sample(campaign,tag,'fit',xtar_coordinate='xtarT');hashes.update(hh)
        a=attach_vertex(campaign,tag,'fit',a,read_tsv);validate_geometry(a)
        x,labels=features(a);t=truth(a)
        rows=matrix_rows(refit);reco=Reconstruction(a,rows)
        # Independent parity against the established afterburner, on every row.
        from sieve_afterburner import reconstruct
        for rg in np.unique(a['rungroup']):
            idx=np.flatnonzero(a['rungroup']==rg);aa=subset(a,idx);aa['ry']=aa['reactycalc']
            _,_,active,target,xt=reconstruct(aa,rows,float(np.rad2deg(aa['angle'][0])),[0,0,0],return_target=True)
            expected=np.column_stack((target[:,0]*1000,target[:,2]*1000,target[:,1]*1000,xt*10))
            np.testing.assert_allclose(expected,reco.base[idx],rtol=1e-13,atol=1e-11)
            np.testing.assert_array_equal(active,reco.base_limited[idx])
        base=reco.base-t;y=a['xtarT']-reco.base[:,3]/10
        folds=np.empty(len(x),int)
        for f,(_,va) in enumerate(GroupKFold(3).split(x,groups=a['rungroup'])):folds[va]=f
        write_tsv(out/'folds.tsv',[dict(rungroup=r,entry=int(e),fold=int(f)) for r,e,f in zip(a['rungroup'],a['entry'],folds)])
        parts=[]
        for f in range(3):
            tr,va=folds!=f,folds==f
            parts.append((tr,va,Metrics(subset(a,va),args.minimum)))
        met=Metrics(a,args.minimum);baseline_score=met.score(base)
        zero=dict(name='unchanged',kind='ridge',mask=(),alpha=.03,l1=.9,amplitude=0.)
        records=[];cache={}
        print(f'{len(x):,} training events; {len(np.unique(a["rungroup"]))} run groups; {x.shape[1]} candidate columns; no intercept',flush=True)

        def assess(spec):
            key=(spec['kind'],tuple(spec['mask']),spec['alpha'],spec['l1'],spec['amplitude'])
            if key in cache:return cache[key]
            residual=np.empty_like(base);conditions=[];active_masks=[];infos=[]
            valid=True;failures=[];worst=-np.inf;maxg=0.;limited=0
            for tr,va,m in parts:
                try:
                    b,info=scaled_fit(x[tr],y[tr],spec['mask'],spec['kind'],spec['alpha'],spec['l1'])
                    b*=spec['amplitude'];g=x[va]@b
                    active_masks.append(np.flatnonzero(b).tolist());conditions.append(info['condition']);infos.append(info)
                    if np.count_nonzero(b)>args.max_terms or np.max(np.abs(g))>.2:valid=False
                    maxg=max(maxg,float(np.max(np.abs(g))))
                    # Subset pre-collapsed polynomial; no transport evaluation at true xtar.
                    child=object.__new__(Reconstruction)
                    child.anchor=reco.anchor[va];child.z=reco.z[va];child.poly=reco.poly[:,va]
                    pred,lim=child.predict(g);residual[va]=pred-t[va]
                    if np.any(lim&~reco.base_limited[va]):valid=False
                    limited+=int(lim.sum())
                    guard=guards(m,base[va],residual[va])
                    valid &= guard['passed'];worst=max(worst,guard['worst_width_increase'])
                    failures.extend(guard['width_violations']+guard['bias_violations'])
                except (ValueError,Warning,np.linalg.LinAlgError) as exc:
                    valid=False;failures.append(str(exc));break
            if len(conditions)==3:
                guard=guards(met,base,residual);valid &= guard['passed']
                worst=max(worst,guard['worst_width_increase']);failures.extend(guard['width_violations']+guard['bias_violations'])
                score=met.score(residual);mse=float(np.mean(residual[:,3]**2))
                # Width reduction cannot be bought with a larger xtar MSE.
                valid &= mse<=float(np.mean(base[:,3]**2))*(1+1e-10)
            else:score=baseline_score*1e6;mse=float(np.mean(base[:,3]**2))*1e6
            r=dict(spec,mask=list(spec['mask']),valid=bool(valid),score=float(score),mse=mse,
                terms=max(map(len,active_masks),default=0),condition=max(conditions,default=1e99),
                worst_width_increase=float(worst) if np.isfinite(worst) else 1e99,max_correction_cm=maxg,
                iteration_limited=limited,violations=len(failures),fold_supports=active_masks)
            cache[key]=r;records.append(r)
            return r

        with threadpool_limits(limits=args.threads):
            base_rec=assess(zero)
            en_candidates=[base_rec]
            for l1 in (.5,.9,1.):
                for alpha in (1.,.3,.1,.03,.01,.003,.001):
                    r=assess(dict(name='EN',kind='enet',mask=tuple(range(x.shape[1])),alpha=alpha,l1=l1,amplitude=1.))
                    if r['valid']:en_candidates.append(r)
            order=lambda r:(r['score'],r['terms'],r['condition'])
            best=min(en_candidates,key=order)
            # Prefer fewer terms within 0.5% of the best validation variance.
            en=min([r for r in en_candidates if r['score']<=best['score']*1.005],key=lambda r:(r['terms'],r['score']))
            write_tsv(out/'enet_cv.tsv',[{k:v for k,v in r.items() if k not in ('mask','fold_supports')} for r in records])
            print(f'EN: {en["terms"]} terms; width ratio={np.sqrt(en["score"]/baseline_score):.6f}',flush=True)
            supports={()}
            supports.update(tuple(s) for s in en['fold_supports'] if len(s)<=args.max_terms)
            # Seed only validated EN supports; arbitrary large bases are excluded.
            frontier=[()];best_beam=base_rec;history=[]
            visited=set();stale=0;anchor_score=baseline_score
            for step in range(args.max_terms+1):
                masks=set(supports) if step==0 else {tuple(sorted(set(m)^{j})) for m in frontier for j in range(x.shape[1]) if len(set(m)^{j})<=args.max_terms}
                masks-=visited;visited |= masks
                accepted=[]
                for mask in sorted(masks):
                    options=[assess(dict(name='Beam',kind='ridge',mask=mask,alpha=.03,l1=.9,amplitude=amp)) for amp in ((0.,) if not mask else (.1,.25,.5,1.))]
                    good=[r for r in options if r['valid']]
                    if good:accepted.append(min(good,key=order))
                if accepted:
                    accepted.sort(key=order);frontier=[tuple(r['mask']) for r in accepted[:30]]
                    best_beam=min([best_beam,*accepted],key=order)
                history.append(dict(step=step,evaluated=len(records),terms=best_beam['terms'],score=best_beam['score'],condition=best_beam['condition']))
                print(f'Beam {step}: {len(records)} candidates, {best_beam["terms"]} terms; width ratio={np.sqrt(best_beam["score"]/baseline_score):.6f}',flush=True)
                (out/'progress.json').write_text(json.dumps(history,indent=2))
                if not accepted or not frontier:break
                if step:
                    if best_beam['score']<anchor_score*.999:anchor_score=best_beam['score'];stale=0
                    else:stale+=1
                    if stale>=2:break
            all_beam=[r for r in records if r['name']=='Beam' and r['valid']]+[base_rec]
            best=min(all_beam,key=order)
            beam=min([r for r in all_beam if r['score']<=best['score']*1.005],key=lambda r:(r['terms'],r['score']))
            # A material width gain is required; zero is always admissible.
            if beam['score']>=baseline_score*.998001:beam=base_rec
            if en['score']>=baseline_score*.998001:en=base_rec
            manifest['selection'] = dict(en=en,beam=beam,baseline_variance=baseline_score,history=history,
                folds='GroupKFold(3), complete runs kept together; not the prior event-level folds',
                ridge='mean squared correction-target loss + 0.03*sum(scaled coefficients squared)',
                maximum_correction_cm=.2,maximum_terms=args.max_terms,minimum_cell_events=args.minimum)
            save()
            proposals={};model_info={}
            for name,spec in [('EN',en),('Beam',beam)]:
                b,info=scaled_fit(x,y,spec['mask'],spec['kind'],spec['alpha'],spec['l1']);b*=spec['amplitude']
                if np.count_nonzero(b)>args.max_terms:raise ValueError('Final fit exceeds term cap')
                model_info[name]=dict(spec=spec,coefficients_cm=b.tolist(),features=labels,intercept_cm=0.,final_fit=info)
                proposals[name]=b
                (out/f'{name}_correction.json').write_text(json.dumps(model_info[name],indent=2)+'\n')
            # Reserved data are not opened until both models and coefficients freeze.
            h,hh=load_sample(campaign,tag,'holdout',xtar_coordinate='xtarT');hashes.update(hh)
            if set(zip(a['rungroup'],a['entry']))&set(zip(h['rungroup'],h['entry'])):raise ValueError('Training/reserved ID overlap')
            h=attach_vertex(campaign,tag,'holdout',h,read_tsv);validate_geometry(h)
            hx,_=features(h);ht=truth(h);hr=Reconstruction(h,rows)
            pred={'oldfit':Reconstruction(h,matrix_rows(old)).base,'refit':hr.base}
            hm=Metrics(h,args.minimum);baseh=hr.base-ht;verdicts={}
            for name,b in proposals.items():
                g=hx@b;pv,lim=hr.predict(g);pred[name]=pv
                verdict=guards(hm,baseh,pv-ht)
                verdict['passed'] &= bool(np.max(np.abs(g))<=.2 and not np.any(lim&~hr.base_limited))
                verdict['xtar_width_change']=float(np.sqrt(hm.score(pv-ht)/hm.score(baseh))-1)
                verdict['xtar_mse_change']=float(np.mean((pv[:,3]-ht[:,3])**2)/np.mean(baseh[:,3]**2)-1)
                verdicts[name]=verdict
            accept=bool(verdicts['Beam']['passed'] and verdicts['Beam']['xtar_width_change']<-.001 and verdicts['Beam']['xtar_mse_change']<=0)
            chosen='Beam' if accept else 'refit'
            pred['final']=pred[chosen]
            final=proposals['Beam'] if accept else np.zeros(x.shape[1])
            (out/'final_correction.json').write_text(json.dumps(dict(schema='hms_xtar_additive_correction_v1',
                coefficients_cm=final.tolist(),features=labels,intercept_cm=0.,selected=chosen,
                transport_sha256=REFIT_HASH,constants='retained, zero external additions',
                application='g added once to fixed-vertex xtar proposal, not cumulatively',
                certification='Offline development experiment, no HCANA installation; sparse cells uncertified'),indent=2)+'\n')
            # Matrix remains identical; the separate correction cannot be encoded
            # as a fourth transport output (that column is delta).
            shutil.copy2(refit,out/'matrices/final_transport.dat')
            manifest.update(verdicts=verdicts,selected=chosen,training_events=len(x),reserved_events=len(hx),
                parity='Zero correction matched established afterburner on all training events',
                historical_closure='Unresolved; same fixed-vertex convention for every matrix, no new replay',
                deployment='No replay installation. final_transport.dat unchanged; final_correction.json is separate.')
            np.savez_compressed(out/'reserved_predictions.npz',truth=ht,entry=h['entry'],rungroup=h['rungroup'],
                delta=h['delta'],zfoil=h['ztarT'],quality=h['quality'],excluded=h['excluded'],**pred)
            summary=[]
            for name,p in pred.items():
                residual=p-ht
                for group,mask in hm.groups:
                    if not mask.any():continue
                    for j,target in enumerate(TARGETS):
                        r=residual[mask,j]
                        summary.append(dict(model=name,group=group,target=target,unit=UNITS[j],n=len(r),
                            bias=float(r.mean()),sd=float(r.std()),rms=float(np.sqrt(np.mean(r*r))),
                            sparse=bool(len(r)<args.minimum)))
            write_tsv(out/'reserved_metrics.tsv',summary)
            write_tsv(out/'search_candidates.tsv',[{k:v for k,v in r.items() if k not in ('mask','fold_supports')}|{'mask':','.join(map(str,r['mask']))} for r in records])
            from xtar_correction_plots import report
            report(out,h,pred,ht,records,en,beam,baseline_score,model_info,verdicts,args.minimum)
            for p,sha in hashes.items():checked(Path(p),sha,{})
            manifest.update(status='complete',inputs=hashes)
            save()
            print(json.dumps(dict(selected=chosen,verdicts=verdicts),indent=2),flush=True)
            print('COMPLETE',out,flush=True)
    except BaseException as exc:
        manifest.update(status='incomplete',error=str(exc));save();raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkout',type=Path,required=True)
    p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--minimum',type=int,default=100)
    p.add_argument('--max-terms',type=int,default=6)
    p.add_argument('--threads',type=int,default=4)
    args=p.parse_args()
    if args.minimum<20 or not 1<=args.max_terms<=6 or args.threads<1:p.error('Invalid conservative limits')
    sys.dont_write_bytecode=True
    run(args)


if __name__=='__main__':main()
