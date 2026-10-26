"""Whole-distribution periodic X discovery without an artificial region split.

The smooth inverse coordinate proposes an ordered polynomial family. Occupancy
is tested separately. Input columns are [xpfp, xfp]; no labels are accepted.
"""
import numpy as np
from scipy.ndimage import maximum_filter, gaussian_filter
from scipy.optimize import minimize
from .x_geometry import _fit, _slot_count, _raw_array, global_centers, global_coordinates

CONFIG=dict(rotation_deg=135.,bins=100,fft_padding=400,starts=18,min_q_cycles=5.,max_q_cycles=24.,max_p_cycles=16.,optimization_maxiter=120,refine_rounds=4)

def fit_x(raw, ndel=None, angle=None, max_slots=9):
    """Discover the complete X family; return a JSON-serializable geometry.

    ``raw`` contains [xpfp, xfp]. The default standardized rotation is 135
    degrees; a finite manual angle replaces it. Occupancy and physical sieve
    labeling are subsequent operations. The HMS/SHMS profile supplies the
    maximum number of possible slots (9/11), without forcing their occupancy.
    """
    raw = _raw_array(raw, minimum=100)
    max_slots = _slot_count(max_slots)
    if angle is not None and not np.isfinite(angle):
        raise ValueError('The X angle override must be finite')
    rotation = CONFIG['rotation_deg'] if angle is None else float(angle) % 180.
    mean=raw.mean(0);scale=raw.std(0,ddof=1)
    if not np.isfinite(scale).all() or np.any(scale <= 0):
        raise ValueError('Both X coordinates must have a positive finite variance')
    z=(raw-mean)/scale
    theta=np.deg2rad(rotation);q=z@np.array([np.cos(theta),np.sin(theta)]);p=z@np.array([-np.sin(theta),np.cos(theta)])
    qlim=np.quantile(q,[.002,.998]);plim=np.quantile(p,[.002,.998]);qc=.5*sum(qlim);pc=.5*sum(plim);qs=.5*np.diff(qlim)[0];ps=.5*np.diff(plim)[0]
    if qs <= 0 or ps <= 0 or not np.isfinite([qs, ps]).all():
        raise ValueError('Degenerate rotated X coordinate range')
    Q=(q-qc)/qs;P=(p-pc)/ps
    H,qe,pe=np.histogram2d(Q,P,bins=CONFIG['bins'],range=[[-1,1],[-1,1]])
    # Smooth only slightly: coarse-bin noise should not supply the oscillation.
    H=gaussian_filter(H,.35);gq=.5*(qe[1:]+qe[:-1]);gp=.5*(pe[1:]+pe[:-1]);QG,PG=np.meshgrid(gq,gp,indexing='ij')
    take=H.ravel()>0;w=H.ravel()[take];w=w/w.sum();QG=QG.ravel()[take];PG=PG.ravel()[take]
    X=np.c_[QG,PG,PG**2,QG*PG,QG*PG**2,QG**2]
    F=np.abs(np.fft.fftshift(np.fft.fft2(H,s=[CONFIG['fft_padding']]*2)))/H.sum();freq=np.fft.fftshift(np.fft.fftfreq(CONFIG['fft_padding'],d=2/CONFIG['bins']))
    good=(freq[:,None]>=CONFIG['min_q_cycles']/2)&(freq[:,None]<=CONFIG['max_q_cycles']/2)&(abs(freq[None,:])<=CONFIG['max_p_cycles']/2)
    peaks=np.argwhere(good&(F==maximum_filter(F,size=5)))
    peaks=sorted(peaks,key=lambda ij:F[tuple(ij)],reverse=True)[:CONFIG['starts']]
    def objective(b):
        E=np.exp(2j*np.pi*(X@b));C=np.sum(w*E);power=abs(C)**2
        grad=2*np.real(np.conj(C)*(2j*np.pi*np.sum((w*E)[:,None]*X,axis=0)))
        # Gentle curvature penalty prevents unconstrained high-frequency folding.
        penalty=2e-5*np.sum(b[2:]**2)
        gg=np.zeros(6);gg[2:]=4e-5*b[2:]
        return -power+penalty,-grad+gg
    solutions=[]
    for a,b in peaks:
        init=np.array([freq[a],freq[b],0,0,0,0.])
        bounds=[(2.5,12.),(-10,10),(-5,5),(-5,5),(-3,3),(-3,3)]
        opt=minimize(objective,init,jac=True,method='L-BFGS-B',bounds=bounds,options={'maxiter':CONFIG['optimization_maxiter'],'ftol':1e-10})
        # Inverse q derivative must remain positive over observed distribution.
        dv=opt.x[0]+opt.x[3]*P+opt.x[4]*P**2+2*opt.x[5]*Q
        if np.quantile(dv,.001)<=.35:continue
        value=-objective(opt.x)[0];solutions.append((value,opt.x,opt.success))
    if not solutions:raise ValueError('No monotone whole-distribution periodic solution')
    solutions.sort(key=lambda r:r[0],reverse=True)
    # Select via concentration of all events in the physical number of consecutive cells, penalizing
    # harmonics that explain sharp peaks by inventing intervening empty columns.
    XP=np.c_[Q,P,P**2,Q*P,Q*P**2,Q**2];candidates=[]
    for coherence,b,success in solutions:
        v=XP@b;phase=np.angle(np.mean(np.exp(2j*np.pi*v)))/(2*np.pi);v-=phase
        offsets=range(int(np.floor(np.quantile(v,.001)))-1,int(np.ceil(np.quantile(v,.999)))+1)
        scores=[(float(np.mean((v>=off-.5)&(v<=off+max_slots-.5))),off) for off in offsets]
        fraction,offset=max(scores);res=abs(v-np.rint(v));sharp=np.mean(np.exp(-.5*(res/.15)**2))
        score=coherence*fraction**4
        candidates.append(dict(score=score,coherence=coherence,b=b,phase=phase,offset=offset,fraction=fraction,sharpness=sharp))
    best=max(candidates,key=lambda c:c['score']);v=XP@best['b']-best['phase']-best['offset'];j=np.rint(v)
    keep=(j>=0)&(j<max_slots)&(abs(v-j)<.24)&(abs(P)<1.03)
    if keep.sum()<100:raise ValueError('Insufficient periodic ridge support')
    p_scale=float(np.quantile(abs(p),.99));w0=np.ones(keep.sum());coef=_fit(j[keep],p[keep],q[keep],w0,2,2,p_scale,max_slots)
    fan=dict(format='provisional_unsplit_periodic_curved_fan_v1',mean=mean.tolist(),scale=scale.tolist(),theta_deg=float(np.rad2deg(theta)),coefficients=coef.tolist(),p_scale=p_scale,j_degree=2,p_degree=2,possible_indices=list(range(max_slots)),max_slots=max_slots,theta_was_overridden=bool(angle is not None),delta_slice=int(ndel or 0),bootstrap='Whole-distribution periodic likelihood; no regions, no seam, no local track bootstrap',discovery_configuration=CONFIG.copy())
    fan['discovery_configuration'].update(rotation_deg=rotation, max_slots=max_slots)
    for iteration in range(CONFIG['refine_rounds']):
        co=global_coordinates(fan,raw);u=co['u'];jj=np.rint(np.nan_to_num(u,nan=-99));ok=np.isfinite(u)&(abs(u-jj)<.28)&(jj>=0)&(jj<max_slots)
        if ok.sum()<100:break
        widths=np.bincount(jj[ok].astype(int),minlength=max_slots);weights=1/np.sqrt(np.maximum(widths[jj[ok].astype(int)],100))
        weights*=np.sqrt(np.median(np.maximum(widths[widths>0],100)))
        coef=_fit(jj[ok],p[ok],q[ok],weights,2,2,p_scale,max_slots);fan['coefficients']=coef.tolist()
    pp=np.linspace(float(p.min()),float(p.max()),200);curves=global_centers(fan,pp);gaps=np.diff(curves,axis=0)
    # A quadratic family can cross between fractional indices even if integer
    # centers do not; check derivative at both ends of the inverse domain.
    dleft=global_centers(fan,pp,indices=[-.499,-.5]);dright=global_centers(fan,pp,indices=[max_slots-.5,max_slots-.501]);der=min(float(np.min((dleft[0]-dleft[1])/.001)),float(np.min((dright[0]-dright[1])/.001)))
    co=global_coordinates(fan,raw)
    fan['diagnostics']=dict(p_support=np.quantile(p,[.001,.999]).tolist(),full_event_p_extent=[float(p.min()),float(p.max())],min_adjacent_band_gap=float(gaps.min()),min_fractional_order_derivative=der,monotone_on_rectangular_support=bool(der>0),full_extent_crossing_pairs=np.flatnonzero(np.any(gaps<=0,axis=1)).tolist(),periodic_coherence=float(best['coherence']),periodic_nine_cell_fraction=float(best['fraction']),finite_inversion_fraction=float(co['finite_geometry'].mean()),periodic_seed_coefficients=best['b'].tolist(),periodic_phase=float(best['phase']),periodic_offset=int(best['offset']),periodic_start_count=len(solutions),joint_support_events=int(keep.sum()),no_split_anywhere=True)
    fan['diagnostics']['periodic_cell_fraction'] = float(best['fraction'])
    if der<=0:raise ValueError('Joint curved family loses order over full observed support')
    if co['finite_geometry'].mean()<.90:raise ValueError('Joint geometry covers under 90 percent of events for %d possible bands' % max_slots)
    return fan
