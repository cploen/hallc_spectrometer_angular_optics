"""Provisional Y segmentation from an unlabeled ordered spacing prior.

API fit_y(raw [ypfp,yfp], angle=None, max_slots=9), predict_y(model,raw).
No campaign labels accepted. All output clusters are local ordered IDs.
Strong projection peaks establish a nearly affine band lattice. Extrapolated
slots count only when local density supports them; missing slots stay empty.
Independent 2D Gaussian ellipses then refine assignments, with covariance
shrinkage toward the smooth ordered family. All operations preserve ddof=1
standardization and an explicit orthogonal rotation.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths
from scipy.special import logsumexp
from .x_geometry import _raw_array, _slot_count

CONFIG=dict(bins=600,smoothing=1.5,strong_prominence=.10,weak_z=5.0,
            max_slots=9,covariance_shrinkage=.15,maximum_em_steps=20)

def _rotate(raw,mean,scale,angle):
    mean, scale = np.asarray(mean, float), np.asarray(scale, float)
    if (mean.shape != (2,) or scale.shape != (2,) or not np.isfinite(mean).all()
            or not np.isfinite(scale).all() or np.any(scale <= 0) or not np.isfinite(angle)):
        raise ValueError('Invalid Y standardization or rotation')
    t=np.deg2rad(angle);return ((np.asarray(raw)-mean)/scale)@np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])

def _seed(q, max_slots=9):
    # A common standardized bin width avoids tuning the KDE bandwidth to
    # the projection angle or to the outermost selected events.
    lo=min(-5.,float(q.min())-.05);hi=max(5.,float(q.max())+.05)
    h,e=np.histogram(q,bins=CONFIG['bins'],range=(lo,hi));x=(e[1:]+e[:-1])/2
    sm=gaussian_filter1d(h.astype(float),CONFIG['smoothing'])
    peaks,props=find_peaks(sm,prominence=CONFIG['strong_prominence']*sm.max(),height=.03*sm.max(),distance=7)
    if len(peaks)<3 or len(peaks)>max_slots:return None
    centers=x[peaks];gaps=np.diff(centers);spacing=np.median(gaps)
    ji=np.r_[0,np.cumsum(np.maximum(1,np.rint(gaps/spacing)).astype(int))]
    if ji[-1]>=max_slots:return None
    fit=np.polyfit(ji,centers,1);base,spacing=fit[1],fit[0]
    residual=np.sqrt(np.mean((centers-(base+spacing*ji))**2))/spacing
    if residual>.13:return None
    width=np.median(peak_widths(sm,peaks,rel_height=.5)[0])*(e[1]-e[0])
    if spacing/width<2.0:return None
    allpeaks,allprops=find_peaks(sm,prominence=.001*sm.max(),distance=3)
    out=[]
    # Extra slots are predictions, not detections. Each needs independent
    # local excess and a real maximum near the predicted position.
    for j in range(-2,int(ji[-1])+3):
        expected=base+spacing*j
        options=allpeaks[np.abs(x[allpeaks]-expected)<.23*spacing]
        if not len(options):continue
        pk=options[np.argmax(sm[options])];c=x[pk]
        core=np.abs(x-c)<.18*spacing
        side=(np.abs(x-c)>.30*spacing)&(np.abs(x-c)<.48*spacing)
        if core.sum()<2 or side.sum()<2:continue
        n=float(h[core].sum());bg=float(h[side].sum())*core.sum()/side.sum()
        excess=n-bg;z=excess/np.sqrt(max(1,n+bg*core.sum()/side.sum()))
        if z<CONFIG['weak_z'] or excess<25:continue
        out.append(dict(slot=j,center=float(c),z=float(z),excess=float(excess),expected=float(expected)))
    if len(out)<3 or len(out)>max_slots:return None
    # More supported bands are preferred only after a resolution and spacing
    # check. Resolution breaks ties and selects a manual-free default angle.
    strength=float(np.median(props['prominences'])/sm.max())
    score=len(out)+.10*np.log(spacing/width)+.10*strength-2*residual
    return dict(peaks=out,spacing=float(spacing),width=float(width),
                strong_centers=centers.tolist(),spacing_rms=float(residual),score=float(score))

def _logp(x,means,covs,weights):
    out=[]
    for mu,cov,w in zip(means,covs,weights):
        z=x-mu;iv=np.linalg.inv(cov)
        out.append(np.log(max(w,1e-12))-.5*(np.einsum('ni,ij,nj->n',z,iv,z)+np.linalg.slogdet(cov)[1]+2*np.log(2*np.pi)))
    return np.array(out).T

def _regularize_cov(covs):
    covs=np.asarray(covs);k=len(covs);j=np.arange(k)
    if (covs.shape != (k, 2, 2) or not np.isfinite(covs).all()
            or np.any(covs[:, 0, 0] <= 0) or np.any(covs[:, 1, 1] <= 0)):
        raise ValueError('Y covariance estimates need positive finite variances')
    values=np.c_[np.log(covs[:,0,0]),np.log(covs[:,1,1]),np.arctanh(np.clip(covs[:,0,1]/np.sqrt(covs[:,0,0]*covs[:,1,1]),-.99,.99))]
    smooth=np.array([np.polyval(np.polyfit(j,v,min(2,k-1)),j) for v in values.T]).T
    values=.85*values+.15*smooth
    out=[]
    for a,b,r in values:
        v1,v2=np.exp([a,b]);cross=np.tanh(r)*np.sqrt(v1*v2)
        out.append([[v1,cross],[cross,v2]])
    return np.asarray(out)+np.eye(2)*1e-6

def fit_y(raw, angle=None, max_slots=9):
    raw = _raw_array(raw, minimum=100)
    max_slots = _slot_count(max_slots)
    if angle is not None:
        if not np.isfinite(angle):
            raise ValueError('The Y angle override must be finite')
        angle = float(angle) % 180.
    mean=raw.mean(0);scale=raw.std(0,ddof=1)
    if not np.isfinite(scale).all() or np.any(scale <= 0):
        raise ValueError('Both Y coordinates must have a positive finite variance')
    scan=[]
    for a in np.arange(0.,180.,1.) if angle is None else [float(angle)]:
        xy=_rotate(raw,mean,scale,a);seed=_seed(xy[:,0], max_slots=max_slots)
        if seed is not None:scan.append((seed['score'],float(a),seed))
    if not scan:raise ValueError('No supported ordered Y geometry')
    score,angle,seed=max(scan,key=lambda a:a[0]);xy=_rotate(raw,mean,scale,angle)
    centers=np.array([p['center'] for p in seed['peaks']]);slots=[p['slot'] for p in seed['peaks']]
    lab=np.searchsorted((centers[:-1]+centers[1:])/2,xy[:,0]);k=len(centers)
    if min(np.bincount(lab, minlength=k)) < 3:
        raise ValueError('Too few events to estimate a supported Y band covariance')
    means=np.array([xy[lab==i].mean(0) for i in range(k)])
    covs=_regularize_cov([np.cov(xy[lab==i].T,ddof=1) for i in range(k)])
    weights=np.array([(lab==i).mean() for i in range(k)])
    for step in range(CONFIG['maximum_em_steps']):
        logs=_logp(xy,means,covs,weights);r=np.exp(logs-logsumexp(logs,axis=1)[:,None]);mass=r.sum(0)
        if not np.isfinite(mass).all() or np.any(mass <= 0):
            raise ValueError('A Y mixture component has no finite effective support')
        newmeans=r.T@xy/mass[:,None]
        # Preserve discovered order and prevent a weak band drifting into its
        # stronger neighbor; the spacing prior supplies its allowed neighborhood.
        newmeans[:,0]=np.clip(newmeans[:,0],centers-.25*seed['spacing'],centers+.25*seed['spacing'])
        newcov=[]
        for i in range(k):
            diff=xy-newmeans[i];newcov.append((diff*r[:,i,None]).T@diff/max(1,mass[i]))
        means=newmeans;covs=_regularize_cov(newcov);weights=mass/len(xy)
    return dict(version='provisional_y_spacing_v1',theta_deg=angle,mean=mean.tolist(),scale=scale.tolist(),
      means=means.tolist(),covariances=covs.tolist(),weights=weights.tolist(),slots=slots,
      band_count=k,seeds=seed,scan_summary=[dict(angle=a,score=s,bands=len(p['peaks'])) for s,a,p in scan],config=dict(CONFIG, max_slots=max_slots))

def predict_y(model, raw, details=False):
    raw = _raw_array(raw)
    if len(raw) == 0:
        if details:
            return dict(labels=np.empty(0, dtype=int), posterior=np.empty(0),
                        mahalanobis2=np.empty(0), coordinates=np.empty((0, 2)))
        return np.empty(0, dtype=int)
    xy=_rotate(raw,np.asarray(model['mean']),np.asarray(model['scale']),model['theta_deg'])
    means=np.asarray(model['means']);covs=np.asarray(model['covariances']);weights=np.asarray(model['weights'])
    if (means.ndim != 2 or means.shape[1:] != (2,) or covs.shape != (len(means), 2, 2)
            or weights.shape != (len(means),) or not len(means)
            or not np.isfinite(means).all() or not np.isfinite(covs).all()
            or not np.isfinite(weights).all() or np.any(weights <= 0)
            or not np.allclose(covs, covs.transpose(0, 2, 1))
            or np.any(np.linalg.eigvalsh(covs) <= 0)):
        raise ValueError('Invalid fitted Y Gaussian model')
    logs=_logp(xy,means,covs,weights);posterior=np.exp(logs-logsumexp(logs,axis=1)[:,None]);lab=np.argmax(posterior,axis=1)
    # Keep classifications reviewable; use a coverage flag rather than silently
    # removing events before their assignment can be inspected.
    diff=xy-means[lab];distance=np.einsum('ni,nij,nj->n',diff,np.linalg.inv(covs)[lab],diff)
    if details:return dict(labels=lab,posterior=posterior.max(1),mahalanobis2=distance,coordinates=xy)
    return lab
