"""Exact algebraic source-kernel closure versus specified reduced recurrence baselines.

All quantities depend on the chosen metric G and projection e. No physical-memory
inference is justified from closure or a residual of an incomplete observer model.
"""
import numpy as np

def project(v,e,w):
    v=np.asarray(v,float);e=np.asarray(e,float);w=np.asarray(w,float)
    if np.any(w<=0):raise ValueError('Metric must be positive diagonal')
    d=np.dot(w*e,e)
    if d <= 1e-28:raise ValueError('Degenerate anchor')
    return v-e*np.sum(v*(w*e),axis=-1,keepdims=True)/d

def paired_decomposition(source,times,lags,w,e,n=64,age_low=0.,age_high=1.):
    """Return kernels and source-only components from same quadrature samples.

    Window s sampled at z-age and z-age-L. Spatial dot products G-weighted.
    closure uses cross-normalized dot product. The reconstructed cosine sum must
    agree with cross-normalized dot product to roundoff (algebraic identity).
    """
    times=np.asarray(times,float);lags=np.asarray(lags,float)
    ages=age_low+(np.arange(n)+.5)*(age_high-age_low)/n
    newer=project(source(times[:,None]-ages[None,:]),e,w)  # (T,N,D)
    older=project(source(times[:,None,None]-lags[None,:,None]-ages[None,None,:]),e,w)
    num=np.sum(newer[:,None,:,:]*older*w,axis=-1)
    aa=np.sqrt(np.sum(newer*newer*w,axis=-1))[:,None,:]
    bb=np.sqrt(np.sum(older*older*w,axis=-1))
    if np.any(aa<1e-10) or np.any(bb<1e-10):
        raise ValueError('Degenerate projected source norm, kernel undefined')
    pair_cos=np.clip(num/(aa*bb),-1.,1.)
    denom=np.sqrt(np.sum(aa**2,axis=-1)*np.sum(bb**2,axis=-1))
    exact=np.sum(num,axis=-1)/denom
    closure=np.sum(aa*bb*pair_cos,axis=-1)/denom
    gate=np.sum(aa*bb,axis=-1)/denom
    directional=pair_cos.mean(axis=-1)
    approx=gate*directional
    weight_cos=np.sum(aa*bb*pair_cos,axis=-1)/np.sum(aa*bb,axis=-1)
    # Correlation amplitude-direction covariance; completely source determined.
    covariance=exact-approx
    # Same lag pairs, but deliberately scrambled temporal pairing of older source
    # ages (circular roll by n//4) to diagnose whether age-registration matters.
    shifted=np.roll(older,max(1,n//4),axis=-2)
    scramble=np.sum(newer[:,None,:,:]*shifted*w,axis=(-1,-2))
    scramble/=np.sqrt(np.sum(newer*newer*w,axis=(-1,-2))[:,None] *
                       np.sum(shifted*shifted*w,axis=(-1,-2)))
    return dict(exact=exact,closure=closure,unit_direction=directional,
                gate=gate,gate_times_unit=approx,covariance=covariance,
                weighted_direction=weight_cos,scrambled=scramble)

def instantaneous(source,times,lags,w,e):
    x=project(source(times),e,w)[:,None,:]
    y=project(source(times[:,None]-lags[None,:]),e,w)
    numer=np.sum(x*y*w,axis=-1)
    denom=np.sqrt(np.sum(x*x*w,axis=-1)*np.sum(y*y*w,axis=-1))
    if np.any(denom<=1e-25):raise ValueError('Degenerate endpoint source')
    return numer/denom

def ridge_lag(lags,k):
    """Grid argmin: do not manufacture subgrid timing confidence from broad ridges."""
    return np.asarray(lags)[np.argmin(k,axis=-1)]

def train_endpoint_calibrator(endpoint,target,times,split=6.):
    """Fit two params to train xi<=split; predict all xi, leave xi>split untouched."""
    mask=np.asarray(times)<=split+1e-12
    x=endpoint[mask].ravel();y=target[mask].ravel()
    design=np.column_stack((np.ones_like(x),x))
    coefs=np.linalg.lstsq(design,y,rcond=None)[0]
    return coefs[0]+coefs[1]*endpoint, coefs
