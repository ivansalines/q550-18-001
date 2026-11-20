"""Source-direction lagged recurrences in a declared fixed 10D gauge.

For source velocity g(t)=dA/dxi, the fixed-support correlation uses paired
ages a in [0,1] for every lag L in [1,2] (matching support independent of L):
   C(t,L) = int_0^1 <Qe g(t-a), Qe g(t-a-L)> da /
      sqrt(int_0^1 ||Qe g(t-a)||^2 da * int_0^1 ||Qe g(t-a-L)||^2 da)
Lag L is age separation, NOT a measured physical delay.
"""
import numpy as np

def project(e,v,weight):
    e=np.asarray(e,float);v=np.asarray(v,float);weight=np.asarray(weight,float)
    return v-e*np.sum(v*(weight*e),axis=-1,keepdims=True)/np.dot(weight*e,e)

def normalized_pairs(left,right,weights):
    """left/right: (..., N, D), output: (...) finite cosine per matched sample set."""
    weights=np.asarray(weights,float)
    if (weights<=0).any():raise ValueError('Metric diagonal must be positive definite')
    num=np.sum(np.sum(left*(weights*right),axis=-1),axis=-1)
    a=np.sum(np.sum(left*(weights*left),axis=-1),axis=-1)
    b=np.sum(np.sum(right*(weights*right),axis=-1),axis=-1)
    den=np.sqrt(a*b)
    return np.where(den>1e-22,num/np.maximum(den,1e-22),np.nan)

def fixed_support(source,xs,lags,n,metric,e,age_lo=0.,age_hi=1.):
    """Vectorize 81 time positions x 41 lags x n age nodes (midpoints)."""
    ages=age_lo+(np.arange(n)+.5)*(age_hi-age_lo)/n
    base=project(e,source(xs[:,None]-ages[None,:]),metric)
    lagged=project(e,source(xs[:,None,None]-lags[None,:,None]-ages[None,None,:]),metric)
    v=normalized_pairs(base[:,None,:,:],lagged,metric)
    return v

def instantaneous(source,xs,lags,metric,e):
    a=project(e,source(xs),metric)
    b=project(e,source(xs[:,None]-lags[None,:]),metric)
    return normalized_pairs(a[:,None,None,:],b[:,:,None,:],metric)

def variable_support(source,xs,lags,metric,e,step=.0125,horizon=3.):
    """Available paired ages a in [0, min(horizon-L, t-1-L)] (age support changes with lag)."""
    out=np.full((len(xs),len(lags)),np.nan)
    for j,L in enumerate(lags):
        ends=np.minimum(horizon-L,xs-1-L)
        for i,end in enumerate(ends):
            if end<.5:continue
            n=max(4,int(np.ceil(end/step)))
            ages=(np.arange(n)+.5)*(end/n)
            a=project(e,source(xs[i]-ages),metric)
            b=project(e,source(xs[i]-ages-L),metric)
            out[i,j]=normalized_pairs(a[None,:,:],b[None,:,:],metric)[0]
    return out

def ridge_row(lags,values,threshold=.02):
    """Return sampled and refined minimum with broadness and edge flags.

    Refined location by local quadratic interpolation is only an estimator;
    boundaries and broad minima must be treated as censored/uncertain.
    """
    lags=np.asarray(lags,float);values=np.asarray(values,float)
    good=np.isfinite(values)
    if good.sum()<3:return {'lag_grid':None,'lag_quad':None,'kernel_min':None,'span_002':None,'edge':True,'contrast':None}
    idx=int(np.nanargmin(values));low=float(values[idx]);refined=float(lags[idx]);n=len(values)
    is_edge=idx==0 or idx==n-1
    if not is_edge and np.all(np.isfinite(values[idx-1:idx+2])):
        h=lags[1]-lags[0];p,q,r=values[idx-1:idx+2]
        denom=p-2*q+r
        if denom>1e-8:
            delta=.5*(p-r)/denom
            if abs(delta)<=1:refined=float(lags[idx]+delta*h)
    support=(values<=low+threshold)&good
    span=float(np.max(lags[support])-np.min(lags[support])) if support.any() else 0.
    return {'lag_grid':float(lags[idx]),'lag_quad':refined,'kernel_min':low,
            'span_002':span,'edge':bool(is_edge),'contrast':float(np.nanmedian(values)-low)}

def event_nearest(xs,z):
    return int(np.argmin(np.abs(np.asarray(xs)-z)))
