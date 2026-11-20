#!/usr/bin/env python3
"""Source-age pair kernel in the DECLARED initial co-moving 10D Euclidean gauge.

Age a=z-u; the input is original cumulative source A(u), represented by PCHIP.
Finite age slabs q_i = Q_e^G[A(z-a_i)-A(z-a_(i+1))] where the integration
is truncated at physical source time u>=1, so no zero-net pre-anchor slabs.

For lag L=k*width, define a cosine of stacked paired slab vectors:
   K(z,L) = sum_i <q_i,q_(i+k)>_G /
     sqrt[sum_i ||q_i||_G^2 * sum_i ||q_(i+k)||_G^2].
The denominator uses the SAME matching pairs, and K(z,0)=1. The overlap
shrinks with lag; always record its duration. This is NOT an observer invariant,
material constitutive kernel, physical propagation time, or mechanical energy.
"""
import numpy as np

def project(e, v,weights):
    e=np.asarray(e,float);v=np.asarray(v,float);weights=np.asarray(weights,float)
    return v - np.sum(v*(weights*e),axis=-1)[...,None]*e/float(np.dot(e*weights,e))

def age_slabs(ai,e,z,width,weights,anchor=1.,horizon=None):
    length=max(0.,z-anchor)
    if horizon is not None:length=min(length,horizon)
    n=int(np.floor((length+1e-12)/width))
    if n<1:return np.empty((0,len(e))),0.
    edges=width*np.arange(n+1)
    q=project(e,ai(z-edges[:-1])-ai(z-edges[1:]),weights)
    return q,n*width

def kernel_from_slabs(q,k,weights):
    n=len(q)
    if k<0 or k>=n:return np.nan,n-k
    a=q[:n-k];b=q[k:]
    # A shift toward older source: paired ages i and i+k, not a full pair matrix.
    numerator=float(np.einsum('ij,j,ij->',a,weights,b))
    aa=float(np.einsum('ij,j,ij->',a,weights,a))
    bb=float(np.einsum('ij,j,ij->',b,weights,b))
    den=np.sqrt(aa*bb)
    return (numerator/den if den>1e-20 else np.nan),n-k

def kernel(ai,e,z,width,lags,weights,horizon=None,min_overlap=.5):
    q,covered=age_slabs(ai,e,z,width,weights,1.,horizon)
    rows=[]
    for L in lags:
        k=int(round(L/width))
        if abs(k*width-L)>1e-9:raise ValueError('lag must be integer multiple of age bin')
        corr,npairs=kernel_from_slabs(q,k,weights)
        overlap=npairs*width
        if overlap<min_overlap-1e-12:corr=np.nan
        rows.append((float(L),corr,float(overlap),int(npairs)))
    return rows,covered

def effective_lag_minimum(rows,lo=1.,hi=2.):
    valid=[r for r in rows if lo-1e-12<=r[0]<=hi+1e-12 and np.isfinite(r[1])]
    if not valid:return None,None
    winner=min(valid,key=lambda r:r[1]);return winner[0],winner[1]

def pairing_identity(q,weights):
    """Numerical identity for all pair contributions including diagonal."""
    gram=(q*weights)@q.T
    return float(abs(np.einsum('i,i,i->',q.sum(0),weights,q.sum(0)) - np.sum(gram)))
