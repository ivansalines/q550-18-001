"""Q550 v18.8: tested 11D metric-response geometry; no physical metric assumed.

xi is support offset.  The stored Jacobian is with respect to original support s;
slope = ds/dxi is read from v18.6, not guessed.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, minimize_scalar

VERSION='2026-10-17-q550-v18.8-material-metric-constraint-sensitivity'


def projectors_from_jacobians(A, rtol=1e-9):
    """Orthogonal projectors onto kernel(DC). A has shape (n,m,11)."""
    A=np.asarray(A,float)
    if A.ndim!=3 or A.shape[-1]!=11 or not np.isfinite(A).all():
        raise ValueError('constraint_jacobians must be finite (n,m,11)')
    out=[];ranks=[]; gaps=[]
    for mat in A:
        u,s,vt=np.linalg.svd(mat,full_matrices=True)
        rank=int(np.count_nonzero(s>rtol*max(float(s[0]) if len(s) else 0.,1e-14)))
        if rank>=11:raise ValueError('constraint kernel is zero dimensional')
        N=vt[rank:].T
        out.append(N@N.T); ranks.append(rank)
        gaps.append(float(s[rank-1]/s[rank]) if rank and rank<len(s) and s[rank]>0 else None)
    if len(set(ranks))!=1:raise ValueError('constraint Jacobian rank varies over nodes')
    return np.asarray(out),ranks[0],gaps


def metric_dot(E,K,W):
    if W.ndim==2:return np.einsum('ni,ij,nj->n',E,W,K)
    return np.einsum('ni,nij,nj->n',E,W,K)


def fixed_metric_dot_and_rate(E,K,J,P,g,slope,W):
    """Signed old-to-total dot and xi-rate for a time-independent SPD W."""
    ep=slope*np.einsum('nij,nj->ni',J,E)
    kp=slope*(np.einsum('nij,nj->ni',J,K)+np.einsum('nij,nj->ni',P,g))
    n=metric_dot(E,K,W)
    derivative=metric_dot(ep,K,W)+metric_dot(E,kp,W)
    return n,derivative


def source_metric_dot_and_rate(early,T,g,slope):
    """Inverse-cocycle pullback, W=P^-T P^-1; tangent strain cancels."""
    old=np.einsum('j,nj->n',early,T)
    derivative=slope*np.einsum('j,nj->n',early,g)
    return old,derivative


def roots_from_samples(xi,N,Np,start=1.,end=8.):
    """Grid-resolved sign-changing roots, cubic Hermite for each bracket.

    NO interval guarantee against even numbers of missed roots inside one cell.
    """
    x=np.asarray(xi,float);v=np.asarray(N,float);d=np.asarray(Np,float)
    use=(x>=start-1e-12)&(x<=end+1e-12); x=x[use];v=v[use];d=d[use]
    if len(x)<2 or not np.isfinite(v).all() or not np.isfinite(d).all():
        raise ValueError('invalid series')
    f=CubicHermiteSpline(x,v,d)
    indices=np.where(v[:-1]*v[1:]<0)[0]
    r=[float(brentq(f,float(x[i]),float(x[i+1]),xtol=1e-13)) for i in indices]
    return r


def fixed_w_from_reference(Pref):
    B=np.linalg.solve(Pref,np.eye(11));return B.T@B


def mixed_metric(Nlab,Dlab,Ninv,Dinv,alpha):
    return (1-alpha)*Nlab+alpha*Ninv,(1-alpha)*Dlab+alpha*Dinv


def dot_derivative_variable_metric(E,K,J,P,g,slope,W,Wprime):
    """Complete differential rate with metric motion; Wprime is dW/dxi."""
    ep=slope*np.einsum('nij,nj->ni',J,E)
    kp=slope*(np.einsum('nij,nj->ni',J,K)+np.einsum('nij,nj->ni',P,g))
    return metric_dot(ep,K,W)+metric_dot(E,kp,W)+metric_dot(E,K,Wprime)


def alpha_crossing_extinction(xi,Nlab,Ninv):
    """Find each negative Euclidean lobe's approximate critical alpha, where
    its most negative N_alpha reaches 0. Used only as sampled/interpolated diagnostic.
    """
    x=np.asarray(xi); q=np.asarray(Nlab); r=np.asarray(Ninv)
    mask=(x>=1)&(x<=8)
    x=x[mask];q=q[mask];r=r[mask]
    lobes=[];sel=q<0;i=0
    while i<len(x):
        if not sel[i]:i+=1;continue
        j=i
        while j+1<len(x) and sel[j+1]:j+=1
        # threshold alpha(x) at which Nalpha(x)=0, for Nlab(x)<0 and Ninv(x)>0
        denom=r-q
        aa=np.where(denom>1e-30,-q/denom,np.nan)
        v=aa[i:j+1]; im=i+int(np.nanargmax(v))
        l=max(i,im-1);u=min(j,im+1)
        if l<u:
            fq=CubicHermiteSpline(x,q,np.gradient(q,x,edge_order=2));fr=CubicHermiteSpline(x,r,np.gradient(r,x,edge_order=2))
            opt=minimize_scalar(lambda z:-float(-fq(z)/(fr(z)-fq(z))),bounds=(x[l],x[u]),method='bounded',options={'xatol':1e-13})
            z=float(opt.x);thr=float(-fq(z)/(fr(z)-fq(z)))
        else: z=float(x[im]);thr=float(aa[im])
        lobes.append({'xi_approx':z,'alpha_critical_approx':thr,'negative_lobe_range':[float(x[i]),float(x[j])]})
        i=j+1
    return lobes
