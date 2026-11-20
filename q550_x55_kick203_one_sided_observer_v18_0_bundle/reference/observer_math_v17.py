#!/usr/bin/env python3
"""Sign-invariant moving-plane observer for Q550 v17.

All geometry lives in the 11D dlocal OBSERVATION space. The support-state
impulse is transferred here ONLY through the historical two-section ODE
observation map (reference/history_bridge_v16.py).
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import expm, null_space


def unit(x, tol=1e-14):
    x=np.asarray(x,float)
    n=np.linalg.norm(x)
    if not np.isfinite(n) or n<tol:raise ValueError('degenerate vector')
    return x/n


def orthonormal_basis(A):
    q,_=np.linalg.qr(np.asarray(A,float))
    return q[:,:A.shape[1]]


def projector(B):
    B=orthonormal_basis(B)
    return B@B.T


def minimal_plane_rotor(A,B,threshold=1e-7):
    """Canonical minimal ambient SO(n) map from span(A) to span(B).

    Unique when all principal angles are < pi/2. No arbitrary O(2)
    sign gauge enters this calculation. Returns R, principal angles.
    """
    A=orthonormal_basis(A);B=orthonormal_basis(B)
    if A.shape!=B.shape:raise ValueError('incompatible planes')
    U,s,Vh=np.linalg.svd(A.T@B)
    if s.min()<threshold: raise ValueError('principal plane nearly orthogonal; minimal lift ambiguous')
    s=np.clip(s,0,1);Ar=A@U;Br=B@Vh.T
    t=[];K=np.zeros((len(A),len(A)))
    for i in range(len(s)):
        a=Ar[:,i];b=Br[:,i]
        c=float(np.clip(a@b,-1.,1.));v=b-c*a
        sn=float(np.linalg.norm(v));theta=float(np.arctan2(sn,c));t.append(theta)
        if sn<1e-13:continue
        q=v/sn
        K+=theta*(np.outer(q,a)-np.outer(a,q))
    t=np.asarray(t)
    R=expm(K)
    return R,t


def weighted_plane(X,window=9,weight='hann'):
    """Causal trailing SIGN-INVARIANT second-moment observer.

    M = sum w x x^T: exact invariance under any sample-wise x->-x.
    This is deliberately NOT the historical centered v8 covariance: v8
    reference postures are compared, not silently identified with M.
    """
    X=np.asarray(X,float);n,d=X.shape
    if window<5 or window%2!=1:raise ValueError('window must be odd >=5')
    if weight=='hann':w=np.hanning(window)
    elif weight=='rect':w=np.ones(window)
    else:raise ValueError(weight)
    w=w/w.sum()
    P=np.full((n,d,d),np.nan);E=np.full((n,d,2),np.nan);L=np.full((n,3),np.nan)
    for k in range(window-1,n):
        sub=X[k+1-window:k+1]
        scale=np.max(np.linalg.norm(sub,axis=1))
        if scale<1e-30:continue
        Z=sub/scale
        C=np.einsum('n,ni,nj->ij',w,Z,Z);C=(C+C.T)/2
        eig,V=np.linalg.eigh(C);idx=np.argsort(eig)[::-1];eig=eig[idx];V=V[:,idx]
        E[k]=V[:,:2];P[k]=V[:,:2]@V[:,:2].T
        L[k]=eig[:3]
    return P,E,L


def acute_line_angle(a,b):
    return np.arccos(np.clip(abs(unit(a)@unit(b)),0,1))


def phase_in_transport(E,P,valid):
    """Intrinsic PCA major-axis phase in each OWN minimally-transported plane.

    Modulo pi, unwrapped. Frames for perturbed and baseline trajectories are
    generated independently; this is the explicit moving-lens compensation.
    """
    phase=np.full(len(P),np.nan);speed=np.full(len(P),np.nan)
    ii=np.where(valid)[0]
    if len(ii)<2:return phase,speed
    i0=int(ii[0]);G=np.column_stack([E[i0],null_space(E[i0].T)])
    phi=[];idx=[]
    prev=i0
    for k in ii:
        if k>i0:
            # Invalid intermediate frames terminate the current segment;
            # interpolation must never invent a transport across missing data.
            if k!=prev+1:break
            R,_=minimal_plane_rotor(E[prev],E[k]);G=R@G
        e=E[k,:,0]
        phi.append(np.arctan2(G[:,1]@e,G[:,0]@e));idx.append(k);prev=k
    phi=np.unwrap(2*np.array(phi))/2
    phase[np.asarray(idx)]=phi
    if len(idx)>1:
        speed[np.asarray(idx[1:])]=np.diff(phi)
    return phase,speed


def causal_exp_smoothing(raw,alpha=.4):
    """Presentation-only phase damper. NEVER feed back into physical projector."""
    x=np.asarray(raw,float);y=np.full_like(x,np.nan);last=None
    for i,v in enumerate(x):
        if not np.isfinite(v):last=None;continue
        if last is None:last=v
        else:last=last+alpha*(v-last)
        y[i]=last
    return y


def guarded_stats(base,plus,minus,eps,window=9,alpha=.4):
    """Full time-series comparison: raw, intrinsic, and lens-compensated.

    No reorientation of dlocal source is assumed. Input sign changes are
    detected and removed by plane projectors or acute principal angles.
    """
    X=[np.asarray(a,float) for a in [base,plus,minus]]
    if X[0].shape!=X[1].shape or X[1].shape!=X[2].shape:raise ValueError('shape mismatch')
    P=[];E=[];L=[]
    for x in X:
        p,e,l=weighted_plane(x,window)
        P.append(p);E.append(e);L.append(l)
    n=len(base)
    g=np.array([np.maximum(l[:,1]-l[:,2],0)/np.maximum(l[:,1],1e-300) for l in L])
    major_gap=np.array([np.maximum(l[:,0]-l[:,1],0)/np.maximum(l[:,0],1e-300) for l in L])
    valid=np.all(g>.03,axis=0)&np.all(major_gap>.02,axis=0)
    valid[:window-1]=False
    angles=np.full((2,n,2),np.nan)
    rot_errors=np.full((2,n),np.nan)
    proj_rate=np.full(n,np.nan)
    for k in np.where(valid)[0]:
        for z in range(2):
            try:
                R,ang=minimal_plane_rotor(E[0][k],E[z+1][k]);angles[z,k]=np.rad2deg(ang)
                rot_errors[z,k]=np.linalg.norm(R@P[0][k]@R.T-P[z+1][k])
            except ValueError:
                valid[k]=False;break
        if valid[k]:proj_rate[k]=np.linalg.norm((P[1][k]-P[2][k])/(2*eps),'fro')
    theta=[];filtered=[]
    for z in range(3):
        p,_=phase_in_transport(E[z],P[z],valid)
        theta.append(p);filtered.append(causal_exp_smoothing(p,alpha))
    # Both + and - phases are measured in their respective own transported lens;
    # their reference initial basis is identical only before the kick.
    pderiv=(theta[1]-theta[2])/(2*eps)
    fdamp=(filtered[1]-filtered[2])/(2*eps)
    return dict(P=P,E=E,L=L,valid=valid,subspace_gap=g,major_gap=major_gap,
                theta=theta,filtered=filtered,phase_sensitivity=pderiv,
                damped_phase_sensitivity=fdamp,principal_angles=angles,
                transport_residual=rot_errors,projector_sensitivity=proj_rate)
