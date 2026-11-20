"""Quadrature, 11D crossing, and objective scalar rate decomposition.

This module is independent of the original Q550 solver; it can be tested with
synthetic linear/nonlinear systems. Terms named 'material' are measured in
full-state Euclidean 11D coordinates, not asserted physical rigid-body angles.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq


def cumulative_simpson(xi, g, slope):
    """Simpson on full cells and half-cell quadratic interpolation at midpoints.

    The grid must be uniform. The value at every other node is exactly the
    v18.6 Simpson accumulation, and the intervening value is an exact integral
    of the same local quadratic interpolant through its left, mid, right nodes.
    """
    x=np.asarray(xi,float); a=np.asarray(g,float)
    if x.ndim!=1 or len(x)<3 or (len(x)-1)%2: raise ValueError('odd number of uniform nodes required')
    if a.shape!=(len(x),11):raise ValueError('expected 11D source integrand at each node')
    dx=np.diff(x);h=float(dx[0])*2.
    if not np.allclose(dx,dx[0],rtol=1e-9,atol=1e-12): raise ValueError('non-uniform grid')
    if slope<=0:raise ValueError('support scale must be positive')
    out=np.zeros_like(a)
    for k in range(0,len(x)-2,2):
        left,mid,right=a[k:k+3]
        out[k+1]=out[k]+slope*h/24.*(5.*left+8.*mid-right)
        out[k+2]=out[k]+slope*h/6.*(left+4.*mid+right)
    return out


def cumulative_coarse(xi,g,slope):
    """Every-fourth-node Simpson re-estimate using 0.05 xi vs 0.025 xi."""
    x=np.asarray(xi,float);a=np.asarray(g,float)
    return x[::4],cumulative_simpson(x[::2],a[::2],slope)[::2]


def continuous_integral(xi, integrand, integral, slope):
    return CubicHermiteSpline(xi,integral,slope*np.asarray(integrand,float),axis=0)


def material_vectors(P, early, total):
    return np.asarray(P,float)@np.asarray(early,float),np.asarray(P,float)@np.asarray(total,float)


def sign_and_metric_rates(P,J,slope,early,total,g):
    """Unnormalized dot rate in xi; exact algebraic decomposition.

    N = early^T G total, G=P^T P, and total'=slope*g.
    N' = early^T G' total + early^T G total', where
    G' = slope*P^T*(J^T+J)*P.
    """
    P=np.asarray(P,float);J=np.asarray(J,float)
    E=np.asarray(early,float);T=np.asarray(total,float);g=np.asarray(g,float)
    G=P.T@P
    n=float(E@(G@T))
    metr=float(slope * E@(P.T@(J.T+J)@P)@T)
    birt=float(slope * E@(G@g))
    return n,metr,birt


def kinematic_rate_check(xi, total_fn, P_fn, J_fn, g_fn, slope, early, delta=1e-4):
    def dot(s):
        P=P_fn(s);early_v=P@early;tot=P@total_fn(s)
        return float(early_v@tot)
    d=(dot(xi+delta)-dot(xi-delta))/(2*delta)
    n,m,b=sign_and_metric_rates(P_fn(xi),J_fn(xi),slope,early,total_fn(xi),g_fn(xi))
    return float(d),float(m+b)


def roots_from_grid(xi, f, start=1., end=8., zero_tol=1e-13):
    x=np.asarray(xi,float); x=x[(x>=start-1e-12)&(x<=end+1e-12)]
    values=np.asarray([f(v) for v in x])
    if not np.all(np.isfinite(values)):raise ValueError('nonfinite crossing function')
    r=[]
    for i in range(len(x)-1):
        a,b=float(x[i]),float(x[i+1]); fa,fb=float(values[i]),float(values[i+1])
        if abs(fa)<=zero_tol:r.append(a)
        if fa*fb<0:r.append(float(brentq(f,a,b,xtol=1e-12,rtol=1e-13)))
    if abs(float(values[-1]))<=zero_tol:r.append(float(x[-1]))
    clean=[]
    for val in r:
        if not clean or abs(val-clean[-1])>1e-8:clean.append(val)
    return clean,values


def polar_rotation(P):
    U,_,Vt=np.linalg.svd(np.asarray(P,float),full_matrices=False)
    return U@Vt


def angle_deg(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    na=float(np.linalg.norm(a));nb=float(np.linalg.norm(b))
    if na<=1e-25 or nb<=1e-25:return float('nan')
    return float(np.rad2deg(np.arccos(np.clip((a@b)/na/nb,-1.,1.))))
