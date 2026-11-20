#!/usr/bin/env python3
"""Exact finite secant decomposition of the two-section descriptor.

All factors A[11,2], B[11,2], M[11,11] are physical-history outputs.
For X+=(A+,B+,M+) and X-=(A-,B-,M-), let C=(X++X-)/2.
Then d(X+) - d(X-) = S_A + S_B + S_M + I_ABC exactly, where
S_A=d(A+,C_B,C_M)-d(A-,C_B,C_M), etc. The interaction residual
I_ABC is NOT presumed zero, and conditional secants are not a
full 11x11 material-observer Jacobian.
"""
from __future__ import annotations
import numpy as np


def factor_secants(minus, plus, descriptor):
    lower = [np.asarray(v, float) for v in minus]
    upper = [np.asarray(v, float) for v in plus]
    if len(lower) != 3 or any(a.shape != b.shape for a,b in zip(lower,upper)):
        raise ValueError('incompatible BA, BB, M factor shapes')
    midpoint = [(a+b)/2 for a,b in zip(lower,upper)]
    contributions = []
    for i in range(3):
        args_l = list(midpoint); args_u = list(midpoint)
        args_l[i] = lower[i]; args_u[i] = upper[i]
        contributions.append(descriptor(*args_u)-descriptor(*args_l))
    secants = np.asarray(contributions)
    total = descriptor(*upper)-descriptor(*lower)
    interaction = total - np.sum(secants,axis=0)
    return secants, interaction, total


def match_grid(source_xi, wanted_xi, tol=1e-11):
    a=np.asarray(source_xi,float); b=np.asarray(wanted_xi,float)
    if np.any(np.diff(a)<=0) or np.any(np.diff(b)<=0):
        raise ValueError('nonmonotone sample grid')
    indices=np.searchsorted(a,b)
    if np.any(indices>=len(a)) or np.max(np.abs(a[indices]-b))>tol:
        raise ValueError('nonmatching xi grid')
    return indices


def unit(v):
    v=np.asarray(v,float); n=np.linalg.norm(v)
    if n<1e-30:raise ValueError('degenerate frozen reference axis')
    return v/n


def angle_deg(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    na=np.linalg.norm(a); nb=np.linalg.norm(b)
    if min(na,nb)<1e-30:return float('nan')
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b)/(na*nb),-1,1))))
