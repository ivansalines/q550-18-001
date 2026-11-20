#!/usr/bin/env python3
"""Exterior-algebra primitives for GL(11)-covariant material plane diagnostics.

All angles/norms require an explicit Euclidean coordinate gauge; only the
vanishing of exterior products and exact exterior transport are intrinsic to
invertible linear coordinate transformations.
"""
from __future__ import annotations
import numpy as np
from itertools import combinations
PAIRS=np.array(list(combinations(range(11),2)),int)
TRIPLES=np.array(list(combinations(range(11),3)),int)


def blade2(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return a[...,PAIRS[:,0]]*b[...,PAIRS[:,1]]-a[...,PAIRS[:,1]]*b[...,PAIRS[:,0]]


def blade3(a,b,c):
    a,b,c=[np.asarray(v,float) for v in (a,b,c)]
    i,j,k=TRIPLES.T
    return (a[...,i]*(b[...,j]*c[...,k]-b[...,k]*c[...,j])
           -a[...,j]*(b[...,i]*c[...,k]-b[...,k]*c[...,i])
           +a[...,k]*(b[...,i]*c[...,j]-b[...,j]*c[...,i]))


def safe_ratio(num,den):
    return np.divide(num,den,out=np.full_like(np.asarray(num,float),np.nan),where=np.asarray(den)>1e-220)


def normalized_bivector(a,b):
    nb=np.linalg.norm(blade2(a,b),axis=-1)
    return safe_ratio(nb,np.linalg.norm(a,axis=-1)*np.linalg.norm(b,axis=-1))


def normalized_trivector(a,b,c):
    nb=np.linalg.norm(blade3(a,b,c),axis=-1)
    return safe_ratio(nb,np.linalg.norm(a,axis=-1)*np.linalg.norm(b,axis=-1)*np.linalg.norm(c,axis=-1))


def wedge_rate(E,K,J,P,g,slope):
    """B'=slope*(J E ^ K + E ^ J K + E ^ P g), with xi time."""
    JE=np.einsum('nij,nj->ni',J,E)
    JK=np.einsum('nij,nj->ni',J,K)
    Pg=np.einsum('nij,nj->ni',P,g)
    strain=slope*(blade2(JE,K)+blade2(E,JK))
    birth=slope*blade2(E,Pg)
    return strain,birth


def rank_from_vectors(columns,rel_tol=1e-8):
    """Scale each vector to unit norm for SVD; numerical, metric-dependent gate."""
    mat=np.column_stack([np.asarray(c,float).ravel() for c in columns])
    norms=np.linalg.norm(mat,axis=0)
    if np.min(norms)<=1e-180:return 0,np.nan,[]
    ss=np.linalg.svd(mat/norms,compute_uv=False)
    return int(np.sum(ss>rel_tol*ss[0])),float(ss[-1]/ss[0]),ss.tolist()


def congruent_blade2(T,a,b):
    """Induced exterior map C2(T)(a^b), built without assuming orthogonality."""
    return blade2(T@a,T@b)


def congruent_blade3(T,a,b,c):
    return blade3(T@a,T@b,T@c)


def normal_component_ratio(a,b,g):
    """Projection magnitude in the EXPLICIT chosen initial-coordinate metric."""
    A=np.column_stack((a,b)); r=np.linalg.norm(g)
    if r<=1e-180:return float('nan')
    return float(np.linalg.norm(g-A@np.linalg.lstsq(A,g,rcond=1e-12)[0])/r)
