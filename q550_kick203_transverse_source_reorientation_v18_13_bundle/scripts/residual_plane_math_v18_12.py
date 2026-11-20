#!/usr/bin/env python3
"""Exact algebra in a DECLARED initial-frame Euclidean 10D gauge.

Let y=(z,theta) and Phi be the 10x10 fixed-phase tangent propagator.
Archived E=Phi e, K=Phi A, A'=s*g, where derivatives are in xi.
Pulling back by Phi eliminates the homogeneous transport of the *vectors*.
The residual line within e^perp is v=(I-ee^T)A / ||(I-ee^T)A||.
Then v' = s (I-ee^T-vv^T) g / ||(I-ee^T)A||.
All norms and angles depend on the Euclidean gauge at the initial frame;
rank, vanishing of wedge and tensor coordinate covariance do not.
"""
from __future__ import annotations
import numpy as np

def source_geometry(e:np.ndarray,a:np.ndarray,g:np.ndarray,s:float,cutoff=1e-11):
    e=np.asarray(e,float); a=np.asarray(a,float); g=np.asarray(g,float)
    en=e/np.linalg.norm(e)
    perp_a=a-np.einsum('ij,j->i',a,en)[:,None]*en
    perp_g=g-np.einsum('ij,j->i',g,en)[:,None]*en
    rho=np.linalg.norm(perp_a,axis=1)
    if np.any(rho<=cutoff): raise ValueError('Residual plane degenerates in requested interval')
    v=perp_a/rho[:,None]
    longitudinal=np.einsum('ij,ij->i',v,perp_g)*s
    transverse=s*(perp_g-np.einsum('ij,ij->i',v,perp_g)[:,None]*v)
    v_dot=transverse/rho[:,None]
    rate=np.linalg.norm(v_dot,axis=1)
    return dict(eunit=en, v=v, radius=rho, radius_rate=longitudinal,
                vdot=v_dot,turn_rate=rate, source_norm=np.linalg.norm(perp_g,axis=1),
                source_outside_plane_norm=np.linalg.norm(transverse,axis=1),
                tangent_source_fraction=np.linalg.norm(transverse,axis=1)/np.maximum(abs(s)*np.linalg.norm(perp_g,axis=1),1e-30))

def principal_line_angle(v,w):
    return np.arccos(np.clip(np.abs(np.einsum('ij,ij->i',v,w)),0,1))

def velocity_from_projector(v,vdot):
    return vdot[:,:,None]*v[:,None,:]+v[:,:,None]*vdot[:,None,:]

def max_tensor_transform_error(T, e,a,g):
    """Covariance of the raw antisymmetric source e wedge g, not of naive gauge angles."""
    base=np.outer(e,g)-np.outer(g,e)
    tr=np.outer(T@e,T@g)-np.outer(T@g,T@e)
    pred=T@base@T.T
    return float(np.linalg.norm(tr-pred)/max(np.linalg.norm(pred),1e-30))
