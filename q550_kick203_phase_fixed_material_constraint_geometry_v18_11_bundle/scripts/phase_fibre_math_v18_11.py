#!/usr/bin/env python3
"""Phase-fixed cocycle splitting and gauge-declared Grassmann velocity.

State y=(z[10],theta), theta'=2pi. If J_theta,*=0 then tangent bundle
V={dtheta=0} is invariant and P=[[Phi,q],[0,1]].
The time-dependent leaf theta(s)=theta0+2pi(s-s0) is NOT a mechanical
constraint. Its 10D flow z'=f_z(z,theta(s)) is nonautonomous.
"""
from __future__ import annotations
import numpy as np


def splitter(J,P):
    return J[:,:10,:10],J[:,:10,10],P[:,:10,:10],P[:,:10,10]


def phase_sector_checks(J,P,vectors):
    return {'max_jacobian_phase_row':float(np.max(np.abs(J[:,10,:]))),
        'max_cocycle_phase_row_to_z':float(np.max(np.abs(P[:,10,:10]))),
        'max_cocycle_phase_diagonal_from_one':float(np.max(np.abs(P[:,10,10]-1))),
        'max_source_phase_component':float(max(np.max(np.abs(a[:,10])) for a in vectors)),
        'max_phase_to_z_jacobian_column_norm':float(np.max(np.linalg.norm(J[:,:10,10],axis=1))),
        'max_phase_to_z_cocycle_column_norm':float(np.max(np.linalg.norm(P[:,:10,10],axis=1)))}


def phase_cocycle_rhs(J,P,slope):
    A,b,F,q=splitter(J,P)
    return slope*np.einsum('nij,njk->nik',A,F),slope*(np.einsum('nij,nj->ni',A,q)+b)


def unit_plane(U,cutoff=1e-12):
    """Return Euclidean-gauge orthonormal basis Q, pseudo-inverse, sv's."""
    u,s,vt=np.linalg.svd(U,full_matrices=False)
    if s[-1]<=max(cutoff*s[0],np.finfo(float).tiny):raise ValueError('plane rank insufficient')
    return u,vt.T@np.diag(1/s)@u.T,s


def projector_derivative(U,Up):
    """Differentiate orthogonal projector for full-rank 10x2 U."""
    Q,pinv,s=unit_plane(U)
    perp=np.eye(U.shape[0])-Q@Q.T
    mixed=perp@Up@pinv
    return mixed+mixed.T, mixed, Q, s


def tangential_plane_rates(E,K,J,P,g,slope):
    U=np.column_stack((E[:10],K[:10]))
    j=J[:10,:10];p=P[:10,:10]
    ut=slope*np.column_stack((j@E[:10],j@K[:10]))
    ub=slope*np.column_stack((np.zeros(10),p@g[:10]))
    total,_,Q,s=projector_derivative(U,ut+ub)
    transport,_,_,_=projector_derivative(U,ut)
    birth,_,_,_=projector_derivative(U,ub)
    return total,transport,birth,Q,s


def principal_step_angle(U0,U1):
    q0,_,_=unit_plane(U0);q1,_,_=unit_plane(U1)
    co=np.linalg.svd(q0.T@q1,compute_uv=False)
    return float(np.max(np.arccos(np.clip(co,-1,1))))
