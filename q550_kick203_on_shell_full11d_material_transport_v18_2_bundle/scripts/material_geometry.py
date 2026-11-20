#!/usr/bin/env python3
"""Pure geometry for v18.2 on-flow material transport.

Two different spaces must NOT be conflated: support-state 11D perturbations
and 11D two-section dlocal observer differences. This module acts ONLY on the
physical support-state coordinates and their explicitly chosen Euclidean metric.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import expm


def unit(v, floor=1e-13):
    v = np.asarray(v, float).reshape(-1)
    n = float(np.linalg.norm(v))
    if not np.isfinite(n) or n < floor:
        raise ValueError(f'degenerate vector: norm={n}')
    return v/n


def frame_from_flow(y, rhs, jac, curvature_floor=1e-11):
    """Tangent-curvature Gram frame from the unperturbed support ODE, no h stencil.

    With f=dy/ds and t=f/|f|, dt/ds=(I-tt') J f / |f|.
    This represents one observational gauge, NOT a rigid physical body frame.
    """
    y = np.asarray(y, float).reshape(11)
    f = np.asarray(rhs(y), float).reshape(11)
    J = np.asarray(jac(y), float).reshape(11,11)
    t = unit(f)
    dt = (J@f - t*float(t@(J@f)))/np.linalg.norm(f)
    curvature = float(np.linalg.norm(dt))
    if not np.isfinite(curvature) or curvature < curvature_floor:
        raise ValueError(f'curvature frame degenerate: {curvature}')
    E = np.column_stack((t, dt/curvature))
    if np.linalg.norm(E.T@E-np.eye(2))>1e-8:
        raise RuntimeError('nonorthogonal material baseline frame')
    return E, curvature


def nearest_plane_rotation(E0, E1, singular_floor=1e-7):
    """Minimal ambient rotation between unoriented subspaces of rank 1 or 2.

    This is a geometrical subtraction convention; it is not a physical SO(11)
    symmetry of the model, nor proof of material torsion.
    """
    E0 = np.asarray(E0, float); E1 = np.asarray(E1, float)
    if E0.shape != E1.shape or E0.shape[0] != 11 or E0.ndim != 2:
        raise ValueError('expected paired 11x1 or 11x2 orthonormal frames')
    U, sigma, Vh = np.linalg.svd(E0.T @ E1)
    if sigma.min() < singular_floor:
        raise ValueError(f'plane transport becomes ambiguous, min cosine={sigma.min()}')
    A=E0@U; B=E1@Vh.T
    K=np.zeros((11,11)); angles=[]
    for j in range(len(sigma)):
        a=A[:,j]; b=B[:,j]
        cosine=np.clip(float(a@b),-1.,1.)
        q=b-cosine*a
        sine=float(np.linalg.norm(q)); angle=float(np.arctan2(sine,cosine))
        angles.append(angle)
        if sine>1e-13:
            q=q/sine
            K += angle*(np.outer(q,a)-np.outer(a,q))
    R=expm(K)
    # Check projector map, since signed frame axes are only observational gauge.
    err=float(np.linalg.norm(R@E0@E0.T@R.T-E1@E1.T,ord='fro'))
    if err>1e-7 or np.linalg.norm(R.T@R-np.eye(11),ord='fro')>1e-8:
        raise RuntimeError(f'invalid minimal plane rotation: err={err}')
    return R, np.asarray(angles)


def build_moving_transport(s_times, state_fn, rhs, jac, max_step, rank=2,
                           curvature_floor=1e-11, singular_floor=1e-7):
    """Accumulate minimum rotations on subdivided exact physical times.

    Returns orthogonal maps R(s): initial-baseline-frame -> baseline-frame at s.
    All requested sample instants must be included in s_times; no interpolation
    of R is performed. Use a second subdivision to test gauge discretization.
    """
    target = np.asarray(sorted(set(float(x) for x in s_times)), float)
    if len(target)<2 or np.min(np.diff(target))<=0:
        raise ValueError('invalid physical grid')
    if not max_step>0:raise ValueError('max_step must be positive')
    grid=[target[0]]
    for a,b in zip(target[:-1],target[1:]):
        count=int(np.ceil((b-a)/max_step))
        for j in range(1,count):grid.append(a+(b-a)*j/count)
        grid.append(b)
    frames=[]; curvature=[]
    for s in grid:
        E,k=frame_from_flow(state_fn(s),rhs,jac,curvature_floor)
        frames.append(E[:,:rank]);curvature.append(k)
    R=np.eye(11); desired={float(k):None for k in target}
    desired[float(grid[0])]=R.copy()
    max_angles=[]
    for i in range(1,len(grid)):
        Q,ang=nearest_plane_rotation(frames[i-1],frames[i],singular_floor)
        R=Q@R
        max_angles.append(float(np.max(ang)))
        if grid[i] in desired:desired[grid[i]]=R.copy()
    assert all(x is not None for x in desired.values())
    max_orth=max(np.linalg.norm(Q.T@Q-np.eye(11),ord='fro') for Q in desired.values())
    if max_orth>1e-7:raise RuntimeError('moving transport accumulated numerical nonorthogonality')
    return desired,dict(steps=len(grid)-1, rank=rank, max_principal_step_deg=float(np.degrees(max(max_angles))),
                        min_curvature=float(min(curvature)), max_orthogonality_error=float(max_orth))


def split_increment(u0,u1,R0,R1):
    """EXACT discrete kinematic decomposition, no derivative estimation.

    u1-u0 = R0*(U1-U0) + (R1-R0)*U1, U_i=R_i.T*u_i.
    First term = co-moving (intrinsic-to-chosen-gauge) increment;
    second term = frame convective contribution. Not a claim of rigid motion.
    """
    U0=R0.T@u0; U1=R1.T@u1
    internal=R0@(U1-U0); convected=(R1-R0)@U1
    closure=float(np.linalg.norm(u1-u0-internal-convected))
    return internal,convected,closure


def angle_deg(a,b,minimum_norm=1e-11):
    a=np.asarray(a,float); b=np.asarray(b,float)
    na=float(np.linalg.norm(a)); nb=float(np.linalg.norm(b))
    if min(na,nb)<=minimum_norm:return float('nan')
    return float(np.degrees(np.arccos(np.clip(float(a@b)/(na*nb),-1.,1.))))
