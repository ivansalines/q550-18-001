#!/usr/bin/env python3
"""Pure geometry for signed, conditional Q550 two-section frame diagnostics.

All angular lifts use the original tangent/curvature Gram gauge. They are
operational frame phases, not gauge-free material torsion or proof of flip.
"""
from __future__ import annotations
import numpy as np


def polar_so(M):
    U, s, Vt = np.linalg.svd(M)
    Q = U @ Vt
    # Lie-group proxy is closest SO(n). Record correction; reflection in the
    # unrestricted O(n) polar factor must not be hidden.
    flipped = float(np.linalg.det(Q)) < 0
    if flipped:
        D = np.eye(len(s))
        D[-1, -1] = -1
        Q = U @ D @ Vt
    return Q, s, flipped


def polar_o2(S):
    U, sv, Vt = np.linalg.svd(S)
    R = U @ Vt
    det = float(np.linalg.det(R))
    angle = float(np.arctan2(R[1, 0], R[0, 0])) if det > 0 else np.nan
    return R, sv, det, angle


def frame_comparison(BA, BB, M, gap=0.20):
    """Cross-section descriptors from BA, BB and exact rebased M.

    BA,BB: 11x2 tangent/curvature frames in historical Gram gauge.
    M:     11x11 tangent cocycle B -> A, recomposed without assuming rigidity.
    """
    from history_bridge_v16 import gram_frame, descriptor
    BA = np.asarray(BA, float); BB = np.asarray(BB, float)
    M = np.asarray(M, float)
    Q, sM, corrected = polar_so(M)
    pulled_full = gram_frame(*((M @ BB).T))
    pulled_so = gram_frame(*((Q @ BB).T))
    Cfull = BA.T @ pulled_full
    Cso = BA.T @ pulled_so
    _, sv_full, parity_full, angle_full = polar_o2(Cfull)
    _, sv_so, parity_so, angle_so = polar_o2(Cso)
    pA = BA @ BA.T
    pB = BB @ BB.T
    pF = pulled_full @ pulled_full.T
    pQ = pulled_so @ pulled_so.T
    distance = float(np.linalg.norm(pF - pQ, 'fro') / np.sqrt(2.))
    dlocal = descriptor(BA, BB, M)
    return dict(
        phase_full_rad=angle_full if min(sv_full) >= gap else np.nan,
        phase_so_rad=angle_so if min(sv_so) >= gap else np.nan,
        parity_full=parity_full, parity_so=parity_so,
        min_overlap_sv_full=float(sv_full[-1]),
        min_overlap_sv_so=float(sv_so[-1]),
        max_principal_angle_full_deg=float(np.degrees(np.arccos(np.clip(sv_full[-1], 0., 1.)))),
        max_principal_angle_so_deg=float(np.degrees(np.arccos(np.clip(sv_so[-1], 0., 1.)))),
        pulled_full_vs_so_grassmann_distance=distance,
        cocycle_smin=float(sM[-1]), cocycle_smax=float(sM[0]),
        cocycle_condition=float(sM[0] / max(sM[-1], 1e-300)),
        closest_SO_det_correction=bool(corrected),
        dlocal_norm=float(np.linalg.norm(dlocal)),
    ), dlocal, (pA, pB, pF, pQ)


def continuous_lift(phase, xi, parity, overlap, min_sv=.20, max_single_step_deg=150.):
    """Conservative lift; reset at invalid parity/overlap and ambiguous jumps.

    A jump exceeding 150 deg on one step is marked as alias-risk and starts a
    new segment rather than claiming a uniquely identified winding.
    """
    p = np.asarray(phase, float); x = np.asarray(xi,float)
    pr = np.asarray(parity,float); ov = np.asarray(overlap,float)
    lift = np.full(len(p),np.nan); step = np.full(len(p),np.nan)
    seg = np.full(len(p),-1,int); alias = np.zeros(len(p),bool)
    sid = -1; last = None
    for i in range(len(p)):
        if not (np.isfinite(p[i]) and pr[i]>0 and ov[i]>=min_sv):
            last = None
            continue
        if last is None:
            sid+=1; seg[i]=sid; lift[i]=p[i]; last=i
            continue
        raw=float(p[i]-p[last]); inc=float(np.arctan2(np.sin(raw),np.cos(raw)))
        if np.degrees(abs(inc))>max_single_step_deg or x[i]<=x[last]:
            alias[i]=True;sid+=1;seg[i]=sid;lift[i]=p[i];last=i
            continue
        lift[i]=lift[last]+inc;step[i]=np.degrees(inc)/(x[i]-x[last]);seg[i]=sid;last=i
    return lift,step,seg,alias


def mark_phase_events(xi, lifted, segments, one_step_deg=45., window_xi=2., total_deg=150.):
    """Predeclared candidate gates; no threshold fitted to Q550 result."""
    xi=np.asarray(xi,float); p=np.asarray(lifted,float); s=np.asarray(segments,int)
    out=[]
    for i in range(1,len(xi)):
        if s[i]<0 or s[i]!=s[i-1]:continue
        dv=np.degrees(p[i]-p[i-1])
        if abs(dv)>=one_step_deg:
            out.append((i,'large_single_step',float(dv),float(xi[i]-xi[i-1])))
        starts=np.where((s[:i]==s[i]) & ((xi[i]-xi[:i])<=window_xi))[0]
        if len(starts)==0:continue
        j=int(starts[0]); total=np.degrees(p[i]-p[j]); width=xi[i]-xi[j]
        if width>0 and abs(total)>=total_deg:
            out.append((i,'short_window_turn',float(total),float(width)))
    return out


def shapley_three_factors(minus, plus):
    """Exact symmetric telescoping allocation of descriptor delta in 11D.

    minus/plus: tuple(BA[11x2],BB[11x2],M[11x11]). Evaluates 8
    counterfactual combinations. The hybrid pairs are NOT full Q550 solutions.
    """
    from history_bridge_v16 import descriptor
    import itertools, math
    vals={}
    for state in itertools.product((0,1),repeat=3):
        aa=tuple(plus[j] if state[j] else minus[j] for j in range(3))
        vals[state]=np.asarray(descriptor(*aa),float)
    target=vals[(1,1,1)]-vals[(0,0,0)]
    contrib=[]
    for j in range(3):
        others=[k for k in range(3) if k!=j]
        delta=np.zeros(11)
        for ss in itertools.product((0,1),repeat=2):
            bit=[0,0,0]
            for k,s in zip(others,ss):bit[k]=s
            n=sum(ss);w=math.factorial(n)*math.factorial(2-n)/math.factorial(3)
            low=tuple(bit)
            bit[j]=1
            high=tuple(bit)
            delta+=w*(vals[high]-vals[low])
        contrib.append(delta)
    error=float(np.linalg.norm(sum(contrib)-target))
    return np.asarray(contrib),target,error
