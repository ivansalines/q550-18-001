#!/usr/bin/env python3
"""Algebra of finite source-history windows in the declared initial Euclidean 10D frame.

A(x)=s*int_0^x g(u) du, e=A(1). For fixed h, Ar=A(x)-A(x-h)
(NEW recent history) and Ao=A(x-h) (OLD history). Both are in the
same *already co-moving* initial frame; the underlying material trajectory
and original source g are NOT changed by this diagnostic ablation.
Ar'=s[g(x)-g(x-h)] and Ao'=s*g(x-h), so Ar'+Ao'=A'=s*g(x).
Projected lines, speeds and angles require a stated SPD metric; the rank of
[e,Ar,Ao] and zero/nonzero wedges are preserved by invertible coordinate maps.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator
from reorientation_math import coherence

HORIZONS=(0.25,0.50,1.,2.)

def orth(e,arr):
    e=e/np.linalg.norm(e)
    return np.asarray(arr)-np.sum(np.asarray(arr)*e,axis=-1)[...,None]*e

def state(e,a,adot,eps=1e-10):
    w=orth(e,a);d=orth(e,adot)
    rho=float(np.linalg.norm(w))
    if rho<=eps: raise ValueError('Degenerate plane: projected radius too small')
    v=w/rho
    radial=float(v@d)
    vdot=(d-v*(v@d))/rho
    speed=float(np.linalg.norm(vdot))
    return {'w':w,'v':v,'rho':rho,'radial':radial,'vdot':vdot,'speed':speed,
            'projector':np.outer(v,v),'pdot':np.outer(v,vdot)+np.outer(vdot,v)}

def eval_at(x,A,g,s,e,z,h=None,ai=None,gi=None):
    if ai is None:ai=PchipInterpolator(x,A,axis=0)
    if gi is None:gi=PchipInterpolator(x,g,axis=0)
    if z<x[0] or z>x[-1]:raise ValueError('Outside original xi support')
    aa=ai(z);dd=s*gi(z)
    if h is None:
        return state(e,aa,dd)
    if z<h:raise ValueError('Window predates zero')
    aa=aa-ai(z-h);dd=dd-s*gi(z-h)
    return state(e,aa,dd)

def event_comparison(x,A,g,s,e,z,h,halfwidths=(.05,.25,.50),ai=None,gi=None):
    if ai is None:ai=PchipInterpolator(x,A,axis=0)
    if gi is None:gi=PchipInterpolator(x,g,axis=0)
    full=eval_at(x,A,g,s,e,z,ai=ai,gi=gi)
    recent=eval_at(x,A,g,s,e,z,h,ai=ai,gi=gi)
    old=ai(z-h)
    oldw=orth(e,old)
    old_exclusive=oldw-recent['v']*np.dot(recent['v'],oldw)
    signed=float(full['v']@recent['v'])
    angle=np.degrees(np.arccos(np.clip(abs(signed),-1.,1.)))
    # For a line, P=vv^T: no sign ambiguity. All associated magnitudes depend on metric.
    arr=np.stack([e/np.linalg.norm(e),
                  recent['w']/np.linalg.norm(recent['w']),
                  oldw/np.linalg.norm(oldw)],axis=1)
    sigma=np.linalg.svd(arr,compute_uv=False)
    record={'xi':float(z),'window_xi':h,'direction_line_angle_deg':float(angle),
        'orientation_signed_cosine':signed,
        'old_exclusive_over_full_radius':float(np.linalg.norm(old_exclusive)/full['rho']),
        'old_radius_over_full_radius':float(np.linalg.norm(oldw)/full['rho']),
        'recent_radius_over_full_radius':float(recent['rho']/full['rho']),
        'full_radius':full['rho'],'recent_radius':recent['rho'],
        'full_turn_speed':full['speed'],'recent_turn_speed':recent['speed'],
        'full_radial_rate':full['radial'],'recent_radial_rate_at_full_event':recent['radial'],
        'full_vs_recent_projector_distance':float(np.linalg.norm(full['projector']-recent['projector'])),
        'full_vs_recent_projector_derivative_distance':float(np.linalg.norm(full['pdot']-recent['pdot'])),
        'relative_rank3_smallest_singular':float(sigma[-1]),
        'rank3':int(np.sum(sigma>1e-9))}
    for hh in halfwidths:
        if z-hh-h<x[0]:
            record[f'coherence_recent_pm{hh}']=None
            record[f'coherence_full_pm{hh}']=None
            continue
        for which,horizon in (('recent',h),('full',None)):
            p=eval_at(x,A,g,s,e,z-hh,horizon,ai,gi)
            q=eval_at(x,A,g,s,e,z+hh,horizon,ai,gi)
            try:
                co=coherence(p['v'],p['vdot'],q['v'],q['vdot'])
            except ValueError:
                co=None
            record[f'coherence_{which}_pm{hh}']=co
    return record
