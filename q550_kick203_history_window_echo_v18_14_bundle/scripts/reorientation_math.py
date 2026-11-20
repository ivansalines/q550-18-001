#!/usr/bin/env python3
"""Q550 v18.13 diagnostics for residual line/sphere tangent, with declared Euclidean initial gauge.

The source-frame plane is span(e,A(xi)). Define u=e/|e|, w=(I-uu^T)A,
 rho=|w|, v=w/rho, g_perp=(I-uu^T)g, and slope ds/dxi=s.
 Then rho'=s v.g_perp and v'=s (I-uu^T-vv^T)g/rho.
A radial turning rho'=0 with v' nonzero is NOT a mechanical flip or a hysteresis certificate.
Directional coherence uses Levi-Civita parallel transport between sphere points v_minus and v_plus;
this comparison depends on the source-frame Euclidean metric (not GL(10)-invariant).
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator,CubicSpline
from scipy.optimize import brentq

def orthoproject(e,arr):
    u=e/np.linalg.norm(e)
    return arr-np.einsum('...i,i->...',arr,u)[...,None]*u

def quantities(e,a,g,s):
    a=np.asarray(a,dtype=float);g=np.asarray(g,dtype=float)
    w=orthoproject(e,a)
    d=orthoproject(e,g)
    rho=np.linalg.norm(w,axis=-1)
    if np.min(rho)<=1e-12: raise ValueError('Degenerate residual line')
    v=w/rho[...,None]
    rad=s*np.sum(v*d,axis=-1)
    t=s*(d-np.sum(v*d,axis=-1)[...,None]*v)
    dv=t/rho[...,None]
    omega=np.linalg.norm(dv,axis=-1)
    return {'rho':rho,'v':v,'radial_rate':rad,'source_transverse':t,
            'vdot':dv,'speed':omega,
            'fraction_outside':np.linalg.norm(t,axis=-1)/np.maximum(np.abs(s)*np.linalg.norm(d,axis=-1),1e-30)}

def roots(x,y,mode='linear'):
    x=np.asarray(x);y=np.asarray(y)
    pairs=[i for i in range(len(x)-1) if y[i]*y[i+1]<0.]
    if mode=='linear':
        return [float(x[i]-(x[i+1]-x[i])*y[i]/(y[i+1]-y[i])) for i in pairs]
    interp={'pchip':PchipInterpolator,'cubic':CubicSpline}[mode](x,y)
    return [float(brentq(interp,float(x[i]),float(x[i+1]),xtol=1e-13)) for i in pairs]

def unit(x):
    a=np.asarray(x,float); n=np.linalg.norm(a)
    if n<1e-14:raise ValueError('Zero tangent vector')
    return a/n

def sphere_tangent_parallel_transport(a,b,t):
    """Shortest geodesic transport on unit sphere, from a to b; requires non-antipodal points."""
    a=unit(a);b=unit(b);t=t-(t@a)*a
    c=float(a@b)
    if c<-0.999999:raise ValueError('Near-antipodal sphere positions')
    y=t-(t@b)/(1.+c)*(a+b)
    return y-(y@b)*b

def coherence(vminus,vdminus,vplus,dvplus):
    t0=unit(vdminus);t1=unit(dvplus)
    moved=unit(sphere_tangent_parallel_transport(vminus,vplus,t0))
    return float(np.clip(moved@t1,-1.,1.))

def event_at(x,e,A,g,s,z):
    # Interpolate raw vectors first, then compute normalized geometry; avoids interpolating normalized vectors.
    ap=PchipInterpolator(x,A,axis=0)(z)
    gp=PchipInterpolator(x,g,axis=0)(z)
    q=quantities(e,ap[None,:],gp[None,:],s)
    return {k:(v[0] if isinstance(v,np.ndarray) else v) for k,v in q.items()}

def root_pair_coherence(x,e,A,g,s,root,halfwidth):
    p=event_at(x,e,A,g,s,root-halfwidth)
    q=event_at(x,e,A,g,s,root+halfwidth)
    cos=coherence(p['v'],p['vdot'],q['v'],q['vdot'])
    endpoint=float(np.arccos(np.clip(p['v']@q['v'],-1.,1.)))
    return {'cosine':cos,'change_angle_rad':float(np.arccos(cos)),
            'endpoint_direction_angle_rad':endpoint,
            'left_speed':float(p['speed']),'right_speed':float(q['speed']),
            'left_radial_rate':float(p['radial_rate']),'right_radial_rate':float(q['radial_rate'])}

def span_three_relative_volume(e,A,g):
    """Gauge-specific normalized 3-volume (its *vanishing* is GL invariant)."""
    T=np.stack([e/np.linalg.norm(e),A/np.linalg.norm(A),g/np.linalg.norm(g)],axis=1)
    return float(np.prod(np.linalg.svd(T,compute_uv=False)))
