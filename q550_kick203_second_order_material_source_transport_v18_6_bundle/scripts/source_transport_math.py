"""Pure 11D source transport algorithms for v18.6.

q'(s) = J(s)q(s) + 1/2 D^2 f(y0(s))[u(s),u(s)]
u(s) = P(s) v, P'(s)=J(s)P(s), P(s0)=I, q(s0)=0.
Thus q(s) = P(s) int_s0^s P(t)^(-1) S2(t) dt.
For the paired echo, r_even(s1+) / amp_abs**2 = q(s1) + O(amp_abs**2).
"""
import numpy as np

def source_hessian_via_jac(jac, s, y, tangent, h):
    """Compute half Hessian contraction with analytic Jacobian central differences.

    h is an ABSOLUTE perturbation of the 11D state along unit tangent.
    Returns S2 = .5*D²f(y)[tangent,tangent] (full 11D).
    """
    y=np.asarray(y,float); u=np.asarray(tangent,float)
    nu=float(np.linalg.norm(u))
    if nu == 0.:return np.zeros_like(y)
    if h<=0:raise ValueError('Hessian step must be positive')
    w=u/nu
    return (nu*nu/(4.*h)) * ((jac(s,y+h*w)-jac(s,y-h*w))@w)

def transported_integrand(P, source, condition_gate=1e10):
    c=float(np.linalg.cond(P))
    if not np.isfinite(c) or c>condition_gate:raise ValueError(f'Unsafe inverse cocycle condition {c}')
    return np.linalg.solve(P,source),c

def simpson_cell(g0,gm,g1, ds):
    return float(ds)/6.*(np.asarray(g0)+4*np.asarray(gm)+np.asarray(g1))

def signed_even(plus,minus,baseline):
    return 0.5*(np.asarray(plus)+np.asarray(minus))-np.asarray(baseline)

def relative_error(a,b, floor=1e-14):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b))/max(float(np.linalg.norm(b)),floor))

def angle_deg(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    na=np.linalg.norm(a);nb=np.linalg.norm(b)
    if min(na,nb)<1e-30:return float('nan')
    return float(np.rad2deg(np.arccos(np.clip(a@b/(na*nb),-1,1))))

def slope_finite_amplitude(K2, amps):
    """Vector extrapolation to zero amplitude, assuming K2 + K4*a^2.

    'amps' are the absolute perturbation magnitudes.
    Rescale for conditioning; this is a diagnostic, not a certificate.
    """
    x=(np.asarray(amps)/max(amps))**2
    y=np.asarray(K2,float)
    return np.linalg.lstsq(np.c_[np.ones(len(x)),x],y,rcond=None)[0][0]
