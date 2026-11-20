"""Pure full-state algebra for the forward write / tangent-matched echo audit.

No geometric-descriptor 11D array is ever treated as a material-state vector.
"""
from __future__ import annotations
import numpy as np

DELAYS_XI=(0.25,0.5,1.0,2.0,4.0,8.0)
READ_AGES_XI=(0.0,0.25,0.5,1.0,2.0)
AMPS=(1e-8,3e-8,1e-7)
TOLS={'standard':(1e-12,1e-14),'tight':(3e-13,3e-15)}


def compensate(state_write, reference, tangent_unit, amp, sign):
    """Apply minus the LINEAR propagated write, leaving nonlinear memory.

    `tangent_unit=Phi(t1,t0) v` has not been normalized a second time.
    """
    return np.asarray(state_write)-float(sign)*float(amp)*np.asarray(tangent_unit)


def transport_from_origin(P_read, P_lag, residual, cond_gate=1e10):
    """Full-state linear transport of residual, conditioning-gated.

    Algebra: Phi(read, lag) residual = P(read) solve(P(lag), residual).
    """
    if not np.all(np.isfinite(P_lag)) or not np.all(np.isfinite(P_read)):
        return np.full_like(residual, np.nan), np.inf
    c=float(np.linalg.cond(P_lag))
    if not np.isfinite(c) or c>cond_gate:
        return np.full_like(residual,np.nan),c
    return np.asarray(P_read)@np.linalg.solve(P_lag,residual),c


def signed_even_odd(plus,minus):
    """Even / odd residual under reversal of FIRST kick sign."""
    return 0.5*(np.asarray(plus)+np.asarray(minus)),0.5*(np.asarray(plus)-np.asarray(minus))


def unit_angle(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    na=float(np.linalg.norm(a));nb=float(np.linalg.norm(b))
    if min(na,nb)<1e-25:return np.nan
    return float(np.rad2deg(np.arccos(np.clip(a@b/(na*nb),-1,1))))


def safe_ratio(a,b,floor=1e-22):
    return float(a/b) if np.isfinite(a) and np.isfinite(b) and abs(b)>floor else float('nan')
