"""Dimensionless odd-polynomial fits for frozen-coordinate 11D scaling audits.

Never mistake exact fit closure for a proof of physical constitutive law.
"""
from __future__ import annotations
import numpy as np

AMPS = np.array([5e-9, 1e-8, 1.5e-8, 2e-8, 3e-8], dtype=float)
EPS_UNIT = 1e-8


def fit_odd(amplitudes, values, powers=(1, 3, 5), eps_unit=EPS_UNIT):
    """Fit amplitudes x 11D or x N x 11D arrays: values=sum(c_p*(eps/eps_unit)^p).

    Returns dimensionless coefficients and residuals using unweighted LS, plus
    leave-one-amplitude-out predictions only if each reduced fit has enough rows.
    """
    amplitudes = np.asarray(amplitudes, float)
    values = np.asarray(values, float)
    powers = tuple(powers)
    if amplitudes.ndim != 1 or values.shape[0] != amplitudes.size:
        raise ValueError('amplitude/value count mismatch')
    if len(amplitudes) < len(powers) or not np.all(amplitudes > 0):
        raise ValueError('insufficient strictly positive amplitudes')
    t = amplitudes / eps_unit
    X = np.stack([t ** p for p in powers], axis=1)
    cond = float(np.linalg.cond(X))
    flat = values.reshape(len(amplitudes), -1)
    coeff_flat, _, _, _ = np.linalg.lstsq(X, flat, rcond=None)
    predicted_flat = X @ coeff_flat
    loo = np.full_like(flat, np.nan)
    if len(amplitudes) - 1 >= len(powers):
        for i in range(len(amplitudes)):
            mask = np.arange(len(amplitudes)) != i
            cf, _, _, _ = np.linalg.lstsq(X[mask], flat[mask], rcond=None)
            loo[i] = X[i] @ cf
    return {'coefficients': coeff_flat.reshape((len(powers),) + values.shape[1:]),
            'predicted': predicted_flat.reshape(values.shape),
            'residual': (flat - predicted_flat).reshape(values.shape),
            'loo_residual': (flat - loo).reshape(values.shape),
            'design_condition': cond, 'powers': powers}


def samplewise_exponent(a, lo, hi):
    a = np.asarray(a, float)
    good = (lo > 0) & (hi > 0)
    return np.where(good, np.log(np.maximum(hi, 1e-300)/np.maximum(lo, 1e-300)) / np.log(a), np.nan)


def finite_or_none(v):
    v = float(v)
    return v if np.isfinite(v) else None


def safe_median(x):
    y = np.asarray(x, float)
    y = y[np.isfinite(y)]
    return finite_or_none(np.median(y)) if y.size else None
