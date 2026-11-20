"""Linear noise and invertible-gauge tests for local jet maps.
Physical jet measurements are *not* being simulated: noise applies directly to
scaled jet coordinates, with an explicitly hypothetical white covariance.
"""
import numpy as np

def orthogonal(n,seed=1821):
    r=np.random.default_rng(seed).normal(size=(n,n))
    q,r=np.linalg.qr(r)
    q=q@np.diag(np.sign(np.diag(r)))
    return q

def gauge_matrices(n,m,k,seed=1821):
    rng=np.random.default_rng(seed)
    # A constant invertible material state gauge and output gauge, both explicit.
    scales=np.geomspace(.05,20.,n)
    T=orthogonal(n,seed)@np.diag(scales)
    out_scales=np.geomspace(.2,5.,m)
    G=orthogonal(m,seed+1)@np.diag(out_scales)
    block=np.kron(np.eye(k+1),G)
    return T,block

def svals(M):return np.linalg.svd(M,compute_uv=False)

def oracle_stats(O,levels,seed=1821,trials=256):
    u,s,vt=np.linalg.svd(O,full_matrices=False)
    if len(s)<10 or s[-1]<1e-12:raise RuntimeError('Rank-deficient final jet map')
    vweak=vt[-1]
    rng=np.random.default_rng(seed)
    out=[]
    for noise in levels:
        N=rng.normal(size=(O.shape[0],trials))*noise
        # Reconstruct a unit perturbation in the weakest (euclidean) direction.
        estimate=vweak[:,None]+vt.T@((u.T@N)/s[:,None])
        err=np.linalg.norm(estimate-vweak[:,None],axis=0)
        # The matched-filter detectability of the weakest unit state perturbation.
        out.append(dict(noise_sigma=noise,weak_gain=float(s[-1]),worst_case_single_sigma=noise/s[-1],
            empirical_rmse=float(np.sqrt(np.mean(err**2))),median_error=float(np.median(err)),
            p95_error=float(np.quantile(err,.95)),weak_direction_snr=float(s[-1]/noise),
            analytical_expected_rmse=float(noise*np.sqrt(np.sum(1./s**2))),
            distinguishability_snr_threshold_1=bool(s[-1]>=noise)))
    return out

def gauge_test(O,m,k,seed=1821):
    n=O.shape[1];T,B=gauge_matrices(n,m,k,seed)
    Ti=np.linalg.inv(T)
    Op=B@O@Ti
    rng=np.random.default_rng(seed+9)
    y=rng.normal(size=n);e=rng.normal(size=O.shape[0])*.0001
    raw=O@y+e
    transformed=B@raw
    # Gauge-aware whitening in output space; infer physical (old-coordinate) state.
    normalized=np.linalg.solve(B,transformed)
    yhat1=np.linalg.lstsq(O,raw,rcond=None)[0]
    yhat2=Ti@np.linalg.lstsq(np.linalg.solve(B,Op),normalized,rcond=None)[0]
    target_pred=O@yhat1; alt_pred=O@yhat2
    return dict(physical_state_error=float(np.linalg.norm(yhat1-yhat2)),
       physical_prediction_error=float(np.linalg.norm(target_pred-alt_pred)),
       raw_cond=float(np.linalg.cond(O)),naively_transformed_cond=float(np.linalg.cond(Op)),
       cond_T=float(np.linalg.cond(T)),cond_B=float(np.linalg.cond(B)),
       mapped_residual_error=float(np.linalg.norm(raw-np.linalg.solve(B,transformed))))
