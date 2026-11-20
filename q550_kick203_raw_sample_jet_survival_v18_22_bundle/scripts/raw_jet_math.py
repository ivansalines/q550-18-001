"""Compare equally noisy direct time observations to reconstructed local derivative jets.

Explicit covariances use the actual polynomial pseudoinverse, including correlations
between polynomial coefficients. No derivative measurement is presumed available.
"""
import math
import numpy as np

def basis_times(halfwidth,sample_count):
    if sample_count<9 or sample_count%2!=1:raise ValueError('Odd number of samples >=9 required')
    return np.linspace(-halfwidth,halfwidth,sample_count)

def polynomial_operator(offsets,degree,halfwidth):
    offsets=np.asarray(offsets,float)
    V=np.vander(offsets/halfwidth,N=degree+1,increasing=True)
    if len(offsets)<=degree:raise ValueError('Polynomial not overdetermined')
    B=np.linalg.pinv(V,rcond=1e-12)
    return V,B

def jet_compression(offsets,halfwidth,degree,output_dimension,jet_matrices):
    V,B=polynomial_operator(offsets,degree,halfwidth)
    # The first four polynomial coefficients are estimates for scaled Taylor jets.
    C=B[:4,:]@B[:4,:].T
    L=np.linalg.cholesky(C)
    whiten=np.linalg.inv(L)
    Hjet=np.vstack([halfwidth**k/math.factorial(k)*jet_matrices[k] for k in range(4)])
    W=np.kron(whiten,np.eye(output_dimension))
    return B[:4],C,W,W@Hjet,Hjet

def direct_observation(P,center_index,indices,H):
    inv=np.linalg.inv(P[center_index])
    blocks=[H@(P[i]@inv) for i in indices]
    return np.vstack(blocks)

def continuous_observation(propagator_spline,P_center,times,H):
    """Evaluate propagator at exact requested times via cubic archived-data interpolation.

    This is a frozen-checkpoint interpolation, not a fresh nonlinear ODE integration.
    """
    inv=np.linalg.inv(P_center)
    return np.vstack([H@(np.asarray(propagator_spline(float(t)))@inv) for t in times])

def direct_estimator(O):
    # Full-state least-squares reconstruction in the specified fixed Euclidean gauge.
    return np.linalg.pinv(O,rcond=1e-12)

def jet_estimator(Ojet):
    return np.linalg.pinv(Ojet,rcond=1e-12)

def compressed_observation(noisy,output_dimension,B4,W):
    y=np.asarray(noisy,float).reshape(-1,output_dimension)
    coef=B4@y
    return W@coef.reshape(-1)

def test_covariance(C,output_dimension,W):
    full=np.kron(C,np.eye(output_dimension))
    return float(np.linalg.norm(W@full@W.T-np.eye(4*output_dimension),ord=np.inf))

def nominal_jet_bias(O,v,B4,W,output_dimension,Ojet):
    compressed=compressed_observation(O@v,output_dimension,B4,W)
    return compressed-Ojet@v
