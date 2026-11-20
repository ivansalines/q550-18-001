"""Jet observability for phase-fixed material variational equation z'=A(xi) z.

For fixed H, M0=H; M1=HA; M2=H(A'+A^2);
M3=H(A''+2A'A+AA'+A^3). Derivatives are taken ALONG base orbit.
"""
import numpy as np
import math

def derivatives(arr,index,dx,jump=1):
    x=np.asarray(arr,float);h=dx*jump;j=int(index);s=int(jump)
    if j-2*s<0 or j+2*s>=len(x):raise ValueError('Stencil outside frozen domain')
    A=x[j]; p=x[j+s];m=x[j-s];pp=x[j+2*s];mm=x[j-2*s]
    D1=(mm-8*m+8*p-pp)/(12*h)
    D2=(-mm+16*m-30*A+16*p-pp)/(12*h*h)
    return A,D1,D2

def motion_matrices(A,D1,D2):
    I=np.eye(A.shape[0]);return [I,A,D1+A@A,D2+2*D1@A+A@D1+A@A@A]

def observability(H,M,tau=.1):
    return [np.vstack([tau**k*(H@M[k]) for k in range(i+1)]) for i in range(len(M))]

def numerical_rank(O,rtol=1e-8):
    s=np.linalg.svd(O,compute_uv=False)
    rank=int(np.sum(s>max(1e-12,rtol*s[0])))
    return rank,s

def null_vectors(O,rtol=1e-8):
    u,s,vt=np.linalg.svd(O,full_matrices=True)
    rank=int(np.sum(s>max(1e-12,rtol*s[0])))
    return vt[rank:].T

def canonical(v):
    u=np.asarray(v,float).copy();u/=np.linalg.norm(u)
    idx=np.argmax(abs(u))
    return u if u[idx]>=0 else -u

def scaled_taylor(H,M,v,t):
    return sum(t**k/math.factorial(k)*(H@M[k]@v) for k in range(len(M)))
