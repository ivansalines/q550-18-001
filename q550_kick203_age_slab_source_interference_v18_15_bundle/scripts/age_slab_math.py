#!/usr/bin/env python3
"""Age-slab algebra for a cumulative source in the initial co-moving 10D gauge.

With A(z)=s int_0^z g(u)du and age a=z-u:
 S[a,b]=A(z-a)-A(z-b), 0<=a<b<=z. e=A(1), Q=I-ee^T/(e.e),
 w=QA(z), v=w/|w|, q_i=QS_i, and sum_i q_i=w.
 L_i=q_i.v, T_i=q_i-L_i*v. sum_i L_i=|w|, sum_i T_i=0.
 Signed L_i, angles, orthogonality and destructive interference are Euclidean-
 gauge dependent. The additive identities are algebraic; this diagnostic does
 not rerun a different ODE history or establish constitutive hysteresis.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator

def orth(e,arr):
    e=np.asarray(e,float);u=e/np.linalg.norm(e)
    arr=np.asarray(arr,float)
    return arr-np.sum(arr*u,axis=-1)[...,None]*u

def slab_edges(z,width):
    if z<=0 or width<=0:raise ValueError('positive z and width required')
    edges=list(np.arange(0,z,width))
    if not edges or edges[0]!=0:edges.insert(0,0.)
    if z-edges[-1]<1e-10:edges[-1]=z
    else:edges.append(z)
    return np.array(edges,float)

def slabs(ai,e,z,width):
    edges=slab_edges(z,width)
    vals=ai(z-edges[:-1])-ai(z-edges[1:])
    return edges,orth(e,vals)

def metrics(ai,e,z,width,gi=None,s=None):
    edge,q=slabs(ai,e,z,width)
    w=orth(e,ai(z));rho=float(np.linalg.norm(w))
    if rho<1e-9:raise ValueError('degenerate full-history residual')
    v=w/rho
    axial=q@v
    lateral=q-axial[:,None]*v
    norms=np.linalg.norm(q,axis=1)
    plus=float(np.maximum(axial,0).sum());minus=float(-np.minimum(axial,0).sum())
    gross=plus+minus
    sum_norm=float(norms.sum())
    signed_identity=float(abs(axial.sum()-rho))
    closure=float(np.linalg.norm(q.sum(axis=0)-w))
    pair=q@q.T
    upper=pair[np.triu_indices(len(q),k=1)]
    pair_identity=float(abs(rho**2-(norms@norms+2*upper.sum())))
    reconstructed=np.zeros_like(w)
    angles=[]
    for k in range(len(q)):
        reconstructed+=q[k]
        nr=np.linalg.norm(reconstructed)
        angles.append(float(np.degrees(np.arccos(np.clip(abs(reconstructed@v/max(nr,1e-30)),0,1)))) if nr>1e-9 else 90.)
    def stable_horizon(thresh):
        # Earliest age after which ALL older prefixes remain within threshold.
        violations=[j for j,t in enumerate(angles) if t>thresh]
        return float(edge[(max(violations)+1) if violations else 1])
    neg_weights=np.maximum(-axial,0)
    pos_weights=np.maximum(axial,0)
    centers=(edge[:-1]+edge[1:])/2
    oldmask=edge[:-1]>=2.-1e-12
    old_w=q[oldmask].sum(axis=0)
    recent_w=w-old_w
    results={'z':float(z),'width':float(width),'rho':rho,'nslabs':len(q),
      'positive_axial':plus,'negative_axial':minus,'gross_axial':gross,
      'destructive_fraction':1.-rho/gross if gross>0 else 0.,
      'vector_cancellation_fraction':1.-rho/sum_norm if sum_norm>0 else 0.,
      'norm_slab_sum':sum_norm,'closure_error':closure,'axial_closure_error':signed_identity,
      'pair_identity_error':pair_identity,'age_horizon_10deg':stable_horizon(10.),
      'age_horizon_25deg':stable_horizon(25.),
      'negative_age_centroid':float(neg_weights@centers/neg_weights.sum()) if neg_weights.sum()>0 else None,
      'positive_age_centroid':float(pos_weights@centers/pos_weights.sum()) if pos_weights.sum()>0 else None,
      'old2_norm':float(np.linalg.norm(old_w)),
      'recent2_norm':float(np.linalg.norm(recent_w)),
      'old_recent_dot':float(old_w@recent_w),
      'old_recent_cos':float(old_w@recent_w/max(np.linalg.norm(old_w)*np.linalg.norm(recent_w),1e-30)),
      'independent_source_quadrature_relative_error':None,
      'max_independent_slab_quadrature_relative_error':None}
    if gi is not None and s is not None:
        integral=gi.antiderivative()
        independent=s*(integral(z-edge[:-1])-integral(z-edge[1:]))
        qi=orth(e,independent)
        err=np.linalg.norm(q-qi,axis=1)
        results['independent_source_quadrature_relative_error']=float(np.linalg.norm(qi.sum(axis=0)-w)/rho)
        results['max_independent_slab_quadrature_relative_error']=float(max(err)/rho)
    rows=[]
    for k,(a,b,vec,l,t,mag) in enumerate(zip(edge[:-1],edge[1:],q,axial,lateral,norms)):
        out=w-vec;outn=float(np.linalg.norm(out))
        if outn>1e-10:angle=float(np.degrees(np.arccos(np.clip(abs(out@v/outn),0.,1.))))
        else:angle=None
        rows.append({'age_start':float(a),'age_end':float(b),
                     'source_xi_start':float(z-b),'source_xi_end':float(z-a),
                     'signed_axial':float(l),'signed_axial_over_full':float(l/rho),
                     'tangential_norm':float(np.linalg.norm(t)),
                     'tangential_over_full':float(np.linalg.norm(t)/rho),
                     'norm':float(mag),'norm_over_full':float(mag/rho),
                     'leave_one_out_line_angle_deg':angle,
                     'role':'opposing' if l<0 else 'reinforcing',
                     'age_mid':float(centers[k])})
    pairs=[]
    for i in range(len(q)):
        for j in range(i+1,len(q)):
            if pair[i,j]<0:
                pairs.append({'age_a_start':float(edge[i]),'age_a_end':float(edge[i+1]),
                              'age_b_start':float(edge[j]),'age_b_end':float(edge[j+1]),
                              'dot':float(pair[i,j]),'interference_over_rho2':float(pair[i,j]/rho**2)})
    pairs.sort(key=lambda r:r['dot'])
    return results,rows,pairs
