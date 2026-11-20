#!/usr/bin/env python3
"""Q550 v18.12: residual motion of phase-fixed memory plane after common transport."""
from __future__ import annotations
import argparse,csv,json,hashlib,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from checkpoint_io_v18_11 import load
from residual_plane_math import source_geometry,principal_line_angle,velocity_from_projector,max_tensor_transform_error

def write_csv(path, rows):
    if rows:
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)

def audit_tol(name,d):
    x=np.asarray(d['xi_nodes']); ids=np.where(x>1.+1e-12)[0]
    F=d['P_full11D'][:,:10,:10]
    J=d['J_full11D'][:,:10,:10]
    E=d['old_source_material_11D'][:,:10]
    K=d['quadratic_total_material_11D'][:,:10]
    A=d['integral_source_frame_11D'][:,:10]
    g=d['source_integrand_11D'][:,:10]
    slope=float(d['slope'])
    e=A[80].copy()
    assert F.shape==(641,10,10) and x.shape==(641,)
    cond=np.linalg.cond(F[ids]);
    if not np.isfinite(cond).all() or np.max(cond)>1e8:raise RuntimeError('Unsafe inverse-cocycle conditioning')
    Ein=np.linalg.solve(F[ids],E[ids,:,None])[:,:,0]
    Kin=np.linalg.solve(F[ids],K[ids,:,None])[:,:,0]
    ee=float(np.max(np.linalg.norm(Ein-e,axis=1)/max(np.linalg.norm(e),1e-30)))
    ek=float(np.max(np.linalg.norm(Kin-A[ids],axis=1)/np.maximum(np.linalg.norm(A[ids],axis=1),1e-30)))
    if max(ee,ek)>1e-8:raise RuntimeError(f'Inverse transport identity failed: {ee} {ek}')
    geom=source_geometry(e,A[ids],g[ids],slope)
    # Independent check from transported vectors WITHOUT referring to archived A:
    backgeom=source_geometry(e,Kin,g[ids],slope)
    v=geom['v']; rate=geom['turn_rate']; dv=geom['vdot']
    if np.max(abs(rate-backgeom['turn_rate']))>1e-8:raise RuntimeError('Transported reconstruction differs')
    # Analytic residual transport: derivative of projectors, source term only.
    projrate=velocity_from_projector(v,dv)
    analytic_id=float(np.max(abs(np.linalg.norm(projrate,axis=(1,2))/np.sqrt(2)-rate)))
    # Finite difference checks on the full A(x) prior to the xi>1 cut.
    a_fd=np.gradient(A,x,axis=0,edge_order=2)
    fd_A=np.linalg.norm(a_fd[ids[1:-1]]-slope*g[ids[1:-1]],axis=1)/np.maximum(np.linalg.norm(slope*g[ids[1:-1]],axis=1),1e-12)
    v_fd=np.gradient(v,x[ids],axis=0,edge_order=2)
    fd_v=np.linalg.norm(v_fd[1:-1]-dv[1:-1],axis=1)/np.maximum(rate[1:-1],1e-12)
    # Continuous travel distance: source derivative + trapezoidal quadrature;
    # compare with sum of adjacent angles (independent of derivative formula).
    length_integral=float(np.trapezoid(rate,x[ids]));
    ds=np.arccos(np.clip(np.sum(v[1:]*v[:-1],axis=1),-1,1))
    length_steps=float(np.sum(ds))
    direct_angle=float(np.arccos(np.clip(np.dot(v[0],v[-1]),-1,1)))
    endpoint_line_angle=float(np.arccos(np.clip(abs(np.dot(v[0],v[-1])),0,1)))
    # The source-forcing direction can have longitudinal and transverse components.
    # A sign of its radial part is not a torsional handedness claim.
    radiusdot=geom['radius_rate'];
    rad_crossings=[]
    for i in range(len(radiusdot)-1):
        if radiusdot[i]*radiusdot[i+1]<0:
            r=x[ids][i]+(x[ids][i+1]-x[ids][i])*abs(radiusdot[i])/(abs(radiusdot[i])+abs(radiusdot[i+1]))
            rad_crossings.append(float(r))
    # Numerical GL(10) covariance of antisymmetric residual forcing.
    T=np.eye(10)
    for j in range(10):T[j,j]=np.exp(-0.8+1.6*j/9)
    for j in range(9):T[j,j+1]=0.025
    cov_err=max(max_tensor_transform_error(T,e,A[i],g[i]) for i in ids[::79])
    if cov_err>1e-11:raise RuntimeError('Residual wedge covariance failed')
    # Record near previous Euclidean zeroes as observers, not invariant events.
    near=[]
    for root in (4.1351989579,4.9973365558,7.1744358136):
        ix=int(np.argmin(abs(x[ids]-root)))
        near.append({'approx_old_euclidean_root':root,'sample_xi':float(x[ids][ix]),'residual_speed':float(rate[ix]),'tangent_source_fraction':float(geom['tangent_source_fraction'][ix]),'radius_rate':float(radiusdot[ix])})
    rows=[]
    for h,i in enumerate(ids):
        rows.append({'tolerance':name,'xi':float(x[i]),'cond_phi':float(cond[h]),
          'source_frame_radius_e_perp_gauge':float(geom['radius'][h]),
          'source_frame_radial_rate_gauge':float(radiusdot[h]),
          'residual_plane_turn_speed_gauge':float(rate[h]),
          'source_turn_outside_plane_gauge':float(geom['source_outside_plane_norm'][h]),
          'source_turn_fraction_gauge':float(geom['tangent_source_fraction'][h]),
          'cumulative_residual_turn_length_to_here_gauge':float(np.trapezoid(rate[:h+1],x[ids][:h+1])) if h else 0.,
          'residual_projector_rate_frobenius':float(np.linalg.norm(projrate[h]))})
    out={'tolerance':name,'samples':len(ids),'xi_range':[float(x[ids[0]]),float(x[ids[-1]])],
      'inverse_transport_E_relative_max':ee,'inverse_transport_K_relative_max':ek,
      'phi_condition_max':float(max(cond)),'phi_condition_median':float(np.median(cond)),
      'source_frame_A_derivative_fd_relative_median':float(np.median(fd_A)),
      'source_frame_A_derivative_fd_relative_max':float(np.max(fd_A)),
      'residual_line_derivative_fd_relative_median':float(np.median(fd_v)),
      'residual_line_derivative_fd_relative_max':float(np.max(fd_v)),
      'projector_speed_identity_max_error':analytic_id,
      'residual_turn_rate_min':float(min(rate)), 'residual_turn_rate_median':float(np.median(rate)),
      'residual_turn_rate_max':float(max(rate)),
      'residual_turn_rate_positive_count_above_1e-7':int(np.sum(rate>1e-7)),
      'residual_line_path_length_from_rate_rad':length_integral,
      'residual_line_path_length_from_steps_rad':length_steps,
      'residual_endpoint_vector_angle_rad':direct_angle,
      'residual_endpoint_unoriented_line_angle_rad':endpoint_line_angle,
      'accumulated_path_to_endpoint_angle_ratio':float(length_integral/max(direct_angle,1e-12)),
      'radial_rate_crossings_linear_interp_xi':rad_crossings,
      'min_source_out_of_plane_fraction_gauge':float(np.min(geom['tangent_source_fraction'])),
      'median_source_out_of_plane_fraction_gauge':float(np.median(geom['tangent_source_fraction'])),
      'covariant_source_wedge_GL10_relative_max':cov_err,
      'old_metric_roots_diagnostics':near,
      'status':'NUMERICAL_COMOVING_RESIDUAL_DETECTED_GAUGE_DECLARED_NO_PHYSICAL_HYSTERESIS_CERTIFICATE'}
    return out,rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=None);arg=ap.parse_args()
    target=Path(arg.output) if arg.output else ROOT/'results';target.mkdir(parents=True,exist_ok=True)
    data,provenance=load()
    reports=[];rows=[]
    for name in ('standard','tight'):
        summary,r=audit_tol(name,data[name]);reports.append(summary);rows+=r
    std,tight=reports
    checks={k:float(abs(std[k]-tight[k])) for k in ('residual_turn_rate_median','residual_line_path_length_from_rate_rad','residual_endpoint_vector_angle_rad')}
    total={'version':'Q550-v18.12-comoving-residual-plane','input_provenance':provenance,
           'by_tolerance':reports,'standard_tight_absolute_disagreements':checks,
           'mathematical_identity':'Phi^-1 E=e_0, Phi^-1 K=A; A_xi=s*g; after homogeneous transport removal the residual plane span(e_0,A) turns iff e_0 wedge A wedge g !=0 (on full rank interval).',
           'metric_and_observer_scope':'Speeds/path-lengths/angles use the Euclidean metric at the INITIAL 10D frame and are not metric-free. The wedge-zero criterion and transformed 2-form covariance are coordinate-invariant under GL(10).',
           'status':'NO_INTRINSIC_MECHANICAL_METRIC_OR_NONLINEAR_CONSTRAINT_CERTIFIED',
           'limits':['Computing accumulated source A already contains history by definition; nonzero turn does NOT alone prove physical hysteresis or material torsion.','Not an interval proof between sampled nodes.','Phase leaf is time dependent, not an autonomous 10D mechanical constraint.']}
    (target/'summary.json').write_text(json.dumps(total,indent=2,allow_nan=False)+'\n')
    write_csv(target/'residual_plane_profile.csv',rows)
    summary=std
    report=f'''# Q550 v18.12 — Co-moving residual plane after common transport

## Result

- Source verified through v18.11 -> v18.10 -> v18.9 -> v18.7 -> v18.6, SHA256.
- Common homogeneous transport removed by the inverse 10x10 fixed-phase propagator.
- Both old and quadratic source identities reconstructed at 1e-8 relative tolerance.
- Source-frame residual plane turn rate median **{summary['residual_turn_rate_median']:.9g}** per xi, min **{summary['residual_turn_rate_min']:.9g}**, max **{summary['residual_turn_rate_max']:.9g}** (Euclidean initial frame ONLY).
- Positive turn rates at **{summary['residual_turn_rate_positive_count_above_1e-7']}/{summary['samples']}** sampled nodes beyond xi=1.
- Integrated turn path **{summary['residual_line_path_length_from_rate_rad']:.9g}** rad versus endpoint angle **{summary['residual_endpoint_vector_angle_rad']:.9g}** rad. The discrepancy indicates bending/backtracking in the declared gauge, not physical hysteresis.
- Max condition number for Phi **{summary['phi_condition_max']:.6g}**; maximum 10D inverse transport errors E/K **{summary['inverse_transport_E_relative_max']:.3g}/{summary['inverse_transport_K_relative_max']:.3g}**.
- Independent finite differences of v(xi) agree with sourced analytical derivative at median relative error **{summary['residual_line_derivative_fd_relative_median']:.3g}**.
- GL(10) covariance of source 2-form: relative error **{summary['covariant_source_wedge_GL10_relative_max']:.3g}**.

## Mechanism

Let e=A(1), K=Phi A, E=Phi e, A_xi=s*g and define ehat=e/||e||, w=(I-ehat ehat^T)A, v=w/||w||.

    v_xi = s * (I-ehat ehat^T-v v^T) g / ||w||.

The homogeneous J transport is no longer present. A nonzero residual means fresh source drives the plane outside its previously occupied 2-plane, **within the chosen fixed initial-frame gauge**.

**Strict boundary**: a changing accumulated source A necessarily remembers earlier integration, but is not automatically hysteresis (multi-valued response to history, loop memory, irreversibility) and is not a certificate of physical body torsion. No source-derived physical metric G(y) or constraints C(y) were provided. The path length is metric dependent, although nonzero raw 3-blade e wedge A wedge g is algebraically invariant under invertible reparametrization.

## Next experimental gate

Compare histories that reach the same verified reduced material state/phase with different kicks, then measure unequal source-frame residual planes/response; alternatively construct a physical energy metric and constraint manifold from the real ODE. Do NOT infer hysteresis merely from cumulative source turning.
'''
    (target/'REPORT.md').write_text(report)
    print(json.dumps({'status':'FROZEN_AUDIT_PASSED','std':{k:std[k] for k in ('residual_turn_rate_median','residual_turn_rate_max','residual_line_path_length_from_rate_rad','residual_endpoint_vector_angle_rad','phi_condition_max')},'standard_tight_abs_diff':checks},indent=2))

if __name__=='__main__':main()
