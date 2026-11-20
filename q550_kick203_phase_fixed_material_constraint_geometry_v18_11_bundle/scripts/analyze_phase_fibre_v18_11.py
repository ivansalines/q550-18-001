#!/usr/bin/env python3
"""Frozen-data v18.11 phase-fixed 10D geometry without physical constraint claims."""
from __future__ import annotations
import csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from checkpoint_io_v18_11 import load
from phase_fibre_math_v18_11 import phase_sector_checks,phase_cocycle_rhs,tangential_plane_rates,principal_step_angle

def save_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def relerr(num,den):
    return np.linalg.norm(num.reshape(len(num),-1),axis=1)/np.maximum(np.linalg.norm(den.reshape(len(den),-1),axis=1),1e-12)

def per_tolerance(label,d):
    x=d['xi_nodes'];J=d['J_full11D'];P=d['P_full11D'];g=d['source_integrand_11D'];E=d['old_source_material_11D'];K=d['quadratic_total_material_11D'];A=d['integral_source_frame_11D'];slope=d['slope']
    if J.shape!=(641,11,11) or P.shape!=(641,11,11) or np.max(np.abs(np.diff(x)-.0125))>1e-10:raise RuntimeError('v18.7 source grid regression')
    check=phase_sector_checks(J,P,[g,E,K,A]);zero_max=max(check[k] for k in ('max_jacobian_phase_row','max_cocycle_phase_row_to_z','max_cocycle_phase_diagonal_from_one','max_source_phase_component'))
    if zero_max>1e-9:raise RuntimeError('phase-fixed sector is NOT preserved: '+str(check))
    F=P[:,:10,:10];q=P[:,:10,10]
    Fder,qder=phase_cocycle_rhs(J,P,slope)
    Ffd=np.gradient(F,x,axis=0,edge_order=2);qfd=np.gradient(q,x,axis=0,edge_order=2)
    ef=relerr(Ffd[1:-1]-Fder[1:-1],Fder[1:-1]);eq=relerr(qfd[1:-1]-qder[1:-1],qder[1:-1])
    # The material vectors E,K are in the fixed-phase invariant tangent bundle.
    # Geometry of span(E,K) is measured using Euclidean 10D only as a DECLARED gauge.
    rows=[];D=[];nfd=[];angles=[];samples=[]
    for i in range(81,641):
        dt,ds,db,Q,sv=tangential_plane_rates(E[i],K[i],J[i],P[i],g[i],slope)
        defect=np.linalg.norm(dt-ds-db)
        D.append(defect)
        U=np.column_stack((E[i,:10],K[i,:10]));pi=Q@Q.T
        samples.append(pi)
        ang=principal_step_angle(U,np.column_stack((E[i+1,:10],K[i+1,:10]))) if i<640 else np.nan
        angles.append(ang)
        rows.append({'tolerance':label,'xi':float(x[i]),
            'plane_sv_min':float(sv[-1]),'plane_condition_gauge':float(sv[0]/sv[-1]),
            'grassmann_total_speed_euclidean_gauge':float(np.linalg.norm(dt)/np.sqrt(2)),
            'grassmann_transport_speed_euclidean_gauge':float(np.linalg.norm(ds)/np.sqrt(2)),
            'grassmann_birth_speed_euclidean_gauge':float(np.linalg.norm(db)/np.sqrt(2)),
            'grassmann_transport_birth_inner_euclidean_gauge':float(np.sum(ds*db)/2),
            'grassmann_decomposition_residual':float(defect),
            'plane_next_step_max_principal_angle_rad_euclidean_gauge':float(ang)})
    # compare finite difference of projector with analytical equation away from initial degeneracy;
    # derivative fd requires neighboring valid projectors.
    numerical=[];analytic=[]
    for j in range(1,len(samples)-1):
        i=81+j;Dpi=(samples[j+1]-samples[j-1])/(x[i+1]-x[i-1])
        dt,_,_,_,_=tangential_plane_rates(E[i],K[i],J[i],P[i],g[i],slope)
        if x[i]>=1.10:
            numerical.append(np.linalg.norm(Dpi-dt));analytic.append(max(1e-12,np.linalg.norm(dt)))
    ratio=np.asarray(numerical)/np.asarray(analytic)
    summary={'tolerance':label,'n_samples_xi_gt_one':len(rows),'clock_coordinate_index_zero_based':10,'phase_sector':check,
      'phase_fixed_subbundle_invariant_on_archived_samples':True,
      'phase_fixed_ten_d_not_automatically_autonomous':True,
      'ten_d_phi_variational_fd_relative_median':float(np.median(ef)),
      'ten_d_phi_variational_fd_relative_max':float(np.max(ef)),
      'phase_to_shape_driven_q_fd_relative_median':float(np.median(eq)),
      'phase_to_shape_driven_q_fd_relative_max':float(np.max(eq)),
      'plane_min_singular_value_gauge':float(min(v['plane_sv_min'] for v in rows)),
      'plane_max_condition_gauge':float(max(v['plane_condition_gauge'] for v in rows)),
      'grassmann_projector_derivative_relative_median_xi_gt_1_1':float(np.median(ratio)),
      'grassmann_projector_derivative_relative_max_xi_gt_1_1':float(np.max(ratio)),
      'grassmann_rate_transport_plus_birth_max_algebraic_error':float(max(D)),
      'max_one_step_principal_angle_rad_gauge':float(np.nanmax(angles)),
      'median_grassmann_total_speed_gauge':float(np.median([z['grassmann_total_speed_euclidean_gauge'] for z in rows])),
      'median_grassmann_transport_speed_gauge':float(np.median([z['grassmann_transport_speed_euclidean_gauge'] for z in rows])),
      'median_grassmann_birth_speed_gauge':float(np.median([z['grassmann_birth_speed_euclidean_gauge'] for z in rows])),
      'phase_clock_is_time_dependent_kinematic_leaf_not_material_constraint':True,
      'identified_nonlinear_material_constraint':False,'identified_physical_metric':False}
    return summary,rows

def main():
    (ROOT/'results').mkdir(exist_ok=True)
    dat,provenance=load(); rows=[];summ=[]
    for name,d in dat.items():
        s,table=per_tolerance(name,d);summ.append(s);rows.extend(table)
    idx={s['tolerance']:s for s in summ}
    reg={k:abs(idx['standard'][k]-idx['tight'][k]) for k in ('median_grassmann_total_speed_gauge','median_grassmann_transport_speed_gauge','median_grassmann_birth_speed_gauge')}
    result={'version':'2026-10-18-Q550-v18.11-phase-fixed-material-transport', 'historical_provenance':provenance,
      'per_tolerance':summ,'standard_tight_median_geometry_disagreement':reg,
      'what_is_proved_exactly':'Algebraic invariance of the fixed-phase tangent sector conditional on the observed block-zero Jacobian structure; samplewise ODE checks are numerical.',
      'constraint_status':'NO_PHYSICAL_HOLONOMIC_CONSTRAINT_IDENTIFIED',
      'metric_status':'EUCLIDEAN_10D_GAUGE_ONLY_NOT_PHYSICAL',
      'higher_order_constraint_curvature_status':'NOT_DERIVABLE_WITHOUT_SOURCE_VERIFIED_C_AND_D2C',
      'limitation':'No interval certification of all phase-space states; phase fibre is a time-dependent kinematic leaf, not a fixed flow invariant hypersurface. The Grassmann speeds and principal angles vary with the metric.'}
    (ROOT/'results/summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    save_csv(ROOT/'results/phase_fixed_plane_profile.csv',rows)
    report=['# Q550 v18.11 — Phase-fixed material transport & constraint geometry','',
      '## Verified source','',f"Nested v18.10 SHA256 {provenance['source_v18_10_sha256']}. SHA provenance checked through historical v18.6.",
      '', '## Mathematical structure','',
      'State y=(z,theta), theta_dot=2pi along original sampled trajectory, J=[[A,b],[0,0]]; the tangent bundle dtheta=0 is preserved by the linear cocycle P=[[Phi,q],[0,1]].',
      'The 10D model is z_dot=f_z(z,theta0+2pi*(s-s0)), GENERALLY NONAUTONOMOUS. q_dot=A q+b is the effect of changing the initial phase. Do not identify the phase time-leaf with a material constraint.',
      'The algebraic rank-2 plane span(E,K) is preserved as 10D-vectors, but its Euclidean projector is not GL-invariant. Projector speed decomposes into homogeneous transport and fresh-source birth. No metric-independent angular rate is claimed.',
      '', '## Numbers','']
    for t in summ:
        report.extend([f"- {t['tolerance']}: phase structural max defect {max(t['phase_sector'][k] for k in ('max_jacobian_phase_row','max_cocycle_phase_row_to_z','max_cocycle_phase_diagonal_from_one','max_source_phase_component')):.3e}; maximal phase-to-z Jacobian column norm {t['phase_sector']['max_phase_to_z_jacobian_column_norm']:.6g}.",
          f"  Median 10D Phi tangent rate finite-difference residual {t['ten_d_phi_variational_fd_relative_median']:.3e}; driven q rate median {t['phase_to_shape_driven_q_fd_relative_median']:.3e}.",
          f"  Median measured Grassmann total/transport/birth speed {t['median_grassmann_total_speed_gauge']:.6g}, {t['median_grassmann_transport_speed_gauge']:.6g}, {t['median_grassmann_birth_speed_gauge']:.6g} (Euclidean gauge).",
          f"  Median projector rate finite difference relative error for xi>=1.1: {t['grassmann_projector_derivative_relative_median_xi_gt_1_1']:.3e}."])
    report.extend(['','## Structural gate','',
     'Original-ODE server test re-integrates the 10D phase-driven variational block independently of the full 11D state and compares it to archived P; source modules are pinned with SHA256. This is a numerical model-consistency audit, not an invariant-manifold or physical-metric certificate.',
     'No independent nonlinear material constraint C(y), D2C, or mechanical G(y) is supplied, so a second-order normal-curvature correction is neither guessed nor claimed. The ODE phase leaf C(s,y)=theta-theta0-2pi*(s-s0) is a kinematic clock identity only.',
     'Every plot has its own declared 10D Euclidean gauge. Invariants: exact linear rank and the invariant phase-fixed variational subbundle relative to the supplied phase split (not arbitrary phase-mixing coordinate transforms).',''])
    (ROOT/'results/REPORT.md').write_text('\n'.join(report)+'\n')
    print(json.dumps({'version':'v18.11','phase_sector':{k:s['phase_sector'] for k,s in idx.items()},'status':'FROZEN_CHECKPOINT_GATE_PASSED'},indent=2))
if __name__=='__main__':main()
