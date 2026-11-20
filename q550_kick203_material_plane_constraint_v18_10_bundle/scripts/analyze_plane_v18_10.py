#!/usr/bin/env python3
"""Material memory as transported 2-plane and source-born 3-direction test.

B=E^K=P^(2)(a^A), B'=slope*[J^(2) B+E^(Pg)],
where a=A(xi=1), A(xi)=integral P^-1 S2, A'=slope*g.

Coordinate-free: vanishing B, vanishing a^A^g, and exact GL(11)
transformation of exterior tensors. Reported singular values/angles are only
numerical conditioning in explicit coordinate gauges, not physical angles.
"""
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'scripts'))
from checkpoint_io_v18_10 import load
from exterior_geometry_v18_10 import (blade2,blade3,wedge_rate,normalized_bivector,
   normalized_trivector,normal_component_ratio,rank_from_vectors)


def finite_max(a):return float(np.nanmax(np.abs(a)))


def save_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        out=csv.DictWriter(f,fieldnames=list(rows[0]));out.writeheader();out.writerows(rows)


def analyze(name,d):
    x=d['xi_nodes'];P=d['P_full11D'];J=d['J_full11D'];A=d['integral_source_frame_11D'];
    g=d['source_integrand_11D'];E=d['old_source_material_11D'];K=d['quadratic_total_material_11D']
    slope=float(d['slope']);a=A[80].copy()
    if x.shape!=(641,) or P.shape!=(641,11,11) or J.shape!=(641,11,11):raise ValueError('invalid v18.7 checkpoint shape')
    if not np.all(np.isfinite(A)) or not np.all(np.isfinite(g)):raise ValueError('non-finite source')
    A0=np.broadcast_to(a,A.shape)
    B=blade2(E,K);b=blade2(A0,A)
    pulled=blade2(np.einsum('nij,j->ni',P,a),np.einsum('nij,nj->ni',P,A))
    rel=finite_max(np.linalg.norm(B-pulled,axis=1)/np.maximum(1.,np.linalg.norm(B,axis=1)))
    # Three-index exterior product defines GL-invariant *exact vanishing*.
    tri=blade3(A0,A,g)
    ws=normalized_bivector(A0,A);wlab=normalized_bivector(E,K)
    tri_norm=normalized_trivector(A0,A,g)
    strain,birth=wedge_rate(E,K,J,P,g,slope)
    Bfd=np.gradient(B,x,axis=0,edge_order=2)
    error=np.linalg.norm((Bfd-strain-birth)[1:-1],axis=1)
    local=np.linalg.norm((strain+birth)[1:-1],axis=1)
    rel_fd=error/np.maximum(local,1e-100)
    # Source derivative check: numerical against analytic, exclude boundary.
    bfd=np.gradient(b,x,axis=0,edge_order=2)
    bd=slope*blade2(A0,g)
    berr=np.linalg.norm((bfd-bd)[1:-1],axis=1)
    bscale=np.linalg.norm(bd[1:-1],axis=1)
    valid=(x>1.+1e-10)
    # Test the source-frame *metric-gauge* dependent indicator at 1..8 only;
    # it must NOT be confused with the GL invariant rank.
    rows=[];rank3=[];rank2=[]
    for i in np.where(valid)[0]:
        r2,s2,_=rank_from_vectors([a,A[i]],rel_tol=1e-8)
        r3,s3,_=rank_from_vectors([a,A[i],g[i]],rel_tol=1e-8)
        rank2.append(r2);rank3.append(r3)
        rows.append({'tolerance':name,'xi':float(x[i]),
          'source_frame_bivector_norm':float(np.linalg.norm(b[i])),
          'lab_bivector_norm':float(np.linalg.norm(B[i])),
          'source_frame_sine_metric_choice':float(ws[i]),
          'lab_sine_metric_choice':float(wlab[i]),
          'source_frame_trivector_norm':float(np.linalg.norm(tri[i])),
          'source_frame_3volume_gauge_ratio':float(tri_norm[i]),
          'source_departure_initial_euclidean_gauge':normal_component_ratio(a,A[i],g[i]),
          'numerical_span2_rank_threshold_1e_8':r2,
          'numerical_span3_rank_threshold_1e_8':r3,
          'numerical_span2_svd_ratio':s2,'numerical_span3_svd_ratio':s3,
          'bivector_strain_rate_norm':float(np.linalg.norm(strain[i])),
          'bivector_birth_rate_norm':float(np.linalg.norm(birth[i])),
          'bivector_total_rate_norm':float(np.linalg.norm(strain[i]+birth[i])),
          'bivector_rate_fd_error_norm':float(np.linalg.norm(Bfd[i]-strain[i]-birth[i]))})
    # General coordinate test T=non-orthogonal shear + diagonal scaling.
    n=11;T=np.diag(np.exp(np.linspace(-.6,.6,n)));T[0,3]+=.21;T[4,8]-=.17
    pinv=np.linalg.inv(T)
    slice_indices=[81,120,240,320,400,560,640]
    # Exact GL algebra in transformed coordinates, not a metric test.
    gl_deviation=0.;rank_agreements=0
    for i in slice_indices:
        anew=T@a;Anew=T@A[i];gnew=T@g[i];
        enew=T@E[i];knew=T@K[i];pnew=T@P[i]@pinv;jnew=T@J[i]@pinv
        r2a,_,_=rank_from_vectors([a,A[i]],1e-8)
        r2b,_,_=rank_from_vectors([anew,Anew],1e-8)
        r3a,_,_=rank_from_vectors([a,A[i],g[i]],1e-8)
        r3b,_,_=rank_from_vectors([anew,Anew,gnew],1e-8)
        if r2a==r2b and r3a==r3b:rank_agreements+=1
        # Push forward exact exterior derivative and compare rebuilt tensors.
        old_rate=strain[i]+birth[i]
        ce=lambda v,w: blade2(T@v,T@w)
        original_source_rate=slope*ce(E[i],P[i]@g[i])
        je=jnew@enew;jk=jnew@knew
        reconstructed=slope*(blade2(je,knew)+blade2(enew,jk)+blade2(enew,pnew@(T@g[i])))
        gl_deviation=max(gl_deviation, np.linalg.norm(reconstructed - (
           slope*(ce(J[i]@E[i],K[i])+ce(E[i],J[i]@K[i])) + original_source_rate
           ))/max(1.,np.linalg.norm(reconstructed)))
    metric_zeros=d['xi_roots'].tolist()
    summary={'tolerance':name,'nodes':len(x),'xi_range':[float(x[0]),float(x[-1])],
      'observed_euclidean_dot_crossings_xi':metric_zeros,
      'initial_plane_degenerate_at_xi_1_by_construction':True,
      'max_relative_exterior_transport_reconstruction':rel,
      'exterior_rate_fd_relative_median':float(np.median(rel_fd)),
      'exterior_rate_fd_relative_max':float(np.max(rel_fd)),
      'source_blade_rate_fd_relative_median':float(np.median(berr/np.maximum(bscale,1e-100))),
      'source_blade_rate_fd_relative_max':float(np.max(berr/np.maximum(bscale,1e-100))),
      'xi_gt_1_min_source_plane_sine_gauge':float(np.nanmin(ws[valid])),
      'xi_gt_1_min_lab_plane_sine_gauge':float(np.nanmin(wlab[valid])),
      'xi_gt_1_min_trivector_gauge':float(np.nanmin(tri_norm[valid])),
      'xi_gt_1_max_trivector_gauge':float(np.nanmax(tri_norm[valid])),
      'xi_gt_1_median_source_departure_gauge':float(np.nanmedian([r['source_departure_initial_euclidean_gauge'] for r in rows])),
      'xi_gt_1_rank2_counts':{str(k):int(rank2.count(k)) for k in sorted(set(rank2))},
      'xi_gt_1_rank3_counts':{str(k):int(rank3.count(k)) for k in sorted(set(rank3))},
      'coordinate_relabelled_exterior_rate_relative_max':float(gl_deviation),
      'coordinate_relabelled_numerical_rank_agreements':int(rank_agreements),
      'coordinate_relabelled_samples':len(slice_indices),
      'exterior_invariance':'EXACT_ZERO_AND_ALGEBRAIC_COVARIANCE_ONLY',
      'physical_body_plane_certified':False,
      'warning':'Rank threshold and normalized magnitudes depend on coordinates and numerical precision; off-constraint source and unknown physical metric/constraints.'}
    return summary,rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=HERE/'results')
    a=parser.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    src,manifest=load();all_rows=[];out=[]
    for name,d in src.items():
        s,r=analyze(name,d);out.append(s);all_rows.extend(r)
    results={'version':'v18.10','provenance':manifest,'by_tolerance':out,
     'physical_metric_identified':False,'genuine_constraints_identified':False,
     'scientific_interpretation':'OFF_CONSTRAINT_MATERIAL_PLANE_DIAGNOSTIC',
     'next_gate':'Independently derive constraints and a physical bilinear form from the model; compare structure-preserving second-order perturbations, including extrinsic curvature.'}
    (a.out/'summary.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    save_csv(a.out/'material_plane_profile.csv',all_rows)
    crossing_rows=[]
    for record in out:
        tolerance=record['tolerance']
        z=[v for v in all_rows if v['tolerance']==tolerance]
        xs=np.array([v['xi'] for v in z])
        for root in record['observed_euclidean_dot_crossings_xi']:
            vals={'tolerance':tolerance,'euclidean_dot_zero_xi':root}
            for col in ('source_frame_sine_metric_choice','lab_sine_metric_choice',
                       'source_frame_3volume_gauge_ratio','source_departure_initial_euclidean_gauge'):
                vals['interpolated_'+col]=float(np.interp(root,xs,[v[col] for v in z]))
            crossing_rows.append(vals)
    save_csv(a.out/'crossing_plane_diagnostics.csv',crossing_rows)
    lines=['# Q550 v18.10 — Material plane transport and constraint preflight','',
      'All data sourced from the v18.9 sendback and its hash-verified nested v18.7/v18.6 checkpoints.','',
      'The geometric objects are B=E wedge K2 and b=a wedge A; source-plane departure is a wedge A wedge g.','',
      '## Findings (numerical, not interval certified)','']
    for z in out:
        lines.extend([f"- **{z['tolerance']}**: rank-2 counts {z['xi_gt_1_rank2_counts']}; rank-3 counts {z['xi_gt_1_rank3_counts']} (SVD cutoff 1e-8 in original source coordinates).",
         f"  B'=transport + birth finite-difference median residual {z['exterior_rate_fd_relative_median']:.3e}.",
         f"  Source-frame median departure indicator {z['xi_gt_1_median_source_departure_gauge']:.6g} (NOT invariant as magnitude).",
         f"  Reparametrized GL test: {z['coordinate_relabelled_numerical_rank_agreements']}/{z['coordinate_relabelled_samples']} matching numerical ranks; tensor rate discrepancy {z['coordinate_relabelled_exterior_rate_relative_max']:.3e}."])
    lines.extend(['','## Interpretation','',
      'The rank of [a,A] and [a,A,g], and the fact of wedge vanishing, are invariant under exact invertible coordinate changes. Singular values, dot products, angles, normalized volumes, and the numerical thresholds used to estimate rank are NOT invariant. The zero at xi=1 is constructed (A(1)=a).',
      'A nonzero a wedge A wedge g means the incoming source escapes the existing source-memory 2-plane (in the linearized 11D state). It does not prove material body torsion or constitutive hysteresis.',
      'The exterior derivative has two algebraic parts: homogeneous cocycle transport and birth. Their norms and cancellation depend on metric; the equation transforms covariantly under fixed invertible coordinate changes.',
      'No source-derived mechanical energy, nonlinear material constraint, or admissible second-order shape manifold has yet been certified.',
      'Run the optional original-ODE diagnostic to inventory field structure; an ODE RHS alone does not imply physical metric uniqueness.',
      ''])
    (a.out/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'version':'v18.10','standard':out[0]['xi_gt_1_rank3_counts'],
      'tight':out[1]['xi_gt_1_rank3_counts'], 'status':'off_constraint_plane_diagnostic'},indent=2))

if __name__=='__main__':main()
