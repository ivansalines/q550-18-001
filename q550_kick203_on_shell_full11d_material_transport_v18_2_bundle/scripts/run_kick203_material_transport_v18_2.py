#!/usr/bin/env python3
"""Q550 kick203 v18.2 — true full11D SUPPORT-STATE material transport audit.

The physically same support-time s is sampled for baseline/+/- trajectories.
The initial kick is identical to v18.1, but it is NOT proven to preserve
constraints unless a verified invariant callback is supplied. Baseline flow
transport is derived from f(y), J(y), with no B finite-difference frames.
The archived v18.1 two-section descriptor is a SEPARATE observable: its
11-vector cannot be subtracted from support-state vectors.
"""
from __future__ import annotations
import argparse, datetime, hashlib, importlib.util, io, json, sys, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
import q550_state_impulse_response_v13 as v13
import history_bridge_v16 as old
from material_geometry import (build_moving_transport, split_increment,
                               frame_from_flow, angle_deg)

VERSION='2026-10-16-Q550-kick203-on-flow-material-transport-v18.2'


def read_inner_npz(zippath, basename):
    with zipfile.ZipFile(zippath) as z:
        found=[name for name in z.namelist() if name.endswith('/'+basename)]
        if len(found)!=1:raise RuntimeError(f'expected exactly one {basename}: {found}')
        with np.load(io.BytesIO(z.read(found[0])),allow_pickle=False) as d:
            return {k:np.asarray(d[k]) for k in d.files}


def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''): h.update(chunk)
    return h.hexdigest()


def constraint_function(path):
    if path is None:return None
    module=v13.load_module(path,'v182_user_constraints')
    function=getattr(module,'invariants',None)
    if not callable(function):raise ValueError('invariant module must export invariants(y) -> vector')
    return lambda y:np.asarray(function(np.asarray(y,float)),float).reshape(-1)


def constraint_distance(fun,y,reference):
    return float(np.max(np.abs(fun(y)-reference)))


def section_parameters(base):
    events=pd.read_csv(base/v13.DEPENDENCIES['events'][0]).set_index('event_n')
    sa0=float(events.loc[390,'support_center'])
    sa1=float(events.loc[433,'support_center'])-sa0
    sb0=float(events.loc[472,'support_center'])
    sb1=float(events.loc[513,'support_center'])-sb0
    if min(sa1,sb1)<=0:raise ValueError('wrong original section orientation')
    return sa0,sa1,sb0,sb1


def get_u(sol0,solplus,solminus,s,s0,amp):
    """Symmetric source-normalized response, evaluated at IDENTICAL support s."""
    baseline=np.asarray(sol0.sol(s),float).reshape(11)
    if s<s0:
        return np.zeros(11),np.zeros(11),baseline
    yp=np.asarray(solplus.sol(s),float).reshape(11)
    ym=np.asarray(solminus.sol(s),float).reshape(11)
    return (yp-ym)/(2*amp),(yp+ym-2*baseline)/(2*amp),baseline


def compute(args):
    args.out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    for p,k in [(args.v181_sendback,'v18_1'),(args.singular_vectors,'singular_vectors')]:
        digest=sha256(p)
        if digest != manifest['sha256'][k]:raise RuntimeError(f'pinned input mismatch: {p}: {digest}')
    provenance=v13.validate_inputs(args.base_root,args.source_bundle_zip)
    prior=read_inner_npz(args.v181_sendback,'second_order_full11D_arrays.npz')
    xi_prior=np.asarray(prior['xi'],float)
    archived_crossing=float(prior['crossing_xi'][0])
    sa0,sa1,sb0,sb1=section_parameters(args.base_root)
    s0=sa0+sa1*args.kick
    crossing=(s0-sb0)/sb1
    if abs(archived_crossing-crossing)>1e-10:raise RuntimeError('v18.1 B crossing differs from source')
    if args.append_exact_crossing:
        xi=np.sort(np.unique(np.r_[xi_prior, crossing]))
    else:xi=xi_prior.copy()
    sA=sa0+sa1*xi; sB=sb0+sb1*xi
    if not np.all(sA>sB):
        raise RuntimeError('A must be chronologically after B for kick203 transport audit')
    stop=float(max(sA.max(),sB.max())+.02)
    observed_support_times=np.unique(np.r_[sA[sA>=s0],sB[sB>=s0],s0])
    if len(observed_support_times)<5:raise RuntimeError('insufficient physical support times')
    invariant=constraint_function(args.invariant_module)
    pa=v13.load_module(args.base_root/v13.DEPENDENCIES['audit'][0],'v182_audit')
    support=v13.load_module(args.base_root/v13.DEPENDENCIES['support'][0],'v182_support')
    tube=v13.load_module(args.base_root/v13.DEPENDENCIES['tube'][0],'v182_tube')
    params=pa.ModelParameters(gamma_rho=.3)
    k0=(55./6.)*(10./550.)
    def rhs_y(y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=k0,params=params,audit=pa),float)
    def jac_y(y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=k0,params=params,audit=pa,base=support),float)
    def rhs(s,y):return rhs_y(y)
    orbit=pa.load_orbit(550.,args.base_root/v13.DEPENDENCIES['orbit'][0])
    tail=pa.collect_uniform_tail(orbit,tail_periods=70,steps_per_period=160,
         samples_per_rotation=720,initial_x_perturbation=1e-6,params=params,progress_every=1000000)
    archived_spline=pa.build_state_spline(tail)
    seed=max(float(tail.s[0]),float(tail.s[-1])-.5)
    if seed>min(min(sA),min(sB),s0):raise RuntimeError('baseline seed is after observed history')
    original_direction=old.columns_for_vectors(args.singular_vectors,args.kick,[1])[1]
    meta=[]; rows=[]; segrows=[]; arrays={'xi':xi,'v18_1_xi':xi_prior,
          'sA':sA,'sB':sB,'kick_support_s':np.array([s0]),
          'xi_B_crossing':np.array([crossing])}
    for tol in args.tolerance_sets:
        rtol,atol={'standard':(1e-12,1e-14),'tight':(3e-13,3e-15)}[tol]
        baseline=old.dense_ivp(rhs,seed,stop,archived_spline(seed),rtol,atol)
        initial_y=np.asarray(baseline.sol(s0),float).reshape(11)
        reference_invariant=invariant(initial_y) if invariant is not None else None
        direction=(v13.invariant_project(original_direction,initial_y,invariant)
                   if invariant is not None else original_direction)
        def state(s):return np.asarray(baseline.sol(float(s)),float).reshape(11)
        # Variational transfer uses the ORIGINAL baseline support ODE, not
        # observer-frame derivatives or 11D descriptor allocations.
        def variational(s,phi):return (jac_y(state(s))@phi.reshape(11,11)).reshape(-1)
        phi_solution=old.dense_ivp(variational,s0,stop,np.eye(11).reshape(-1),rtol,atol)
        P={float(s):np.asarray(phi_solution.sol(float(s)),float).reshape(11,11)
           for s in observed_support_times}
        mappings={}; diagnostics={}
        for rank in (1,2):
            for scale,max_step in [('base',args.transport_step),('refined',args.transport_step/2)]:
                mats,summary=build_moving_transport(observed_support_times,state,rhs_y,jac_y,
                       max_step=max_step,rank=rank,curvature_floor=args.curvature_floor)
                mappings[(rank,scale)]=mats
                diagnostics[f'rank{rank}_{scale}']=summary
        geometry_refinement=float(max(np.linalg.norm(mappings[(2,'base')][float(s)]-
                    mappings[(2,'refined')][float(s)],ord='fro') for s in observed_support_times))
        for label, times in [('A', sA), ('B', sB)]:
            arrays[f'{tol}__baseline_{label}_11D']=np.asarray([state(float(s)) for s in times])
            for rank in (1,2):
                for refinement in ('base','refined'):
                    arrays[f'{tol}__R_{label}_rank{rank}_{refinement}_11x11']=np.asarray([
                        mappings[(rank,refinement)][float(s)] if s>=s0 else np.eye(11)
                        for s in times])
            arrays[f'{tol}__Phi_{label}_11x11']=np.asarray([
                P[float(s)] if s>=s0 else np.full((11,11),np.nan) for s in times])
        for eps in args.eps:
            amp=float(eps)*max(1.,np.linalg.norm(initial_y))
            left=initial_y-amp*direction;right=initial_y+amp*direction
            if invariant is not None:
                left=v13.retract_invariant(left,reference_invariant,invariant,rtol=args.constraint_tol)
                right=v13.retract_invariant(right,reference_invariant,invariant,rtol=args.constraint_tol)
                initial_constraint_error=max(constraint_distance(invariant,left,reference_invariant),
                                             constraint_distance(invariant,right,reference_invariant))
                if initial_constraint_error>args.constraint_tol*max(1.,np.max(abs(reference_invariant))):
                    raise RuntimeError('constraint retraction of kick fails')
            else:initial_constraint_error=None
            plus=old.dense_ivp(rhs,s0,stop,right,rtol,atol)
            minus=old.dense_ivp(rhs,s0,stop,left,rtol,atol)
            uA=[];uB=[];wA=[];wB=[];secondA=[];secondB=[]
            constraint_max=0.
            by_branch={'A':[],'B':[]}
            pre_B_count=0
            for j,x in enumerate(xi):
                record={'tolerance':tol,'eps_rel':eps,'xi':float(x),
                        'relative_xi':float(x-crossing),'sA':float(sA[j]),'sB':float(sB[j])}
                for label,s in [('A',sA[j]),('B',sB[j])]:
                    s=float(s); causal=s>=s0
                    if label=='B' and not causal:pre_B_count+=1
                    u,even,_=get_u(baseline,plus,minus,s,s0,amp)
                    w=P[s]@direction if causal else np.zeros(11)
                    if not causal and np.linalg.norm(u)>1e-30:
                        raise RuntimeError('acausal B response before actual restart')
                    E,_=frame_from_flow(state(s),rhs_y,jac_y,args.curvature_floor)
                    tangent=E[:,0]
                    transverse=u-(u@tangent)*tangent
                    q=dict(record,section=label,post_restart=bool(causal),
                        u_norm=float(np.linalg.norm(u)),linear_norm=float(np.linalg.norm(w)),
                        linear_absolute_error=float(np.linalg.norm(u-w)),
                        linear_relative_error=float(np.linalg.norm(u-w)/max(np.linalg.norm(w),1e-14)) if causal else np.nan,
                        even_defect=float(np.linalg.norm(even)),
                        phase_tangent_component=float(u@tangent),
                        transverse_norm=float(np.linalg.norm(transverse)),
                        direction_vs_linear_deg=angle_deg(u,w))
                    for rank in (1,2):
                        R=mappings[(rank,'refined')][s] if causal else np.eye(11)
                        Rcoarse=mappings[(rank,'base')][s] if causal else np.eye(11)
                        U=R.T@u
                        q[f'body_rank{rank}_frozen_signed']=float(U@direction)
                        q[f'body_rank{rank}_rotation_versus_coarse_norm']=float(np.linalg.norm((R-Rcoarse).T@u))
                        q[f'body_rank{rank}_plane_fraction']=float(np.linalg.norm(E[:,:rank].T@u)/max(np.linalg.norm(u),1e-30)) if np.linalg.norm(u)>1e-14 else np.nan
                    by_branch[label].append((s,u))
                    if label=='A':uA.append(u);wA.append(w);secondA.append(even)
                    else:uB.append(u);wB.append(w);secondB.append(even)
                    rows.append(q)
                    if invariant is not None and causal:
                        for traj in (plus,minus):
                            constraint_max=max(constraint_max,constraint_distance(invariant,
                                      np.asarray(traj.sol(s),float).reshape(11),reference_invariant))
                u_a,u_b=uA[-1],uB[-1]
                wb=wB[-1];wa=wA[-1]
                if sB[j]>=s0:
                    phiA=P[float(sA[j])];phiB=P[float(sB[j])]
                    cond=float(np.linalg.cond(phiB))
                    if not np.isfinite(cond) or cond>args.transfer_condition_gate:
                        cross_error=np.nan;cross_status='condition_gate'
                    else:
                        M=np.linalg.solve(phiB.T,phiA.T).T
                        cross_error=float(np.linalg.norm(u_a-M@u_b))
                        cross_status='B_to_A_variational_transport'
                else:
                    cond=np.nan
                    cross_error=float(np.linalg.norm(u_a-wa))
                    cross_status='kick_source_at_s0_not_homogeneous_B_to_A'
                record.update(cross_section_transport_error=cross_error,
                      transfer_condition=cond,transport_status=cross_status,
                      source_predicted_A_error=float(np.linalg.norm(u_a-wa)))
                # The v18.1 dlocal scalar is contextual only; no cross-space
                # 11D dot product or subtraction of u(s) and Delta_dlocal.
                prior_key=f'{tol}__kick{args.kick:g}__eps{eps:g}__A_hquarter__B_hsixteenth__non_straddling2__dlocal_delta'
                prior_delta=prior.get(prior_key)
                if prior_delta is None:raise RuntimeError(f'missing v18.1 reference: {prior_key}')
                match=np.where(np.isclose(xi_prior,x,rtol=0,atol=5e-12))[0]
                record['prior_dlocal_delta_norm']=float(np.linalg.norm(prior_delta[match[0]])) if len(match)==1 and np.isfinite(prior_delta[match[0]]).all() else np.nan
                rows.append(dict(record,section='A_vs_B_causal_audit',post_restart=True))
            prefix=f'{tol}__eps{eps:g}'
            for field,data in [('uA',uA),('uB',uB),('linear_A',wA),('linear_B',wB),
                               ('even_A',secondA),('even_B',secondB)]:
                arrays[prefix+'__'+field]=np.asarray(data,float)
            arrays[prefix+'__direction_initial']=np.asarray(direction,float)
            # Exact 11D discrete subtraction: successive samples of the SAME
            # physical support branch; DO NOT bridge the instantaneous kick.
            max_closure=0.
            for label,pairs in by_branch.items():
                for rank in (1,2):
                    # Store original full11D physical vectors in addition to scalar CSVs.
                    # Pre-kick support samples carry NaN in body coordinates because
                    # the baseline frame transport is anchored AT the kick time.
                    data_vectors=uA if label=='A' else uB
                    support_times=sA if label=='A' else sB
                    arrays[f'{prefix}__body_{label}_rank{rank}_11D']=np.asarray([
                        mappings[(rank,'refined')][float(s)].T@u
                        if s>=s0 else np.full(11,np.nan)
                        for s,u in zip(support_times,data_vectors)])
                    arrays[f'{prefix}__interval_{label}_rank{rank}_xi_end']=np.asarray([
                        (after-(sa0 if label=='A' else sb0))/(sa1 if label=='A' else sb1)
                        for (before,_),(after,_) in zip(pairs[:-1],pairs[1:])
                        if before>=s0 and after>s0 and after>before],float)
                    intrinsic_history=[];frame_history=[]
                    for (before,ub),(after,ua) in zip(pairs[:-1],pairs[1:]):
                        if before<s0 or after<=s0 or after<=before:continue
                        R0=mappings[(rank,'refined')][before];R1=mappings[(rank,'refined')][after]
                        internal,convected,closure=split_increment(ub,ua,R0,R1)
                        intrinsic_history.append(internal);frame_history.append(convected)
                    arrays[f'{prefix}__interval_{label}_rank{rank}_intrinsic_11D']=np.asarray(intrinsic_history,dtype=float).reshape(-1,11)
                    arrays[f'{prefix}__interval_{label}_rank{rank}_convected_11D']=np.asarray(frame_history,dtype=float).reshape(-1,11)
            for label,pairs in by_branch.items():
                for (before,ub),(after,ua) in zip(pairs[:-1],pairs[1:]):
                    if before<s0 or after<=s0 or after<=before:continue
                    for rank in (1,2):
                        R0=mappings[(rank,'refined')][before];R1=mappings[(rank,'refined')][after]
                        intrinsic,frame,err=split_increment(ub,ua,R0,R1)
                        max_closure=max(max_closure,err)
                        norm_i=float(np.linalg.norm(intrinsic));norm_f=float(np.linalg.norm(frame))
                        segrows.append(dict(tolerance=tol,eps_rel=eps,section=label,rank=rank,
                            physical_s_before=before,physical_s_after=after,
                            ds=after-before,
                            xi_before=(before-(sa0 if label=='A' else sb0))/(sa1 if label=='A' else sb1),
                            xi_after=(after-(sa0 if label=='A' else sb0))/(sa1 if label=='A' else sb1),
                            intrinsic_increment_norm=norm_i,frame_increment_norm=norm_f,
                            observed_increment_norm=float(np.linalg.norm(ua-ub)),
                            intrinsic_fraction_noncancelled=norm_i/max(norm_i+norm_f,1e-300),
                            body_direction_change_deg=angle_deg(R0.T@ub,R1.T@ua),
                            identity_closure=err))
            constraint_status='SAMPLED_CONSTRAINT_GATE_PASSED' if invariant is not None and constraint_max<=args.constraint_tol*max(1.,np.max(abs(reference_invariant))) else ('SAMPLED_CONSTRAINT_GATE_FAILED' if invariant is not None else 'NOT_CERTIFIED_OFF_CONSTRAINT')
            meta.append(dict(tolerance=tol,eps_rel=eps,amp_absolute=amp,source_s=s0,
                  kick_vector_rank=1,kick_initial_constraint_error=initial_constraint_error,
                  constraint_drift_max=constraint_max if invariant is not None else None,
                  constraint_status=constraint_status,
                  pre_restart_B_points=pre_B_count,
                  worst_discrete_identity_error=max_closure,
                  geometry_step_refinement_fro=geometry_refinement,
                  geometry_diagnostics=diagnostics))
            print(f'computed {tol} eps={eps:g}  material physical states A/B: {len(xi)}  '
                  f'constraints={constraint_status}  identity={max_closure:.2e}',flush=True)
    df=pd.DataFrame(rows)
    seg=pd.DataFrame(segrows)
    if df.empty or seg.empty:raise RuntimeError('empty audit')
    df.to_csv(args.out/'material_responses.csv',index=False,float_format='%.17g')
    seg.to_csv(args.out/'moving_frame_increment_decomposition.csv',index=False,float_format='%.17g')
    np.savez_compressed(args.out/'full11D_material_arrays.npz',**arrays)
    # Internal camera never substitutes for physical-state records. Vertical
    # marker indicates B's arrival at the material source support s0.
    fig,ax=plt.subplots(figsize=(10,5))
    subset=df[(df.tolerance==args.tolerance_sets[0]) & (df.eps_rel==args.eps[0]) &
              (df.section.isin(['A','B']))]
    for label in ('A','B'):
        p=subset[subset.section==label]
        ax.plot(p.relative_xi,p.u_norm,label=f'support material response at {label}')
    ax.axvline(0,ls=':',c='gray',label='B reaches kick physical support')
    ax.set_xlabel('xi - xi_B');ax.set_ylabel('|| (y+ - y-) / (2 kick amplitude) || in support 11D')
    ax.legend();ax.grid(alpha=.22);fig.tight_layout()
    fig.savefig(args.out/'material_response_A_B.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5))
    plot=seg[(seg.tolerance==args.tolerance_sets[0])&(seg.eps_rel==args.eps[0])&(seg.rank==2)&(seg.section=='B')]
    ax.plot(plot.xi_after-crossing,plot.intrinsic_increment_norm,label='co-moving increment')
    ax.plot(plot.xi_after-crossing,plot.frame_increment_norm,label='reference transport increment')
    ax.axvline(0,ls=':',c='gray')
    ax.set_xlabel('B xi - xi_B (post-restart only)');ax.set_ylabel('11D increment norm')
    ax.legend();ax.grid(alpha=.22);fig.tight_layout()
    fig.savefig(args.out/'B_increment_intrinsic_vs_frame.png',dpi=170);plt.close(fig)
    summary={'version':VERSION,'s0':s0,'xi_B_crossing':crossing,
             'sample_count':len(xi),'previous_sample_count':len(xi_prior),
             'constraint_module':str(args.invariant_module) if args.invariant_module else None,
             'status':'FLOW_ON_SHELL_ONLY_WITHOUT_CONSTRAINT_CERTIFICATE' if invariant is None else
                      ('SAMPLED_CONSTRAINT_GATE_PASSED_NOT_INTERVAL_CERTIFIED' if all(x['constraint_status']=='SAMPLED_CONSTRAINT_GATE_PASSED' for x in meta) else 'SAMPLED_CONSTRAINT_GATE_FAILED'),
             'max_discrete_identity_closure':max(x['worst_discrete_identity_error'] for x in meta),
             'max_geometry_step_refinement_fro':max(x['geometry_step_refinement_fro'] for x in meta),
             'max_causal_B_response_before_source':float(df[(df.section=='B')&(~df.post_restart)].u_norm.max()),
             'invariant_conservation_gate':args.constraint_tol,
             'experiments':meta,
             'source_provenance':provenance,
             'limitations':[
                 'Support 11D and dlocal descriptor 11D are different observables; no raw vector subtraction.',
                 'The two-plane moving frame is a baseline Euclidean gauge, not an established SO(11) body symmetry.',
                 'Identity closure is algebraic and does not by itself prove real torsion or a flip.',
                 'Without an independent invariant callback, impulse points are off-constraint diagnostics. A callback passing at sampled times is not an interval enclosure certificate.',
                 'The B-before-source material signal vanishes by causal restart construction; not evidence of delayed force.',
             ]}
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (args.out/'REPORT.md').write_text(
        '# Q550 kick203 v18.2 — 11D material-state transport audit\n\n'
        f'- B crossing xi={crossing:.15f}, kick support s0={s0:.15f}.\n'
        f'- Status: **{summary["status"]}**.\n'
        '- A/B compare identical physical support-time samples across the +/- histories, not derivative frames.\n'
        '- The instantaneous impulse at kick is treated as a source; homogeneous B->A transfer is used ONLY when B is post-restart.\n'
        '- Baseline moving frames are generated by analytic flow curvature and transported with minimum ambient rotations.\n'
        '- Rank-one tangent frame, rank-two tangent/curvature frame and half-step transport refinements are controls.\n'
        f'- Max exact increment-decomposition closure: {summary["max_discrete_identity_closure"]:.5g}.\n'
        f'- Max fine/coarse frame-Frobenius discrepancy: {summary["max_geometry_step_refinement_fro"]:.5g}.\n'
        '- Descriptor v18.1 values are stored only as event-aligned context, not mixed into state vector arithmetic.\n'
        '- See material_responses.csv, moving_frame_increment_decomposition.csv, full11D_material_arrays.npz and the figures.\n\n'
        '**Interpretive boundary:** neither a co-moving directional change nor a negative signed projection alone certifies torsional reversal.\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('source_provenance','experiments')},indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base-root',type=Path,default=Path('/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle'))
    ap.add_argument('--source-bundle-zip',type=Path,default=Path('/root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip'))
    ap.add_argument('--v181-sendback',type=Path,default=ROOT/'inputs/q550_kick203_second_order_observer_v18_1_sendback.zip')
    ap.add_argument('--singular-vectors',type=Path,default=ROOT/'inputs/singular_direction_vectors.csv')
    ap.add_argument('--out',type=Path,default=ROOT/'results/kick203_v18_2')
    ap.add_argument('--kick',type=float,default=203.)
    ap.add_argument('--eps',nargs='+',type=float,default=[1e-8,3e-8])
    ap.add_argument('--tolerance-sets',nargs='+',choices=['standard','tight'],default=['standard','tight'])
    ap.add_argument('--transport-step',type=float,default=.01,help='Maximum step in physical support s, independently halved for refinement')
    ap.add_argument('--curvature-floor',type=float,default=1e-11)
    ap.add_argument('--transfer-condition-gate',type=float,default=1e10)
    ap.add_argument('--invariant-module',type=Path,default=None,help='Optional independent model-certified invariants(y) plugin')
    ap.add_argument('--constraint-tol',type=float,default=1e-10)
    ap.add_argument('--append-exact-crossing',action=argparse.BooleanOptionalAction,default=True)
    compute(ap.parse_args())
