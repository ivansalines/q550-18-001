#!/usr/bin/env python3
"""Q550 v17.9 event-micro lens, functionally derived from v17.5 exact 11D trajectories.

WARNING: OFF_CONSTRAINT finite-amplitude diagnostic. Canonical tangent-curvature
frames are NOT body-fixed coordinates; M and its SO polar factor are mathematical
transport operators, not rigid-body rotation measurements.
"""
from __future__ import annotations
import argparse, hashlib, io, json, sys, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'reference'), str(ROOT/'scripts')]
import history_bridge_v16 as old
import q550_state_impulse_response_v13 as v13
from frame_geometry import frame_comparison, continuous_lift, mark_phase_events, shapley_three_factors
VERSION='2026-10-16-Q550-kick203-event-micro-lens-v17.9-derived-from-v17.5'


def reconstruct(history, xi, sA, sB, h_by_stencil, rtol, atol, block_xi, max_cond=1e10):
    """Recompute M exactly as v16/v17.4, caching it for all stencil heights.

    Adaptive subdivision preserves original max_cond guard. Samples after a
    kick always derive from the appropriate common-restarted nonlinear branch.
    """
    n=len(xi)
    matrices=np.full((n,11,11),np.nan)
    cond=np.full(n,np.nan)
    sizes=np.zeros(n,dtype=int)
    # Nonuniform event-aligned grid: group points by actual xi span, not by
    # the minimum grid gap. This changes computational batching only, never
    # the original cocycle, frame, or Shapley definitions.
    if not np.all(np.diff(xi)>0):
        raise ValueError('xi grid must be strictly increasing')
    def work(idx):
        lo=float(min(np.min(sB[idx]),np.min(sA[idx])))
        hi=float(max(np.max(sB[idx]),np.max(sA[idx])))
        try:
            q=np.concatenate([sB[idx],sA[idx]])
            phi=history._piecewise_transfer(lo,hi,q,rtol,atol)
            m=len(idx);pb=phi[:m];pa=phi[m:]
            cc=np.array([max(np.linalg.cond(pb[j]),np.linalg.cond(pa[j])) for j in range(m)])
            if not np.all(np.isfinite(cc)) or np.max(cc)>max_cond: raise np.linalg.LinAlgError('block conditioning')
            for j,k in enumerate(idx):
                M=np.linalg.solve(pb[j].T,pa[j].T).T
                if not np.isfinite(M).all():raise np.linalg.LinAlgError('nonfinite M')
                matrices[k]=M
                cond[k]=cc[j];sizes[k]=len(idx)
        except (ValueError,RuntimeError,np.linalg.LinAlgError,FloatingPointError):
            if len(idx)==1:raise
            c=len(idx)//2;work(idx[:c]);work(idx[c:])
    i=0
    while i<n:
        end=int(np.searchsorted(xi,xi[i]+block_xi,side='right'))
        end=max(i+1,end)
        work(np.arange(i,min(end,n)))
        i=end
    frames={}
    for st,h in h_by_stencil.items():
        BA=np.stack([history.frame(float(t),h) for t in sA])
        BB=np.stack([history.frame(float(t),h) for t in sB])
        frames[st]=(BA,BB)
    return matrices,cond,sizes,frames


def source_check(args):
    actual=v13.validate_inputs(args.base_root,args.source_bundle_zip)
    if v13.hash_file(args.source_bundle_zip)!=v13.SOURCE_SHA:raise RuntimeError('x500 SHA changed')
    return actual


def make_histories(args, kick, baseline, sA0, sA1, stop, rt, at, direction, rhs, jac):
    s0=float(sA0+sA1*kick)
    y0=np.asarray(baseline.sol(s0),float).reshape(11)
    nscale=max(1.,float(np.linalg.norm(y0)))
    for eps in [0.]+args.eps:
        signs=[0] if eps==0 else [+1,-1]
        for sign in signs:
            sy=old.dense_ivp(rhs,s0,stop,y0+sign*eps*nscale*direction,rt,at)
            yield (eps,sign,eps*nscale),old.History(baseline,sy,s0,rhs,jac)


def serialize(obj):
    if isinstance(obj,dict):return {str(k):serialize(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):return [serialize(x) for x in obj]
    if isinstance(obj,(np.integer,)):return int(obj)
    if isinstance(obj,(np.bool_,)):return bool(obj)
    if isinstance(obj,(np.floating,)):return serialize(float(obj))
    if isinstance(obj,float) and not np.isfinite(obj):return None
    return obj


def run(args):
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    provenance=source_check(args)
    base=args.base_root
    audit=v13.load_module(base/v13.DEPENDENCIES['audit'][0],'v175_audit')
    support=v13.load_module(base/v13.DEPENDENCIES['support'][0],'v175_support')
    tube=v13.load_module(base/v13.DEPENDENCIES['tube'][0],'v175_tube')
    params=audit.ModelParameters(gamma_rho=.3)
    k0=(55./6.)*(10./550.)
    ev=pd.read_csv(base/v13.DEPENDENCIES['events'][0]).set_index('event_n')
    a0=float(ev.loc[390,'support_center']);a1=float(ev.loc[433,'support_center'])-a0
    b0=float(ev.loc[472,'support_center']);b1=float(ev.loc[513,'support_center'])-b0
    if min(a1,b1)<=0:raise ValueError('nonpositive section scales')
    shared=json.loads((base/v13.DEPENDENCIES['shared'][0]).read_text())
    heights={'hhalf':.125*float(shared['p_line'])/(int(shared['arc_samples'])-1),
             'hquarter':.0625*float(shared['p_line'])/(int(shared['arc_samples'])-1)}
    def rhs(s,y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=k0,params=params,audit=audit),float)
    def jac(s,y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=k0,params=params,audit=audit,base=support),float)
    orb=audit.load_orbit(550.,base/v13.DEPENDENCIES['orbit'][0])
    tail=audit.collect_uniform_tail(orb,tail_periods=70,steps_per_period=160,samples_per_rotation=720,
                                    initial_x_perturbation=1e-6,params=params,progress_every=1000000)
    arch=audit.build_state_spline(tail)
    seed=max(float(tail.s[0]),float(tail.s[-1])-.5)
    b_cross={kick:float((a0+a1*kick-b0)/b1) for kick in args.kicks}
    # Preserve the historical xi lattice for projector regression, add a
    # small original 0.025-lattice overlap, and place a dense event-relative
    # window around the exact B crossing. No physical interpolation is used.
    probes=np.asarray([-0.025,-0.01,-0.005,-0.001,-0.00025,-0.00005,-0.00001,
                       0.0,0.00001,0.00005,0.00025,0.001,0.005,0.01,0.025],float)
    grids={}
    for k in args.kicks:
        xb=b_cross[k]
        historical=np.arange(k-args.pre_kick,xb+args.post_B+1e-9,args.step)
        old_fine=np.arange(args.overlap_xi_min,args.overlap_xi_max+1e-9,.025)
        micro=xb+np.arange(-args.micro_pre,args.micro_post+0.5*args.micro_step,args.micro_step)
        special=xb+probes
        candidates=np.sort(np.concatenate((historical,old_fine,micro,special)))
        dedup=[]
        for x in candidates:
            if dedup and abs(float(x)-dedup[-1]) < 1e-9:
                # Retain the exact crossing node when floating-point siblings
                # from independent grids differ only by roundoff.
                if abs(float(x)-xb)<1e-12: dedup[-1]=float(x)
            else:
                dedup.append(float(x))
        grids[k]=np.array(dedup,float)
        if np.min(grids[k])<k-args.pre_kick-1e-9:
            raise ValueError('overlap extends before historical grid')
        print('event grid kick',k,'crossing',xb,'n=',len(grids[k]),
              'nearest_left=',float(grids[k][grids[k]<xb][-1]),
              'nearest_right=',float(grids[k][grids[k]>xb][0]),flush=True)
    if any(arr[-1]>250 or arr[0]<0 for arr in grids.values()):raise ValueError('outside frozen x250 domain')
    endxi=max(v.max() for v in grids.values())
    stop=float(max(a0+a1*endxi,b0+b1*endxi)+max(.1,2*max(heights.values())))
    if seed>=min(b0+b1*np.min(v) for v in grids.values()):raise ValueError('seed exceeds section B')
    print('source validated, grids', {str(k):[len(v),v[0],v[-1]] for k,v in grids.items()},flush=True)
    cache={};rows=[];events=[];branches=[];regression=[]
    for tol in args.tolerance_sets:
        rt,at={'standard':(1e-12,1e-14),'tight':(3e-13,3e-15)}[tol]
        baseline=old.dense_ivp(rhs,seed,stop,arch(seed),rt,at)
        for kick in args.kicks:
            xi=grids[kick]; sa=a0+a1*xi;sb=b0+b1*xi;s0=float(a0+a1*kick)
            direction=old.columns_for_vectors(args.singular_vectors,kick,[1])[1]
            for (eps,sign,absolute),history in make_histories(args,kick,baseline,a0,a1,stop,rt,at,direction,rhs,jac):
                M,cond,sizes,frames=reconstruct(history,xi,sa,sb,heights,rt,at,args.block)
                branch=f'{tol}__kick{kick:g}__eps{eps:g}__sign{sign:+d}'
                cache[branch+'__xi']=xi
                cache[branch+'__M']=M
                cache[branch+'__local_cond']=cond
                for st,(BA,BB) in frames.items():
                    subkey=branch+'__'+st
                    cache[subkey+'__BA']=BA
                    cache[subkey+'__BB']=BB
                    results=[];dv=[]
                    for i in range(len(xi)):
                        vals,dlocal,_=frame_comparison(BA[i],BB[i],M[i],args.min_overlap)
                        vals.update(xi=float(xi[i]),B_relative=float(xi[i]-b_cross[kick]),kick_xi=float(kick),
                            B_crossing_xi=b_cross[kick],stencil=st,tolerance=tol,
                            eps_rel=eps,eps_absolute=absolute,sign=int(sign),block_size=int(sizes[i]),
                            local_transfer_condition=float(cond[i]),
                            near_kick_stencil=bool(min(abs(sa[i]-s0),abs(sb[i]-s0))<=2*heights[st]))
                        results.append(vals);dv.append(dlocal)
                    cache[subkey+'__dlocal']=np.asarray(dv)
                    # Pointwise frame phases are computed before deciding event candidates.
                    for kind in ('full','so'):
                        phase=np.array([r['phase_'+kind+'_rad'] for r in results])
                        pr=np.array([r['parity_'+kind] for r in results])
                        ov=np.array([r['min_overlap_sv_'+kind] for r in results])
                        lift,rate,seg,alias=continuous_lift(phase,xi,pr,ov,args.min_overlap,args.max_unaliased_step_deg)
                        for i,r in enumerate(results):
                            r[f'lift_{kind}_deg']=float(np.degrees(lift[i]))
                            r[f'rate_{kind}_deg_per_xi']=float(rate[i]);r[f'segment_{kind}']=int(seg[i])
                            r[f'alias_risk_{kind}']=bool(alias[i])
                        if eps==0:
                            for i,g,v,width in mark_phase_events(xi,lift,seg,args.candidate_step_deg,
                                     args.candidate_window_xi,args.candidate_window_deg):
                                events.append(dict(kick_xi=kick,stencil=st,tolerance=tol,kind=kind,
                                   xi=float(xi[i]),B_relative=float(xi[i]-b_cross[kick]),
                                   criterion=g,turn_deg=v,width_xi=width,
                                   min_overlap_sv=float(ov[i]),source_status='FRAME_PHASE_PROXY_ONLY'))
                    rows.extend(results)
                    # Exact saved v17.4 projection regression on overlapping
                    # original half-step grid, not an interpolation.
                    if eps==0:
                        with zipfile.ZipFile(args.v17_4_sendback) as z:
                            with np.load(io.BytesIO(z.read('server_x500/common_restart_arrays.npz')),allow_pickle=False) as previous:
                                key=f'{tol}__kick{kick:g}__{st}'
                                oldxi=previous[key+'__xi']
                                indexes=[(i,int(np.argmin(abs(oldxi-x)))) for i,x in enumerate(xi) if np.min(abs(oldxi-x))<1e-9]
                                if len(indexes)<8:raise RuntimeError('insufficient exact archive overlap')
                                ia=np.array([i for i,j in indexes]);ib=np.array([j for i,j in indexes]);
                                for label,B in [('A',BA),('B',BB)]:
                                    newproj=B[ia]@B[ia].transpose(0,2,1)
                                    oldproj=previous[key+f'__frame{label}_zero'][ib]
                                    rms=float(np.linalg.norm(newproj-oldproj)/max(np.linalg.norm(oldproj),1e-300))
                                    regression.append(dict(kick_xi=kick,tolerance=tol,stencil=st,
                                     section=label,n_exact_samples=len(ia),relative_projector_error=rms))
                                    if rms>args.regression_gate:
                                        raise RuntimeError(f'original frame projector mismatch {key} {label} {rms:g}')
                branches.append(dict(kick_xi=kick,tolerance=tol,eps_rel=eps,sign=sign,n=len(xi)))
                print('reconstructed',branch,'n=',len(xi),flush=True)
    # v17.5 addendum: mathematically exact, signed nonlinear observer allocation.
    # Shapley triplet averages all orderings of substitutions BA,BB,M. This does
    # NOT provide causal component histories because its hybrids are off-shell.
    allocations=[];alloc_arrays={};hybrid_invalid=0
    for tol in args.tolerance_sets:
        for kick in args.kicks:
            xx=grids[kick]
            for er in args.eps:
                ptag=f'{tol}__kick{kick:g}__eps{er:g}__sign+1'
                mtag=f'{tol}__kick{kick:g}__eps{er:g}__sign-1'
                for st in heights:
                    pp=ptag+'__'+st;mm=mtag+'__'+st
                    contributions=np.full((len(xx),3,11),np.nan)
                    full_delta=np.full((len(xx),11),np.nan)
                    for i,x in enumerate(xx):
                        minus=(cache[mm+'__BA'][i],cache[mm+'__BB'][i],cache[mtag+'__M'][i])
                        plus=(cache[pp+'__BA'][i],cache[pp+'__BB'][i],cache[ptag+'__M'][i])
                        try:
                            share,target,abs_error=shapley_three_factors(minus,plus)
                        except (ValueError,np.linalg.LinAlgError,FloatingPointError):
                            hybrid_invalid+=1
                            allocations.append(dict(tolerance=tol,kick_xi=kick,eps_rel=er,stencil=st,xi=x,
                                B_relative=x-b_cross[kick],status='DEGENERATE_HYBRID_UNIDENTIFIED'))
                            continue
                        if abs_error>1e-10:
                            raise RuntimeError(f'Shapley algebraic closure failed: {abs_error}')
                        direct=cache[pp+'__dlocal'][i]-cache[mm+'__dlocal'][i]
                        # Exact source branches must match direct difference.
                        if not np.allclose(target,direct,atol=1e-10,rtol=1e-7):
                            raise RuntimeError('counterfactual endpoints mismatch actual dlocal')
                        contributions[i]=share;full_delta[i]=target
                        norms=np.linalg.norm(share,axis=1)
                        magnitude=float(np.linalg.norm(target))
                        directional=(share@target)/max(magnitude,1e-300)
                        allocations.append(dict(tolerance=tol,kick_xi=kick,eps_rel=er,stencil=st,xi=x,
                           B_relative=x-b_cross[kick],status='EXACT_MATH_ALLOCATION_OFF_SHELL_HYBRIDS',
                           observed_delta_norm=magnitude,component_A_norm=float(norms[0]),
                           component_B_norm=float(norms[1]),component_M_norm=float(norms[2]),
                           signed_A_along_final=float(directional[0]),signed_B_along_final=float(directional[1]),
                           signed_M_along_final=float(directional[2]),
                           directed_sum_vs_final=float(np.sum(directional)-magnitude),
                           norm_compensation_proxy=float(magnitude/max(sum(norms),1e-300)),
                           closure_absolute_error=abs_error))
                    tag=f'{tol}__kick{kick:g}__eps{er:g}__{st}'
                    alloc_arrays[tag+'__shapley_ABC']=contributions
                    alloc_arrays[tag+'__dlocal_delta']=full_delta
    pd.DataFrame(allocations).to_csv(out/'oriented_shapley_allocation.csv',index=False,float_format='%.16g')
    np.savez_compressed(out/'oriented_shapley_arrays.npz',**alloc_arrays)
    data=pd.DataFrame(rows)
    data.to_csv(out/'oriented_flip_timeseries.csv',index=False,float_format='%.16g')
    pd.DataFrame(events).to_csv(out/'oriented_flip_candidates.csv',index=False)
    pd.DataFrame(regression).to_csv(out/'historical_frame_regression.csv',index=False)
    np.savez_compressed(out/'full_frames_cocycle_arrays.npz',**cache)
    # Cross-observer candidate consensus only if BOTH stencil and tol display
    # predeclared candidate near same B-relative xi, with same kind and kick.
    candidates=pd.DataFrame(events)
    consensus=[]
    for kick in args.kicks:
        for kind in ('full','so'):
            sub=candidates[(candidates.kick_xi==kick)&(candidates.kind==kind)] if len(candidates) else pd.DataFrame()
            for pos in sorted(sub.B_relative.unique()) if len(sub) else []:
                near=sub[abs(sub.B_relative-pos)<=args.consensus_halfwidth]
                observers=set(zip(near.tolerance,near.stencil))
                required={(t,s) for t in args.tolerance_sets for s in heights}
                if required.issubset(observers):
                    consensus.append(dict(kick_xi=kick,kind=kind,B_relative=float(pos),
                            observers=len(observers),candidate_count=int(len(near)),
                            status='ROBUST_OBSERVER_PHASE_PROXY_CANDIDATE_NOT_PHYSICAL_FLIP'))
    pd.DataFrame(consensus,columns=['kick_xi','kind','B_relative','observers','candidate_count','status']).to_csv(out/'observer_consensus.csv',index=False)
    report=dict(version=VERSION,source_status='OFF_CONSTRAINT_DIAGNOSTIC',
       question='Where, if anywhere, does oriented A/B observer frame phase show a reproducible rapid turn?',
       interpretation='11D Gram tangent/curvature frame and locally-rebased cocycle phase proxy; never physical material-body flip.',
       source_sha=v13.hash_file(args.source_bundle_zip),provenance=provenance,
       steps=dict(xi_step='nonuniform',historical_xi_step=args.step,micro_step=args.micro_step, micro_pre=args.micro_pre,micro_post=args.micro_post, event_exact_and_symmetric_probes=True,pre_kick=args.pre_kick,post_B=args.post_B,block=args.block,eps=args.eps,
           tolerances=args.tolerance_sets,stencils=list(heights),min_overlap_sv=args.min_overlap,
           max_unaliased_step_deg=args.max_unaliased_step_deg,candidate_single_step_deg=args.candidate_step_deg,
           candidate_window_xi=args.candidate_window_xi,candidate_window_deg=args.candidate_window_deg),
       B_crossings=b_cross,branches=branches,frame_regression=regression,
       n_proxy_candidates=len(events),n_observer_consensus=len(consensus),n_hybrid_degenerate=hybrid_invalid,
       observed_flip_certified=False,
       cautions=['No invariant-preserving kick is used.',
         '11D tangent/curvature basis is a canonical probe construction, not a material body-attached frame.',
         'Signed SO(2) phase is chart/gauge-dependent; parity-negative and poorly overlapping samples are not assigned a phase.',
         'SO polar(M) is a numerical factor of a local variational cocycle; not a physical rigid-body rotation.',
         'Shapley factor substitutions BA/BB/M close algebraically exactly on observed dlocal, but hybrid histories are not physical Q550 trajectories.',
         'Event-relative nonuniform xi grid diagnoses one-sided allocation geometry only, not a physical impulse or material flip.',
         'Rebasing and transfer can depend on block width; repeat with changed block before physics interpretation.'])
    (out/'oriented_flip_report.json').write_text(json.dumps(serialize(report),indent=2,allow_nan=False)+'\n')
    print('Completed',VERSION,'events',len(events),'consensus',len(consensus),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-root',type=Path,default=v13.DEFAULT_BASE)
    p.add_argument('--source-bundle-zip',type=Path,default=v13.DEFAULT_SOURCE)
    p.add_argument('--singular-vectors',type=Path,default=ROOT/'inputs/singular_direction_vectors.csv')
    p.add_argument('--v17-4-sendback',type=Path,default=ROOT/'inputs/v17_4_original_sendback.zip')
    p.add_argument('--out',type=Path,default=ROOT/'results/server_x500')
    p.add_argument('--kicks',nargs='+',type=float,default=[203.,219.])
    p.add_argument('--eps',nargs='+',type=float,default=[1e-8,3e-8])
    p.add_argument('--tolerance-sets',nargs='+',choices=['standard','tight'],default=['standard','tight'])
    p.add_argument('--pre-kick',type=float,default=1.0)
    p.add_argument('--post-B',type=float,default=12.)
    p.add_argument('--step',type=float,default=.25)
    p.add_argument('--block',type=float,default=.75)
    p.add_argument('--micro-step',type=float,default=.0025)
    p.add_argument('--micro-pre',type=float,default=.15)
    p.add_argument('--micro-post',type=float,default=.35)
    p.add_argument('--overlap-xi-min',type=float,default=210.5)
    p.add_argument('--overlap-xi-max',type=float,default=211.25)
    p.add_argument('--min-overlap',type=float,default=.20)
    p.add_argument('--max-unaliased-step-deg',type=float,default=150.)
    p.add_argument('--candidate-step-deg',type=float,default=45.)
    p.add_argument('--candidate-window-xi',type=float,default=2.)
    p.add_argument('--candidate-window-deg',type=float,default=150.)
    p.add_argument('--consensus-halfwidth',type=float,default=.5)
    p.add_argument('--regression-gate',type=float,default=1e-3)
    a=p.parse_args()
    if a.step<=0 or abs(.5/a.step-round(.5/a.step))>1e-10:p.error('step must exactly refine old half-step grid')
    if a.step>.25:p.error('addendum requires step <=0.25xi')
    if a.pre_kick<1 or a.post_B<10 or a.block<a.step:p.error('the frozen test needs >=1xi lead, >=10xi tail, block >= step')
    if not (0<a.micro_step<=.025 and 0<a.micro_pre<1 and 0<a.micro_post<1):p.error('invalid micro window')
    if a.overlap_xi_max<=a.overlap_xi_min:p.error('invalid overlap range')
    if len(set(a.eps))!=len(a.eps) or min(a.eps)<=0:p.error('positive distinct eps')
    run(a)
