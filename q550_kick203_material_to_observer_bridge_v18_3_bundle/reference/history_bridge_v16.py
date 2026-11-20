#!/usr/bin/env python3
"""Q550 ×500 v16 — nonlinear history-to-dlocal directional pushforward.

This program DOES NOT identify support coordinates with dlocal. It differentiates
an explicitly reconstructed, TWO-SECTION observation operator by paired
finite differences of original Q550 nonlinear trajectories and cocycles.

Scientific status is OFF_CONSTRAINT_DIAGNOSTIC without a verified invariant
plugin. Probe-frame metrics are descriptors, NOT material body rotations.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, sys, zipfile
from pathlib import Path
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
import q550_state_impulse_response_v13 as v13

VERSION='2026-10-13-q550-support-dlocal-history-bridge-v16'
SOURCE_SHA=v13.SOURCE_SHA


def unit(x):
    x=np.asarray(x,float); n=float(np.linalg.norm(x))
    if n<1e-14 or not np.isfinite(n):raise ValueError('degenerate direction')
    return x/n


def gram_frame(a,b):
    a=unit(a); b=b-a*(a@b)
    return np.column_stack([a,unit(b)])


def best_o2(src,tgt):
    u,_,vt=np.linalg.svd(src.T@tgt)
    q=u@vt
    return q


def descriptor(Ba,Bb,M,sign_ref=None):
    """Exact BODY difference rule of the historical multiloop extender."""
    pulled=gram_frame(*(M@Bb).T)
    aligned=pulled@best_o2(pulled,Ba)
    dv=aligned[:,1]-Ba@(Ba.T@aligned[:,1])
    if sign_ref is not None and dv@sign_ref<0:dv=-dv
    return dv


def span_intervals(kick,centers,sA,sB):
    """Whether A/B sections are causally reached by an impulse at sA(kick)."""
    return np.asarray(sA)>=kick, np.asarray(sB)>=kick


def dense_ivp(fun,lo,hi,y0,rtol=1e-12,atol=1e-14,max_step=1/720):
    if hi<=lo:raise ValueError('bad interval')
    sol=solve_ivp(fun,(float(lo),float(hi)),np.asarray(y0,float),method='DOP853',
                  rtol=rtol,atol=atol,max_step=max_step,dense_output=True)
    if not sol.success or sol.sol is None:raise RuntimeError(sol.message)
    return sol


@dataclass
class History:
    baseline:object
    perturbed:object|None
    kick_s:float|None
    rhs:object
    jac:object

    def state(self,s):
        if self.perturbed is None or self.kick_s is None or s<self.kick_s:
            return np.asarray(self.baseline.sol(float(s)),float).reshape(11)
        return np.asarray(self.perturbed.sol(float(s)),float).reshape(11)

    def tangent(self,s):return unit(self.rhs(float(s),self.state(s)))

    def frame(self,s,h):
        t=self.tangent(s)
        d=(self.tangent(s+h)-self.tangent(s-h))/(2*h)
        return gram_frame(t,d)

    def jac_at(self,s):return self.jac(float(s),self.state(s))

    def _piecewise_transfer(self,slo,shi,squery,rtol,atol):
        """Locally rebased 11x11 fundamental solution; split at kick discontinuity."""
        q=np.asarray(squery,float)
        if np.any(q<slo-1e-9) or np.any(q>shi+1e-9):raise ValueError('evaluation outside interval')
        def tangent_eq(s,p):return (self.jac_at(s)@p.reshape(11,11)).reshape(-1)
        I=np.eye(11).reshape(-1)
        kick=self.kick_s if self.perturbed is not None else None
        if kick is None or not (slo+1e-12<kick<shi-1e-12):
            sol=dense_ivp(tangent_eq,slo,shi,I,rtol,atol)
            return np.asarray(sol.sol(q),float).T.reshape(-1,11,11)
        first=dense_ivp(tangent_eq,slo,kick,I,rtol,atol)
        at=first.sol(kick)
        second=dense_ivp(tangent_eq,kick,shi,at,rtol,atol)
        ans=[]
        for s in q:
            ans.append(np.asarray((first.sol(s) if s<kick else second.sol(s)),float).reshape(11,11))
        return np.asarray(ans)

    def samples(self,xi,sA,sB,h,block_xi=.75,rtol=1e-12,atol=1e-14,max_cond=1e10):
        """Reconstruct true dlocal, rebasing the relative variational cocycle."""
        if not(len(xi)==len(sA)==len(sB)):raise ValueError('grid length mismatch')
        out=np.zeros((len(xi),11));conds=np.zeros(len(xi))
        def work(ids):
            if len(ids)==0:return
            lo=float(min(np.min(sB[ids]),np.min(sA[ids]))); hi=float(max(np.max(sB[ids]),np.max(sA[ids])))
            try:
                alltimes=np.concatenate([sB[ids],sA[ids]])
                phis=self._piecewise_transfer(lo,hi,alltimes,rtol,atol)
                n=len(ids);pb=phis[:n];pa=phis[n:]
                cond=np.array([max(np.linalg.cond(pb[i]),np.linalg.cond(pa[i])) for i in range(n)])
                if not np.isfinite(cond).all() or cond.max()>max_cond:raise np.linalg.LinAlgError('local cocycle condition')
                for j,k in enumerate(ids):
                    M=np.linalg.solve(pb[j].T,pa[j].T).T
                    BA=self.frame(sA[k],h);BB=self.frame(sB[k],h)
                    out[k]=descriptor(BA,BB,M)
                    conds[k]=cond[j]
            except (np.linalg.LinAlgError,FloatingPointError,RuntimeError,ValueError):
                if len(ids)==1:raise
                cut=len(ids)//2
                work(ids[:cut]);work(ids[cut:])
        dx=np.diff(xi)
        if np.any(dx<=0):raise ValueError('non-increasing xi')
        block=max(1,int(np.floor(block_xi/np.min(dx))))
        for i in range(0,len(xi),block):work(np.arange(i,min(i+block,len(xi))))
        return out,conds


def original_reference(source,stencil,centers):
    if v13.hash_file(source)!=SOURCE_SHA:raise RuntimeError('wrong x500 input source')
    with zipfile.ZipFile(source) as z:
        names=[n for n in z.namelist() if n.endswith(f'/inputs/local_y2_driver_{stencil}_extended.csv')]
        if len(names)!=1:raise RuntimeError('missing unique historical dlocal driver')
        wanted={round(float(x),6) for x in centers};rows={}
        with z.open(names[0]) as f:
            for chunk in pd.read_csv(f,usecols=['xi','sA','sB']+[f'dlocal_{j}' for j in range(11)],chunksize=100000):
                part=chunk[chunk.xi.round(6).isin(wanted)]
                for _,r in part.iterrows():
                    xi=round(float(r.xi),6)
                    if xi in rows:raise RuntimeError('duplicate source xi')
                    rows[xi]=r
                if len(rows)==len(wanted):break
    if len(rows)!=len(wanted):raise RuntimeError(f'original dlocal rows missing: {len(wanted)-len(rows)}')
    data=[rows[round(float(x),6)] for x in centers]
    return (np.array([r.sA for r in data],float),np.array([r.sB for r in data],float),
            np.array([[r[f'dlocal_{j}'] for j in range(11)] for r in data],float))


def read_dlocal_frames(path,variant,xi):
    with np.load(path,allow_pickle=False) as d:
        centers=d[variant+'__centers_xi']
        frame=d[variant+'__global_parallel_frames_11d']
        proj=d[variant+'__projector_11d']
    ix=np.searchsorted(centers,xi)
    if np.any(ix>=len(centers)) or np.max(np.abs(centers[ix]-xi))>1e-8:
        raise RuntimeError('v8 posture reference does not cover exactly these xi centers')
    return frame[ix],proj[ix]


def columns_for_vectors(path,kick,ranks):
    df=pd.read_csv(path)
    results={}
    for rank in ranks:
        chunk=df[(np.isclose(df.kick_xi,kick))&(df.metric=='euclidean')&(df['rank']==rank)].sort_values('channel')
        if len(chunk)!=11:raise RuntimeError(f'no 11-vector singular direction rank={rank}, kick={kick}')
        results[rank]=unit(chunk.initial_physical_unit.to_numpy(float))
    return results


def make_plot(df,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for name,group in df[(df.stencil=='hhalf')&(df.eps_rel==1e-7)].groupby(['kick_xi','rank']):
        plt.figure(figsize=(8,4.5));
        plt.plot(group.xi,group.observed_pushforward_norm,marker='.',label='|delta dlocal| / |delta y0|')
        plt.plot(group.xi,group.body_plane_fraction,label='fraction in dlocal rank-2 plane')
        plt.axvline(group.B_reaches_kick_at_xi.iloc[0],color='gray',linestyle=':',label='B catches impulse')
        plt.xlabel('xi');plt.title(f'Q550 v16 support-to-dlocal response: origin {name[0]:g}, rank {name[1]}')
        plt.grid(alpha=.25);plt.legend();plt.tight_layout()
        plt.savefig(out/f'bridge_xi_{name[0]:g}_rank_{name[1]}.png',dpi=150);plt.close()


def run(args):
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    base=args.base_root;source=args.source_bundle_zip
    hashes=v13.validate_inputs(base,source)
    audit=v13.load_module(base/v13.DEPENDENCIES['audit'][0],'v16_audit')
    support=v13.load_module(base/v13.DEPENDENCIES['support'][0],'v16_support')
    tube=v13.load_module(base/v13.DEPENDENCIES['tube'][0],'v16_tube')
    params=audit.ModelParameters(gamma_rho=.3);kappa0=(55./6.)*(10./550.)
    events=pd.read_csv(base/v13.DEPENDENCIES['events'][0]).set_index('event_n')
    s390,s433,s472,s513=[float(events.loc[n,'support_center']) for n in [390,433,472,513]]
    b0=s472;b1=s513-s472;a0=s390;a1=s433-s390
    if not(a1>0 and b1>0):raise RuntimeError('section slopes must be positive')
    xbeg=float(min(args.kicks));xend=float(max(args.kicks)+args.horizon_xi)
    stamps={kick:kick+np.arange(1,round(args.horizon_xi/args.step_xi)+1)*args.step_xi for kick in args.kicks}
    required_xi=np.unique(np.r_[tuple(stamps.values())]);sA=a0+a1*required_xi;sB=b0+b1*required_xi
    refs={}
    for st in args.stencils:refs[st]=original_reference(source,st,required_xi)
    if any(np.max(abs(refs[st][0]-sA))>1e-8 or np.max(abs(refs[st][1]-sB))>1e-8 for st in refs):
        raise RuntimeError('section-map mismatch against historical x500 original dlocal data')
    shared=json.loads((base/v13.DEPENDENCIES['shared'][0]).read_text())
    p_line=float(shared['p_line']);arcs=int(shared['arc_samples'])
    heights={'hhalf':.125*p_line/(arcs-1),'hquarter':.0625*p_line/(arcs-1)}
    def rhs(s,y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=kappa0,params=params,audit=audit),float)
    def jac(s,y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=kappa0,params=params,audit=audit,base=support),float)
    orb=audit.load_orbit(550.,base/v13.DEPENDENCIES['orbit'][0])
    tail=audit.collect_uniform_tail(orb,tail_periods=70,steps_per_period=160,samples_per_rotation=720,initial_x_perturbation=1e-6,params=params,progress_every=1000000)
    arch=audit.build_state_spline(tail)
    seed=max(float(tail.s[0]),float(tail.s[-1])-.5)
    sstop=float(max(sA.max(),sB.max())+.02)
    if seed>=min(sB):raise RuntimeError('tail seed newer than observation sections')
    baseline=dense_ivp(rhs,seed,sstop,arch(seed),args.rtol,args.atol)
    baseline_history=History(baseline,None,None,rhs,jac)
    # Reject missing or mismatched historical baseline before considering perturbed trajectories.
    base_dlocal={};reference_sign={};base_condition={};regression=[]
    for st in args.stencils:
        original=refs[st][2]
        reconstructed,cond=baseline_history.samples(required_xi,sA,sB,heights[st],args.block_xi,args.rtol,args.atol)
        # Historical sign was transported on a *much denser* xi driver; at our
        # sparse centers the Z2 sheet may differ. Anchor each sample to the
        # verified original driver rather than inventing sign continuity.
        overlap=np.sum(original*reconstructed,axis=1)
        denom=np.linalg.norm(original,axis=1)*np.linalg.norm(reconstructed,axis=1)
        cosine=overlap/np.maximum(denom,1e-300)
        if np.any(np.abs(cosine)<.9):raise RuntimeError('descriptor orientation not identifiable against original source')
        sign=np.where(overlap>=0,1.,-1.)
        reconstructed*=sign[:,None]
        errors=np.linalg.norm(reconstructed-original,axis=1)
        norm=np.sqrt(np.mean(np.sum(original**2,axis=1)))
        rel=float(np.sqrt(np.mean(errors**2))/max(norm,1e-15))
        corr=float(np.corrcoef(reconstructed.ravel(),original.ravel())[0,1])
        regression.append(dict(stencil=st,relative_rms_error=rel,flattened_correlation=corr,max_local_condition=float(cond.max()),reference_per_sample_sign_flips=int(np.sum(sign<0)),minimum_abs_gauge_cosine=float(np.min(abs(cosine)))))
        if rel>args.baseline_rms_tolerance or corr<.999:
            raise RuntimeError(f'baseline dlocal reconstruction FAILS {st}, relative RMS={rel}, correlation={corr}; do NOT interpret pushed directions')
        base_dlocal[st]=reconstructed;reference_sign[st]=sign;base_condition[st]=cond
    pd.DataFrame(regression).to_csv(out/'baseline_regression.csv',index=False)
    rows=[];array_save={};orth_errors=[];causal=[]
    for kick in args.kicks:
        xi=stamps[kick];jidx=np.searchsorted(required_xi,xi)
        ssA=sA[jidx];ssB=sB[jidx]
        s0=float(a0+a1*kick);y0=np.asarray(baseline.sol(s0),float).reshape(11)
        tangent0=unit(rhs(s0,y0))
        directions=columns_for_vectors(args.singular_vectors,float(kick),args.ranks)
        directions[0]=tangent0  # control; rank 0 labels the phase tangent.
        bcatch=(s0-b0)/b1
        causal.append(dict(kick_xi=kick,sA_at_kick=s0,B_reaches_kick_at_xi=bcatch,delay_xi=bcatch-kick,
            check_B_before=(b0+b1*kick)<s0))
        for rank,vec in directions.items():
            for eps in (args.eps_rel if rank==1 else args.eps_rel[:1]):
                amp=float(eps)*max(1.,np.linalg.norm(y0))
                variants=[]
                for sign in (+1,-1):
                    pert=dense_ivp(rhs,s0,sstop,y0+sign*amp*vec,args.rtol,args.atol)
                    history=History(baseline,pert,s0,rhs,jac)
                    # Impulse crosses the A/B sections at different xi; retain history, no 11x11 coordinate identification.
                    by_st={}
                    for st in args.stencils:
                        raw,cond=history.samples(xi,ssA,ssB,heights[st],args.block_xi,args.rtol,args.atol)
                        # Align each perturbed frame with the baseline frame, then restore global historical sign.
                        source_base=base_dlocal[st][jidx]
                        overlaps=np.sum(raw*source_base,axis=1)
                        gauge_cos=overlaps/(np.maximum(np.linalg.norm(raw,axis=1)*np.linalg.norm(source_base,axis=1),1e-300))
                        raw=raw*np.where(overlaps>=0,1.,-1.)[:,None]
                        by_st[st]=(raw,gauge_cos)
                    variants.append(by_st)
                for st in args.stencils:
                    rplus,gauge_plus=variants[0][st];rminus,gauge_minus=variants[1][st];raw_delta=(rplus-rminus)/(2*amp)
                    # Odd FD is derivative per original-state unit, before any support-to-dlocal metric scale.
                    fixed_G,proj=read_dlocal_frames(args.v8_covariances,f'hann_matched__{st}',xi)
                    locked=np.einsum('tji,tj->ti',fixed_G,raw_delta)
                    # These G frames are fixed from unperturbed baseline; do not refit under the impulse.
                    proj2=np.einsum('tij,tj->ti',proj,raw_delta)
                    norms=np.linalg.norm(raw_delta,axis=1);active=np.linalg.norm(proj2,axis=1)
                    fractions=(active/np.maximum(norms,1e-200))**2
                    diff2=np.linalg.norm((rplus+rminus)*.5-base_dlocal[st][jidx],axis=1)
                    # Critical guard: do not interpret sample where finite-difference body stencil crosses discontinuous impulse.
                    min_near=np.minimum(np.minimum(abs(ssA-s0),abs(ssB-s0)),1e3)
                    near=(min_near<=2*heights[st]);
                    for i,x in enumerate(xi):
                        rows.append(dict(kick_xi=kick,rank=int(rank),eps_rel=float(eps),stencil=st,xi=float(x),delta_xi=float(x-kick),
                             section_A_perturbed=bool(ssA[i]>=s0),section_B_perturbed=bool(ssB[i]>=s0),causal_phase=('both_sections' if ssB[i]>=s0 else 'A_only'),B_reaches_kick_at_xi=bcatch,
                             impulse_stencil_nearby=bool(near[i]),
                             gauge_plus_cos_abs=float(abs(gauge_plus[i])),gauge_minus_cos_abs=float(abs(gauge_minus[i])),
                             gauge_ambiguous=bool(min(abs(gauge_plus[i]),abs(gauge_minus[i]))<0.8),
                             observed_pushforward_norm=float(norms[i]),
                             body_plane_fraction=float(fractions[i]),body_plane_gain=float(active[i]),
                             body_normal_gain=float(np.sqrt(max(0,norms[i]**2-active[i]**2))),
                             even_nonlinear_defect_per_amplitude=float(diff2[i]/amp),
                             dlocal_locked_c0=float(locked[i,0]),dlocal_locked_c1=float(locked[i,1]),
                             dlocal_locked_normal_norm=float(np.linalg.norm(locked[i,2:])),
                             source='paired-original-ODE-history-reconstruction',status='OFF_CONSTRAINT_DIAGNOSTIC'))
                    key=f'kick_{str(kick).replace(".","p")}_rank_{rank}_eps_{str(eps).replace(".","p")}_{st}'
                    array_save[key+'_xi']=xi;array_save[key+'_dlocal_derivative']=raw_delta
                    array_save[key+'_dlocal_derivative_body_locked']=locked
                    if st=='hhalf':print(f'kick {kick:g} rank {rank} eps {eps:g}: max |d dlocal/d y0|={norms.max():.4g}',flush=True)
    result=pd.DataFrame(rows);result.to_csv(out/'support_to_dlocal_directional_pushforward.csv',index=False)
    pd.DataFrame(causal).to_csv(out/'section_causality_timeline.csv',index=False)
    # Compare eps pairs where available: linearity gate of derivative estimate.
    linearity=[]
    for keys,sub in result.groupby(['kick_xi','rank','stencil','xi']):
        if len(sub)<2:continue
        a=sub.sort_values('eps_rel').observed_pushforward_norm.to_numpy(float)
        linearity.append(dict(kick_xi=keys[0],rank=keys[1],stencil=keys[2],xi=keys[3],
                              ratio_bigger_to_smaller=float(a[-1]/a[0]) if a[0]>0 else np.nan))
    pd.DataFrame(linearity).to_csv(out/'amplitude_linearity.csv',index=False)
    linearity_v=[]
    for k in args.kicks:
        for st in args.stencils:
            x=stamps[k]
            trial=[]
            for e in args.eps_rel:
                key=f'kick_{str(k).replace(".","p")}_rank_1_eps_{str(e).replace(".","p")}_{st}_dlocal_derivative'
                if key in array_save:trial.append(array_save[key])
            if len(trial)>1:
                mask=(b0+b1*x >= a0+a1*k+2*heights[st])
                aa=np.stack([v[mask] for v in trial])
                relative=float(np.linalg.norm(aa[-1]-aa[0])/max(np.linalg.norm(aa[0]),1e-15)) if aa[0].size else float('nan')
                linearity_v.append(dict(kick_xi=k,stencil=st,relative_difference_of_rank1_finite_differences=relative,n_samples_after_B=int(mask.sum())))
    pd.DataFrame(linearity_v).to_csv(out/'directional_linearity_gate.csv',index=False)
    np.savez_compressed(out/'directional_observation_derivatives.npz',**array_save)
    make_plot(result,out)
    report=dict(version=VERSION,status='OFF_CONSTRAINT_DIAGNOSTIC',
        gauge_ambiguous_rows=int(result.gauge_ambiguous.sum()),
        rank1_linearity_max_relative_difference=max([x['relative_difference_of_rank1_finite_differences'] for x in linearity_v if np.isfinite(x['relative_difference_of_rank1_finite_differences'])],default=None),
        provenance=hashes,origins=args.kicks,ranks=args.ranks,eps_rel=args.eps_rel,
        step_xi=args.step_xi,horizon_xi=args.horizon_xi,stencils=args.stencils,baseline_regression=regression,
        x500_source_sha=v13.hash_file(source),section_causality=causal,linearity_gate=linearity_v,rows=int(len(rows)),
        interpretation=['Nonlinear finite differences through two-section BODY/variational operator; NOT direct identification of both 11D coordinate spaces.',
        'Perturbations applied at support section A: body descriptor section B is affected only after its own later crossing.',
        'Covariance posture lens G comes from independent unperturbed v8; it is not refitted after kicks.',
        'Directional gains depend on support and dlocal coordinate metrics; no physical invariants certified.',
        'No material Dzhanibekov inference without inertia and constraint-preserving trajectory comparison.'])
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['version','status','rows','baseline_regression','section_causality']},indent=2))


def parse():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--base-root',type=Path,default=v13.DEFAULT_BASE)
    a.add_argument('--source-bundle-zip',type=Path,default=v13.DEFAULT_SOURCE)
    a.add_argument('--singular-vectors',type=Path,default=ROOT/'inputs/singular_direction_vectors.csv')
    a.add_argument('--v8-covariances',type=Path,default=ROOT/'inputs/body_locked_frames_and_covariances_v8.npz')
    a.add_argument('--kicks',type=float,nargs='+',default=[203.,219.])
    a.add_argument('--ranks',type=int,nargs='+',default=[1,2,3,10])
    a.add_argument('--eps-rel',type=float,nargs='+',default=[1e-7,3e-7],help='Both levels on rank 1; first level on all other ranks')
    a.add_argument('--horizon-xi',type=float,default=16.)
    a.add_argument('--step-xi',type=float,default=.5)
    a.add_argument('--block-xi',type=float,default=.75)
    a.add_argument('--stencils',nargs='+',choices=['hhalf','hquarter'],default=['hhalf','hquarter'])
    a.add_argument('--baseline-rms-tolerance',type=float,default=1e-5)
    a.add_argument('--rtol',type=float,default=1e-12)
    a.add_argument('--atol',type=float,default=1e-14)
    a.add_argument('--out',type=Path,default=ROOT/'results/server_x500')
    args=a.parse_args()
    if args.horizon_xi<=0 or args.step_xi<=0 or args.block_xi<=0: a.error('positive horizon/step/block required')
    if any((x < 190 or x>225.5) for x in args.kicks):a.error('origins must lie within v14 reference range')
    if any(r<1 or r>10 for r in args.ranks):a.error('singular ranks are 1..10')
    if any(e<=0 for e in args.eps_rel):a.error('eps must be positive')
    if any(abs(args.horizon_xi/args.step_xi-round(args.horizon_xi/args.step_xi))>1e-8 for _ in [0]):a.error('horizon/step not integer')
    if max(args.kicks)+args.horizon_xi>240.:a.error('v8 reference frame ends at xi=240')
    return args

if __name__=='__main__':run(parse())
