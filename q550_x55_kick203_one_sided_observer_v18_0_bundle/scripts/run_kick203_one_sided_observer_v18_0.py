#!/usr/bin/env python3
"""Q550 kick203 v18.0: centered versus unilateral B-section frame.

Reuses the exact v17.9 full11D tangent cocycle M and frozen section A frames.
The nonlinear +/- branches are regenerated from the same original source and
checked against the v17.9 archived central B frames, without recomputing M.
Only the tangent difference used to construct the B frame is changed.

No claim of physical impulsive flip, causal Shapley attribution, or uniquely
signed 11D torsion is made. All comments, reports, and plots are in English.
"""
from __future__ import annotations
import argparse, hashlib, io, json, sys, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
import history_bridge_v16 as old
import q550_state_impulse_response_v13 as v13
from frame_geometry import shapley_three_factors

VERSION='2026-10-16-Q550-kick203-centered-vs-one-sided-observer-v18.0'


def read_v179_arrays(zpath):
    with zipfile.ZipFile(zpath) as z:
        def load(name):
            hits=[x for x in z.namelist() if x.endswith('/'+name)]
            if len(hits)!=1:raise RuntimeError(f'expected one {name}, found {hits}')
            with np.load(io.BytesIO(z.read(hits[0])),allow_pickle=False) as npz:
                return {k:np.asarray(npz[k]) for k in npz.files}
        return load('full_frames_cocycle_arrays.npz'),load('oriented_shapley_arrays.npz')


def three_frames(history,s,h):
    """Central, forward, backward Gram frames from exactly the same 11D history."""
    t=history.tangent(float(s))
    tm=history.tangent(float(s-h))
    tp=history.tangent(float(s+h))
    return {
        'centered':old.gram_frame(t,(tp-tm)/(2*h)),
        'forward':old.gram_frame(t,(tp-t)/h),
        'backward':old.gram_frame(t,(t-tm)/h),
    }


def angle_deg(u,v):
    a=np.linalg.norm(u);d=np.linalg.norm(v)
    return float(np.degrees(np.arccos(np.clip((u@v)/(a*d),-1,1)))) if a*d>1e-25 else float('nan')


def safe_unit(v):
    s=np.linalg.norm(v)
    return v/s if s>1e-25 else np.zeros_like(v)


def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    from scipy.integrate import solve_ivp  # explicit dependency check
    # Immutable source identity (the same exact source used by v17.9).
    original_hash=v13.hash_file(args.source_bundle_zip)
    if original_hash!=v13.SOURCE_SHA:
        raise RuntimeError('Original x500 source SHA is different: refusing altered genealogy')
    provenance=v13.validate_inputs(args.base_root,args.source_bundle_zip)
    cached,orig_shap=read_v179_arrays(args.v179_sendback)
    audit=v13.load_module(args.base_root/v13.DEPENDENCIES['audit'][0],'v180_audit')
    support=v13.load_module(args.base_root/v13.DEPENDENCIES['support'][0],'v180_support')
    tube=v13.load_module(args.base_root/v13.DEPENDENCIES['tube'][0],'v180_tube')
    params=audit.ModelParameters(gamma_rho=.3)
    k0=(55./6.)*(10./550.)
    ev=pd.read_csv(args.base_root/v13.DEPENDENCIES['events'][0]).set_index('event_n')
    a0=float(ev.loc[390,'support_center']);a1=float(ev.loc[433,'support_center'])-a0
    b0=float(ev.loc[472,'support_center']);b1=float(ev.loc[513,'support_center'])-b0
    if min(a1,b1)<=0:raise RuntimeError('invalid section scales')
    shared=json.loads((args.base_root/v13.DEPENDENCIES['shared'][0]).read_text())
    heights={'hhalf':.125*float(shared['p_line'])/(int(shared['arc_samples'])-1),
             'hquarter':.0625*float(shared['p_line'])/(int(shared['arc_samples'])-1)}
    crossing=(a0+a1*args.kick-b0)/b1
    key0=f'standard__kick{args.kick:g}__eps0__sign+0'
    xi=np.asarray(cached[key0+'__xi'],float)
    n=len(xi)
    idx=np.where((np.abs(xi-crossing)<=args.micro_radius)|
                 np.isclose(xi,211.,atol=1e-9,rtol=0)|
                 np.isclose(xi,211.25,atol=1e-9,rtol=0)|
                 np.isclose(xi,212.10,atol=1e-9,rtol=0))[0]
    if not len(idx):raise RuntimeError('empty probe selection')
    si=a0+a1*xi;sj=b0+b1*xi
    s0=float(a0+a1*args.kick)
    stop=float(max(a0+a1*xi.max(),b0+b1*xi.max())+max(.1,2*max(heights.values())))
    def rhs(s,y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=k0,params=params,audit=audit),float)
    def jac(s,y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=k0,params=params,audit=audit,base=support),float)
    orb=audit.load_orbit(550.,args.base_root/v13.DEPENDENCIES['orbit'][0])
    tail=audit.collect_uniform_tail(orb,tail_periods=70,steps_per_period=160,samples_per_rotation=720,
                                    initial_x_perturbation=1e-6,params=params,progress_every=1000000)
    arch=audit.build_state_spline(tail)
    seed=max(float(tail.s[0]),float(tail.s[-1])-.5)
    if seed>=float(min(sj)):raise RuntimeError('start outside original baseline support')
    direction=old.columns_for_vectors(args.singular_vectors,args.kick,[1])[1]
    rows=[];regression=[];npz_outputs={}
    for tol in args.tolerance_sets:
        rt,at={'standard':(1e-12,1e-14),'tight':(3e-13,3e-15)}[tol]
        baseline=old.dense_ivp(rhs,seed,stop,arch(seed),rt,at)
        y0=np.asarray(baseline.sol(s0),float).reshape(11)
        nscale=max(1.,float(np.linalg.norm(y0)))
        for er in args.eps:
            branches={};original_branches={}
            for sign in (+1,-1):
                sy=old.dense_ivp(rhs,s0,stop,y0+sign*er*nscale*direction,rt,at)
                history=old.History(baseline,sy,s0,rhs,jac)
                tag=f'{tol}__kick{args.kick:g}__eps{er:g}__sign{sign:+d}'
                assert np.array_equal(xi,cached[tag+'__xi'])
                # M is taken unchanged from the original v17.9 archive.
                matrix=cached[tag+'__M'][idx]
                bs={}
                for st,h in heights.items():
                    baseB=cached[tag+'__'+st+'__BB'][idx]
                    recon=[three_frames(history,float(sj[i]),h) for i in idx]
                    cent=np.stack([f['centered'] for f in recon])
                    # Numerical check on the central observer before any claims.
                    basis_diff=np.max(np.linalg.norm(cent-baseB,axis=(1,2)))
                    projector_diff=np.max(np.linalg.norm(
                        cent@cent.transpose(0,2,1)-baseB@baseB.transpose(0,2,1),axis=(1,2)))
                    regression.append(dict(tolerance=tol,eps=er,sign=sign,stencil=st,
                        central_frame_max_abs=basis_diff,central_projector_max_abs=projector_diff))
                    if basis_diff>args.frame_gate:
                        raise RuntimeError(f'central frame reproduction fails {tag}/{st}: {basis_diff:.6g}')
                    bs[st]={
                       'BA':cached[tag+'__'+st+'__BA'][idx],
                       'M':matrix,
                       'BB':{'centered':baseB,
                             'forward':np.stack([f['forward'] for f in recon]),
                             'backward':np.stack([f['backward'] for f in recon])}
                    }
                branches[sign]=bs
            for st,h in heights.items():
                op={}
                for mode in ('centered','forward','backward','non_straddling'):
                    tri=np.full((len(idx),3,11),np.nan)
                    delt=np.full((len(idx),11),np.nan)
                    for j,i in enumerate(idx):
                        if mode=='non_straddling':
                            m='backward' if xi[i]<crossing else 'forward'
                            if abs(float(xi[i])-crossing)<1e-11:
                                # The exact event point has a branch-selection convention.
                                continue
                        else:m=mode
                        neg=branches[-1][st];pos=branches[+1][st]
                        minus=(neg['BA'][j],neg['BB'][m][j],neg['M'][j])
                        plus=(pos['BA'][j],pos['BB'][m][j],pos['M'][j])
                        share,target,error=shapley_three_factors(minus,plus)
                        if error>1e-10:raise RuntimeError(f'Shapley closure error {error}')
                        tri[j]=share;delt[j]=target
                    op[mode]=(tri,delt)
                orig_key=f'{tol}__kick{args.kick:g}__eps{er:g}__{st}'
                # Exact central-frame Shapley regression at all retained probes.
                original=orig_shap[orig_key+'__shapley_ABC'][idx]
                valid=np.isfinite(original).all(axis=(1,2))
                central_err=np.max(np.linalg.norm(op['centered'][0][valid]-original[valid],axis=2))
                regression.append(dict(tolerance=tol,eps=er,sign='pair',stencil=st,
                    centered_shapley_max_abs=central_err))
                if central_err>args.shapley_gate:
                    raise RuntimeError(f'central Shapley mismatch {orig_key}: {central_err:.6g}')
                # Frozen 11D direction is extracted once from historical central A,
                # at the first eligible point strictly after the B crossing.
                ia=int(np.where((xi[idx]>crossing+0.025)&(xi[idx]<crossing+0.035))[0][0])
                frozenA=safe_unit(op['centered'][0][ia,0])
                for mode,(tri,delt) in op.items():
                    npz_outputs[orig_key+'__'+mode+'__shapley_ABC']=tri
                    npz_outputs[orig_key+'__'+mode+'__dlocal_delta']=delt
                    for j,i in enumerate(idx):
                        if not np.isfinite(tri[j]).all():continue
                        A,B,M=tri[j];S=A+B;D=delt[j]
                        An=np.linalg.norm(A);Bn=np.linalg.norm(B)
                        Sdir=safe_unit(A)
                        rows.append(dict(tolerance=tol,eps_rel=er,stencil=st,observer=mode,
                            xi=float(xi[i]),relative_xi=float(xi[i]-crossing),
                            angle_AB_deg=angle_deg(A,B),
                            A_norm=float(An),B_norm=float(Bn),S_norm=float(np.linalg.norm(S)),
                            M_norm=float(np.linalg.norm(M)),delta_norm=float(np.linalg.norm(D)),
                            AB_cancel=float(np.linalg.norm(S)/max(An+Bn,1e-300)),
                            S_along_A=float(S@Sdir),M_along_A=float(M@Sdir),
                            D_along_A=float(D@Sdir),D_along_frozenA=float(D@frozenA),
                            S_transverse_to_A=float(np.linalg.norm(S-(S@Sdir)*Sdir)),
                            D_transverse_to_frozenA=float(np.linalg.norm(D-(D@frozenA)*frozenA)),
                            closure=float(np.linalg.norm(A+B+M-D)),
                            support_h_in_s=float(h),support_halfwidth_xi=float(h/b1)))
                print('completed',tol,er,st, 'n',len(idx), 'regression',central_err,flush=True)
    df=pd.DataFrame(rows).sort_values(['tolerance','eps_rel','stencil','observer','xi'])
    reg=pd.DataFrame(regression)
    if df.empty:raise RuntimeError('no valid Shapley points')
    df.to_csv(args.out/'observer_comparison.csv',index=False,float_format='%.16g')
    reg.to_csv(args.out/'central_regression.csv',index=False,float_format='%.16g')
    npz_outputs['xi']=xi[idx]
    np.savez_compressed(args.out/'observer_shapley_arrays.npz',**npz_outputs)
    # Compare both stencils & all tolerances at the same event offsets.
    primary=df[(df.tolerance=='standard')&(df.eps_rel==1e-8)]
    for st in heights:
        piece=primary[primary.stencil==st]
        fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True)
        for mode in ('centered','forward','backward','non_straddling'):
            sl=piece[piece.observer==mode]
            axes[0].plot(sl.relative_xi,sl.B_norm,label=mode,lw=1.2)
            axes[1].plot(sl.relative_xi,sl.angle_AB_deg,label=mode,lw=1.2)
            axes[2].plot(sl.relative_xi,sl.D_along_frozenA,label=mode,lw=1.2)
        for ax in axes:
            ax.axvline(0,color='gray',linestyle=':',linewidth=0.8)
            ax.axvspan(-heights[st]/b1,heights[st]/b1,color='gray',alpha=.09)
            ax.grid(alpha=.17)
            ax.set_xlim(-args.micro_radius,args.micro_radius)
        axes[0].set_ylabel('||Shapley B||')
        axes[1].set_ylabel('angle(A,B) [deg]')
        axes[2].set_ylabel('Delta along frozen A')
        axes[2].set_xlabel('xi - B crossing')
        axes[0].legend(loc='best',ncol=2,fontsize='small')
        fig.suptitle(f'Q550 kick203 — B-only stencil change, {st} / standard / eps=1e-8')
        fig.tight_layout()
        fig.savefig(args.out/f'B_only_observer_comparison_{st}.png',dpi=180)
        plt.close(fig)
    # A tabulation at symmetric event offsets, not a claim about continuum limits.
    probes=[]
    for st in heights:
      for off in (-.025,-.01,-.005,-.001,-.00001,.00001,.001,.005,.01,.025):
        for mode in ('centered','forward','backward','non_straddling'):
          q=primary[(primary.stencil==st)&(primary.observer==mode)]
          q=q.iloc[np.argmin(abs(q.relative_xi-off))] if len(q) else None
          if q is not None and abs(q.relative_xi-off)<3e-7:
            probes.append(dict(stencil=st,observer=mode,requested_relative_xi=off,
                B_norm=float(q.B_norm),angle_AB_deg=float(q.angle_AB_deg),
                delta_norm=float(q.delta_norm),delta_signed_frozen=float(q.D_along_frozenA)))
    pd.DataFrame(probes).to_csv(args.out/'symmetric_offset_comparison.csv',index=False)
    # Simple cross-stencil descriptive diagnosis; no induced physics assertion.
    gates={st:float(h/b1) for st,h in heights.items()}
    summary=dict(version=VERSION,crossing_xi=float(crossing),selected_points=len(idx),
       full_archive_points=n,heights_support_s=heights,centered_support_halfwidth_xi=gates,
       stencil_halfwidth_ratio=gates['hhalf']/gates['hquarter'],
       max_central_frame_difference=float(reg.central_frame_max_abs.max()),
       max_central_projector_difference=float(reg.central_projector_max_abs.max()),
       max_central_shapley_difference=float(reg.centered_shapley_max_abs.max()),
       max_shapley_closure=float(df.closure.max()),source_zip_sha=original_hash,
       source_provenance=provenance,
       reading=('B-only centered/forward/backward/non-straddling comparison reuses the '
                'same archived full11D M and section-A frames. The non-straddling '
                'observer is backward for xi<xiB, forward for xi>xiB and undefined '
                'at xiB. Shapley is off-shell. No material torsion or physical impulse '
                'can be inferred from an observer-dependent band.'))
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (args.out/'REPORT.md').write_text(f"""# Q550 kick203 centered/one-sided B-observer audit v18.0

B crossing `xi_B = {crossing:.15f}`. Selected `{len(idx)}` of original `{n}`
points. No re-integration of the saved 11D variational cocycle M.

Observer differences (t is the unit 11D state tangent):

- centered: `d=(t(s+h)-t(s-h))/(2h)`;
- forward: `d=(t(s+h)-t(s))/h`;
- backward: `d=(t(s)-t(s-h))/h`;
- non-straddling: backward strictly before B, forward strictly after B.

Only the **B frame** is modified. A and M remain the original v17.9
arrays in all allocations. Same +/- nonlinear histories, same perturbations,
same Q550 parameters, same source provenance.

Expected centered support halfwidth in xi:
- h/2: `{gates['hhalf']:.12g}`
- h/4: `{gates['hquarter']:.12g}`

Strict reconstruction checks:
- max central B frame difference: `{summary['max_central_frame_difference']:.3e}`
- max centered Shapley vector difference: `{summary['max_central_shapley_difference']:.3e}`
- max algebraic closure error: `{summary['max_shapley_closure']:.3e}`

## Boundaries

The archived v17.9 M and A frames are frozen **by design**, so any effect in
this audit is attributable to the B-observer construction, not to a changed
transport. The Shapley hybrid histories remain off-shell. This is not a
certification of physical flip or a force. The exact event node is undefined
for the non-straddling observer; one-sided offsets are compared instead.

See `observer_comparison.csv`, `symmetric_offset_comparison.csv`, figures,
`summary.json` and `central_regression.csv`.
""")
    print(json.dumps({k:v for k,v in summary.items() if k!='source_provenance'},indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base-root',type=Path,default=Path('/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle'))
    ap.add_argument('--source-bundle-zip',type=Path,default=Path('/root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip'))
    ap.add_argument('--v179-sendback',type=Path,default=ROOT/'inputs/q550_kick203_event_micro_v17_9_sendback.zip')
    ap.add_argument('--singular-vectors',type=Path,default=ROOT/'inputs/singular_direction_vectors.csv')
    ap.add_argument('--out',type=Path,default=ROOT/'results/kick203_v18_0')
    ap.add_argument('--kick',type=float,default=203.)
    ap.add_argument('--micro-radius',type=float,default=.04)
    ap.add_argument('--eps',type=float,nargs='+',default=[1e-8,3e-8])
    ap.add_argument('--tolerance-sets',nargs='+',choices=['standard','tight'],default=['standard','tight'])
    ap.add_argument('--frame-gate',type=float,default=1e-4)
    ap.add_argument('--shapley-gate',type=float,default=5e-4)
    args=ap.parse_args()
    run(args)
