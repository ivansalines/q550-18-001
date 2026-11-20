#!/usr/bin/env python3
"""Q550 v18.7: dense source-birth sign crossing and tangent co-moving transport.

The Q550 11D baseline and its analytic-Jacobian cocycle are replayed from
original pinned sources. The source Hessian and source-frame integrand are
REUSED from v18.6, not refit from a sparse sign chart. This is explicitly an
OFF_CONSTRAINT_DIAGNOSTIC, not a proof of material hysteresis or rigid torsion.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,os,sys,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
import q550_state_impulse_response_v13 as v13
import history_bridge_v16 as old
from sign_crossing_math import (cumulative_simpson,continuous_integral,roots_from_grid,
 polar_rotation,angle_deg,sign_and_metric_rates)

VERSION='2026-10-16-q550-v18.7-material-birth-crossing'
TOLS={'standard':(1e-12,1e-14),'tight':(3e-13,3e-15)}
DELAYS=[.25,.5,1.,2.,4.,8.]
HISTORY=ROOT/'inputs/v18_6_original_sendback.zip'


def digest(f):
    h=hashlib.sha256()
    with open(f,'rb') as inp:
        for piece in iter(lambda:inp.read(1<<20),b''):h.update(piece)
    return h.hexdigest()


def verify_history():
    doc=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    expected=doc['v18_6_sha256']
    if digest(HISTORY)!=expected:raise RuntimeError('v18.6 frozen sendback SHA256 mismatch')
    with zipfile.ZipFile(HISTORY) as z:
        bad=z.testzip()
        if bad:raise RuntimeError('v18.6 ZIP member damaged: '+bad)
        man=json.loads(z.read('SHA256SUMS_OUTPUT.json'))['files']
        for name,val in man.items():
            if hashlib.sha256(z.read(name)).hexdigest()!=val:raise RuntimeError('v18.6 manifest mismatch '+name)
        for tol in TOLS:
            with np.load(io.BytesIO(z.read(f'checkpoint/{tol}_source_transport.npz')),allow_pickle=False) as d:
                if d['source_integrand_11D'].shape!=(641,11):raise ValueError('wrong source-integrand dimensions')
    return {'history_sha256':expected,'source_manifest_entries':len(man)}


def read_source(tol):
    with zipfile.ZipFile(HISTORY) as z:
        with np.load(io.BytesIO(z.read(f'checkpoint/{tol}_source_transport.npz')),allow_pickle=False) as a:
            return {k:a[k].copy() for k in a.files}


def atomic_save(path,**arrays):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.stem+'.incomplete.npz')
    np.savez_compressed(temp,**arrays);os.replace(temp,path)


def written_csv(file,rows):
    if not rows:return
    file.parent.mkdir(parents=True,exist_ok=True)
    with open(file,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def row_data(tol,src,xi,P,J,integral,rate_h):
    """Dense invariants and co-moving/polar angles, evaluated on saved source grid."""
    slope=float(src['source_scale'][0]);I1=integral[int(round(1./(xi[1]-xi[0])))]
    old=np.einsum('nij,j->ni',P,I1)
    total=np.einsum('nij,nj->ni',P,integral)
    new=total-old
    closure=np.linalg.norm(old+new-total,axis=1)
    ref_idx=int(round(1./(xi[1]-xi[0])))
    R1=polar_rotation(P[ref_idx]);ref_l=old[ref_idx];ref_r=R1.T@ref_l
    rows=[]; old_corot=[];condition=[];normed=[]
    for i,t in enumerate(xi):
        R=polar_rotation(P[i]);oo=R.T@old[i];old_corot.append(oo)
        condition.append(float(np.linalg.cond(P[i])))
        if i<ref_idx:continue
        tn=float(np.linalg.norm(total[i]));en=float(np.linalg.norm(old[i]))
        n,metric,birth=sign_and_metric_rates(P[i],J[i],slope,I1,integral[i],src['source_integrand_11D'][i])
        nom=float(old[i]@total[i]);
        if abs(n-nom)>1e-8*max(1,abs(nom)):raise RuntimeError('metric decomposition inconsistent with physical dot')
        if min(en,tn)<=1e-25:raise RuntimeError('degenerate source contribution')
        co_angle=angle_deg(I1,integral[i])
        angles=angle_deg(old[i],total[i])
        # Gauge-free under a simultaneous orthogonal rotation of material state.
        rows.append(dict(tolerance=tol,xi=float(t),old_norm=en,new_norm=float(np.linalg.norm(new[i])),
            total_norm=tn,old_dot_total=nom,old_signed_percent=100*nom/tn**2,
            angle_old_total_degrees=angles,angle_early_late_degrees=angle_deg(old[i],new[i]),
            old_vs_at_one_lab_angle_degrees=angle_deg(old[i],ref_l),
            old_vs_at_one_polar_corot_angle_degrees=angle_deg(oo,ref_r),
            original_source_frame_early_total_angle_degrees=co_angle,
            physical_vs_source_frame_angle_delta=angles-co_angle,
            dot_rate_metric=metric,dot_rate_new_source=birth,dot_rate_total=metric+birth,
            cocycle_condition=condition[-1],material_11D_closure=closure[i],
            source_integrand_norm=float(np.linalg.norm(src['source_integrand_11D'][i]))))
        normed.append(nom/tn**2)
    return rows,old,total,new,np.asarray(old_corot)


def bracket_and_rates(tol,src,xis,mat,state,jac,branch,early,g_nodes=None,limit=8):
    slope=float(src['source_scale'][0]);s0=float(src['source_s'][0]);
    g=src['source_integrand_11D'] if g_nodes is None else g_nodes
    interp=continuous_integral(xis,g,branch,slope)
    def P_fn(x):return mat(s0+float(x)*slope)
    def dot_fn(x):
        P=P_fn(x);return float((P@early)@(P@interp(x)))
    roots,values=roots_from_grid(xis,dot_fn,start=1.,end=limit)
    output=[]
    for index,r in enumerate(roots):
        s=s0+r*slope;P=P_fn(r);T=interp(r);gg=interp.derivative()(r)/slope
        n,metric,birth=sign_and_metric_rates(P,jac(s,state(s)),slope,early,T,gg)
        dxi=0.0005
        fd=(dot_fn(r+dxi)-dot_fn(r-dxi))/(2*dxi)
        if abs(fd-(metric+birth))>3e-3*max(1.,abs(fd)):
            raise RuntimeError(f'Crossing derivative rate closure failed: {fd=} {metric+birth=}')
        total=P@T;old=P@early;num_old=float(np.linalg.norm(old));num_total=float(np.linalg.norm(total))
        def pct(x):
            M=P_fn(x);q=M@interp(x);e=M@early;return 100.*float(e@q)/float(q@q)
        output.append(dict(tolerance=tol,crossing_index=index+1,
            xi_root=float(r),absolute_s=float(s),dot_at_root=n,
            old_norm=num_old,total_norm=num_total,
            old_total_angle_deg=angle_deg(old,total),
            old_signed_pct_left_0_05=pct(max(1.,r-.05)),
            old_signed_pct_right_0_05=pct(min(8.,r+.05)),
            dot_rate_metric=float(metric),dot_rate_source=float(birth),
            dot_rate_sum=float(metric+birth),dot_rate_fd=float(fd),
            metric_fraction_of_net_derivative=(float(metric/(metric+birth)) if abs(metric+birth)>1e-12 else None),
            source_fraction_of_net_derivative=(float(birth/(metric+birth)) if abs(metric+birth)>1e-12 else None)))
    return output


def analyze_tol(tol,src,Pmat,state,jac,args):
    s0=float(src['source_s'][0]);slope=float(src['source_scale'][0]);xi=src['xi_nodes'];g=src['source_integrand_11D']
    if not np.allclose(np.diff(xi),.0125,atol=1e-12):raise RuntimeError('unexpected xi source grid')
    mat=lambda s:np.asarray(Pmat.sol(float(s)),float).reshape(11,11)
    P=np.asarray([mat(s0+float(t)*slope) for t in xi]);J=np.asarray([jac(s0+float(t)*slope,state(s0+float(t)*slope)) for t in xi])
    if not np.all(np.isfinite(P)) or not np.all(np.isfinite(J)):raise RuntimeError('nonfinite replay')
    cond=np.asarray([np.linalg.cond(t) for t in P])
    if np.max(cond)>args.condition_gate:raise RuntimeError('cocycle condition gate exceeded')
    if np.max(abs(cond-src['cocycle_condition_nodes']))>1e-4*max(1,float(cond.max())):
        raise RuntimeError('v18.6 cocycle conditioning mismatch')
    def reconstruct(gsub,source_xi=xi):
        return cumulative_simpson(source_xi,gsub,slope)
    I=reconstruct(g);I_half=reconstruct(src['source_integrand_halfstep_11D']);
    I_coarse=reconstruct(g[::2],xi[::2]);xcoarse=xi[::2]
    reference_delays=src['predicted_quadratic_at_delays_11D']
    lag_errs=[]
    for j,d in enumerate(DELAYS):
        i=int(round(d/.0125));estimate=P[i]@I[i]
        rel=float(np.linalg.norm(estimate-reference_delays[j])/max(1e-15,np.linalg.norm(reference_delays[j])))
        lag_errs.append(rel)
        if rel>args.replay_relative_gate:raise RuntimeError(f'{tol} 11D v18.6 anchor regression failed at {d}: {rel}')
    rows,early_v,total_v,late_v,corot_v=row_data(tol,src,xi,P,J,I,None)
    early_idx=int(round(1./.0125));early=I[early_idx]
    standard_root=bracket_and_rates(tol,src,xi,mat,state,jac,I,early)
    root_half=bracket_and_rates(tol,src,xi,mat,state,jac,I_half,I_half[early_idx],g_nodes=src['source_integrand_halfstep_11D'])
    # 0.05-xi Simpson cumulant control: same original 641 source evaluations, nested coarse grid.
    root_coarse=bracket_and_rates(tol,src,xcoarse,mat,state,jac,I_coarse,I_coarse[int(round(1./.025))],g_nodes=g[::2])
    # NOTE root_coarse use x[::2] so spacing .025 and source g is still original 641-grid.
    # Its integral used coarser .05-xi cells.
    compare={
       'half_hessian_step_roots_xi':[r['xi_root'] for r in root_half],
       'coarse_simpson_roots_xi':[r['xi_root'] for r in root_coarse],
       'fine_roots_xi':[r['xi_root'] for r in standard_root]
    }
    if len(standard_root)!=len(root_half) or len(standard_root)!=len(root_coarse):
        raise RuntimeError('root count changes under quadrature or Hessian-step controls')
    diffs=[{'fine_root_xi':r['xi_root'],
             'half_hessian_difference_xi':float(root_half[k]['xi_root']-r['xi_root']),
             'coarse_quadrature_difference_xi':float(root_coarse[k]['xi_root']-r['xi_root'])}
            for k,r in enumerate(standard_root)]
    sample_i=np.asarray([int(round(d/.0125)) for d in DELAYS])
    def signed(olde,tot):return float(100.*(olde@tot)/(tot@tot))
    anchor_sign=[{'xi':float(xi[i]),'signed_percent':signed(early_v[i],total_v[i])} for i in sample_i if xi[i]>=1]
    # The historical 4/8 xi signs must be reproduced by the regenerated material cocycle.
    if not(anchor_sign[-2]['signed_percent']>0 and anchor_sign[-1]['signed_percent']<0):
        raise RuntimeError('v18.6 sign-change anchor did not reproduce')
    # Source-frame old vector is exactly fixed in the initial tangent coordinate.
    drift=max(float(np.linalg.norm(np.linalg.solve(P[i],early_v[i])-early))
              /max(1.,float(np.linalg.norm(early))) for i in range(early_idx,len(xi),10))
    summary={'version':VERSION,'tolerance':tol,'constraint_status':'OFF_CONSTRAINT_DIAGNOSTIC',
             'full_11D_anchor_max_relative_error':max(lag_errs),'max_cocycle_condition':float(cond.max()),
             'old_co_moving_inverse_drift_max':drift,
             'dense_nodes_from_xi1_to_xi8':len(rows),'roots':standard_root,
             'root_controls':diffs,'anchors':anchor_sign,
             'warning':'Crossing defined by an 11D Euclidean inner product; not a rigid-body flip or constitutive hysteresis.'}
    atomic_save(ROOT/'checkpoint'/f'{tol}_sign_crossing.npz',
          xi_nodes=xi,P_full11D=P,J_full11D=J,integral_source_frame_11D=I,
          source_integrand_11D=g,old_source_material_11D=early_v,
          late_source_material_11D=late_v,quadratic_total_material_11D=total_v,
          old_source_polar_corotated_11D=corot_v,xi_roots=np.asarray([r['xi_root'] for r in standard_root]),
          v186_anchor_relative_errors=np.asarray(lag_errs),
          version=np.asarray([VERSION]))
    written_csv(ROOT/'results'/f'{tol}_dense_projection.csv',rows)
    written_csv(ROOT/'results'/f'{tol}_crossing_roots.csv',standard_root)
    (ROOT/'results'/f'{tol}_crossing_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return summary


def run_tol(tol,args):
    checkpoint=ROOT/'checkpoint'/f'{tol}_sign_crossing.npz'
    if checkpoint.exists() and not args.force:
        with np.load(checkpoint) as z:
            if str(z['version'][0])!=VERSION:raise RuntimeError('stale checkpoint version; use --force')
        if not (ROOT/'results'/f'{tol}_crossing_summary.json').is_file():raise RuntimeError('missing result summary; --force')
        print(f'REUSE {tol}: {checkpoint}',flush=True);return
    rt,at=TOLS[tol];src=read_source(tol)
    proven=v13.validate_inputs(args.base_root,args.source_bundle_zip)
    pa=v13.load_module(args.base_root/v13.DEPENDENCIES['audit'][0],f'v187_audit_{tol}')
    sup=v13.load_module(args.base_root/v13.DEPENDENCIES['support'][0],f'v187_support_{tol}')
    tube=v13.load_module(args.base_root/v13.DEPENDENCIES['tube'][0],f'v187_tube_{tol}')
    params=pa.ModelParameters(gamma_rho=.3);kappa0=(55./6.)*(10./550.)
    def rhs(s,y):return np.asarray(sup.support_rhs_float(np.asarray(y,float),kappa0=kappa0,params=params,audit=pa),float)
    def jac(s,y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=kappa0,params=params,audit=pa,base=sup),float)
    orbit=pa.load_orbit(550.,args.base_root/v13.DEPENDENCIES['orbit'][0])
    tail=pa.collect_uniform_tail(orbit,tail_periods=70,steps_per_period=160,
          samples_per_rotation=720,initial_x_perturbation=1e-6,params=params,progress_every=1000000)
    spline=pa.build_state_spline(tail)
    seed=max(float(tail.s[0]),float(tail.s[-1])-.5)
    s0=float(src['source_s'][0]);slope=float(src['source_scale'][0]);finish=s0+8.001*slope
    if seed>=s0:raise RuntimeError('baseline seed occurs after source')
    base=old.dense_ivp(rhs,seed,finish,np.asarray(spline(seed),float),rt,at)
    def state(s):return np.asarray(base.sol(float(s)),float).reshape(11)
    state_err=float(np.linalg.norm(state(s0)-src['source_state_11D'])/max(1.,float(np.linalg.norm(src['source_state_11D']))))
    if state_err>args.state_relative_gate:raise RuntimeError(f'v18.6 initial state mismatch: {state_err}')
    def tangent_rhs(s,p):return (jac(s,state(s))@p.reshape(11,11)).ravel()
    tangent=old.dense_ivp(tangent_rhs,s0,finish,np.eye(11).ravel(),rt,at)
    result=analyze_tol(tol,src,tangent,state,jac,args)
    (ROOT/'checkpoint'/f'{tol}_provenance.json').write_text(json.dumps({
        'version':VERSION,'mode':'OFF_CONSTRAINT_DIAGNOSTIC','original_reference_hashes':proven,
        'source_replay_relative_error':state_err,'input_v18_6_sha256':digest(HISTORY),
        'tolerance':tol,'solver_rtol':rt,'solver_atol':at},indent=2)+'\n')
    print(f'{tol}: found {len(result["roots"])} roots: {[round(r["xi_root"],9) for r in result["roots"]]}',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base-root',type=Path,default=v13.DEFAULT_BASE)
    ap.add_argument('--source-bundle-zip',type=Path,default=v13.DEFAULT_SOURCE)
    ap.add_argument('--tolerances',nargs='+',default=['standard','tight'],choices=list(TOLS))
    ap.add_argument('--replay-relative-gate',type=float,default=5e-8)
    ap.add_argument('--state-relative-gate',type=float,default=5e-9)
    ap.add_argument('--condition-gate',type=float,default=1e8)
    ap.add_argument('--force',action='store_true')
    args=ap.parse_args()
    print('VERIFIED v18.6 history:',verify_history(),flush=True)
    for t in args.tolerances:run_tol(t,args)

if __name__=='__main__':main()
