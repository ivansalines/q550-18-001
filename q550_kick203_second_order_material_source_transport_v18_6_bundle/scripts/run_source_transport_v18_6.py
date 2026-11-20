#!/usr/bin/env python3
"""Compute the material quadratic source, Duhamel transport, and v18.5 echo closure.

The Q550 model is verified against the original pinned SHA256 history.
This is an OFF-CONSTRAINT material-diagnostic, not physical hysteresis proof.
"""
from __future__ import annotations
import argparse,hashlib,io,json,sys,zipfile,os
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
import q550_state_impulse_response_v13 as v13
import history_bridge_v16 as old
from source_transport_math import (source_hessian_via_jac,transported_integrand,
    simpson_cell,signed_even,relative_error,angle_deg,slope_finite_amplitude)
VERSION='2026-10-16-q550-second-order-material-source-v18.6'
TOLS={'standard':(1e-12,1e-14),'tight':(3e-13,3e-15)}
DELAYS=(0.25,0.5,1.,2.,4.,8.)
AGES=(0.,0.25,0.5,1.,2.)
AMPS=(1e-8,3e-8,1e-7)


def read_v185(zip_path,tol,eps,delay):
    name=f'checkpoint/{tol}_eps{eps:g}_delay{delay:g}.npz'
    with zipfile.ZipFile(zip_path) as z:
        with np.load(io.BytesIO(z.read(name)),allow_pickle=False) as q:
            return {key:q[key].copy() for key in q.files}

def atomic_save(path,**kwargs):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.stem+'.incomplete.npz')
    np.savez_compressed(temp,**kwargs)
    os.replace(temp,path)

def verify_history():
    m=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    f=ROOT/'inputs'/m['input_file']
    h=hashlib.sha256(f.read_bytes()).hexdigest()
    if h!=m['input_sha256']:raise RuntimeError('v18.5 pinned sendback SHA256 mismatch')
    with zipfile.ZipFile(f) as z:
        bad=z.testzip()
        if bad:raise RuntimeError(f'v18.5 damaged member: {bad}')
    return f

def layout(base):
    events=pd.read_csv(base/v13.DEPENDENCIES['events'][0]).set_index('event_n')
    s0=float(events.loc[390,'support_center'])+(float(events.loc[433,'support_center'])-float(events.loc[390,'support_center']))*203.
    slope=float(events.loc[433,'support_center'])-float(events.loc[390,'support_center'])
    if slope<=0:raise RuntimeError('nonpositive historical xi scale')
    return s0,slope

def process_tol(tol,args,history):
    path=ROOT/'checkpoint'/f'{tol}_source_transport.npz'
    if path.exists() and not args.force:
        print(f'REUSE {path}',flush=True);return
    rt,at=TOLS[tol]
    # Every dependency hash is checked at runtime before importing original code.
    prov=v13.validate_inputs(args.base_root,args.source_bundle_zip)
    pa=v13.load_module(args.base_root/v13.DEPENDENCIES['audit'][0],f'v186_audit_{tol}')
    sup=v13.load_module(args.base_root/v13.DEPENDENCIES['support'][0],f'v186_support_{tol}')
    tube=v13.load_module(args.base_root/v13.DEPENDENCIES['tube'][0],f'v186_tube_{tol}')
    params=pa.ModelParameters(gamma_rho=.3)
    kappa0=(55./6.)*(10./550.)
    def rhs(s,y):return np.asarray(sup.support_rhs_float(np.asarray(y,float),kappa0=kappa0,params=params,audit=pa),float)
    def jac(s,y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=kappa0,params=params,audit=pa,base=sup),float)
    s0,slope=layout(args.base_root)
    reference=read_v185(history,tol,AMPS[1],DELAYS[0]); y_ref=reference['source_state_11D']; v=reference['kick_direction_11D']
    if abs(float(reference['source_s'][0])-s0)>1e-8:raise RuntimeError('source support mismatch')
    if abs(float(reference['a_section_slope'][0])-slope)>1e-10:raise RuntimeError('xi support slope mismatch')
    # Original v18.5 Q550 reconstruction exactly, except only needed up to xi=8.
    orbit=pa.load_orbit(550.,args.base_root/v13.DEPENDENCIES['orbit'][0])
    tail=pa.collect_uniform_tail(orbit,tail_periods=70,steps_per_period=160,
          samples_per_rotation=720,initial_x_perturbation=1e-6,params=params,progress_every=1000000)
    spline=pa.build_state_spline(tail)
    seed=max(float(tail.s[0]),float(tail.s[-1])-.5)
    if seed>=s0:raise RuntimeError('baseline seed after kick')
    finish=s0+10.*slope
    base=old.dense_ivp(rhs,seed,finish,np.asarray(spline(seed),float),rt,at)
    y0=np.asarray(base.sol(s0),float).reshape(11)
    state_err=relative_error(y0,y_ref)
    print(f'{tol}: source state relative replay difference={state_err:.4e}',flush=True)
    if state_err>args.source_state_relative_gate:raise RuntimeError('v18.5 source state replay mismatch')
    # Variational fundamental matrix P(s) anchored to identical material source.
    def tangent_rhs(s,p):
        y=np.asarray(base.sol(s),float).reshape(11)
        return (jac(s,y)@p.reshape(11,11)).ravel()
    tangent=old.dense_ivp(tangent_rhs,s0,finish,np.eye(11).ravel(),rt,at)
    def mat(s):return np.asarray(tangent.sol(float(s)),float).reshape(11,11)
    def state(s):return np.asarray(base.sol(float(s)),float).reshape(11)
    # h is an absolute directional step in material state units; h/2 is independent check.
    h=args.hessian_step
    halfcount=int(round(8./(args.cell_xi/2)))
    if abs(halfcount*args.cell_xi/2-8)>1e-10 or halfcount%2:raise RuntimeError('cell-xi must divide 8')
    dx=args.cell_xi/2
    xis=np.arange(halfcount+1,dtype=float)*dx
    g_main=[];g_half=[];source_main=[];source_half=[];conds=[]
    for i,xi in enumerate(xis):
        s=s0+xi*slope;y=state(s);P=mat(s);u=P@v
        S=source_hessian_via_jac(jac,s,y,u,h)
        S2=source_hessian_via_jac(jac,s,y,u,h/2)
        g,c=transported_integrand(P,S,args.condition_gate)
        g2,_=transported_integrand(P,S2,args.condition_gate)
        g_main.append(g);g_half.append(g2)
        source_main.append(S);source_half.append(S2);conds.append(c)
        if i%80==0:print(f'  {tol}: source xi={xi:g} cond(P)={c:.3g}',flush=True)
    g=np.asarray(g_main); gh=np.asarray(g_half)
    # Simpson fine cells have endpoints xi=cell*n, midpoints=cell*(n+.5).
    nf=halfcount//2; fine=[];hhalf=[]
    for c in range(nf):
        k=2*c
        fine.append(simpson_cell(g[k],g[k+1],g[k+2],slope*args.cell_xi))
        hhalf.append(simpson_cell(gh[k],gh[k+1],gh[k+2],slope*args.cell_xi))
    fine=np.asarray(fine);hhalf=np.asarray(hhalf)
    # Simpson coarse cells double the interval width and share source nodes.
    if nf%2:raise RuntimeError('fine cell count must be even for nested coarse Simpson')
    coarse=[]
    for c in range(nf//2):
        k=4*c
        coarse.append(simpson_cell(g[k],g[k+2],g[k+4],slope*2*args.cell_xi))
    coarse=np.asarray(coarse)
    targets=[];binned=[];norms=[];cosines=[]
    band_edges=np.asarray((0.,0.25,0.5,1.,2.,4.,8.))
    for lag in DELAYS:
        j=int(round(lag/args.cell_xi));
        if abs(j*args.cell_xi-lag)>1e-10:raise RuntimeError('delay must be integer multiple of cell spacing')
        Pend=mat(s0+lag*slope)
        k2=Pend@np.sum(fine[:j],axis=0)
        k2h=Pend@np.sum(hhalf[:j],axis=0)
        k2c=Pend@np.sum(coarse[:j//2],axis=0)
        piece=[]
        for l,r in zip(band_edges[:-1],band_edges[1:]):
            start=int(round(l/args.cell_xi));stop=min(j,int(round(r/args.cell_xi)))
            piece.append(Pend@fine[start:stop].sum(axis=0) if stop>start else np.zeros(11))
        pieces=np.asarray(piece)
        if relative_error(np.sum(pieces,axis=0),k2)>1e-10:raise RuntimeError('source bin reconstruction failed')
        for eps in AMPS:
            q=read_v185(history,tol,eps,lag)
            amp=float(q['amp_abs'][0])
            if relative_error(q['baseline_pulse_11D'],state(s0+lag*slope))>args.source_state_relative_gate:
                raise RuntimeError('v18.5 baseline pulse does not match replay')
            for ai,age in enumerate(AGES):
                s_read=s0+(lag+age)*slope
                propagated=mat(s_read)@np.linalg.solve(Pend,k2)
                baseline=q['baseline_read_11D'][ai]
                observed=signed_even(q['plus_echo_read_11D'][ai],q['minus_echo_read_11D'][ai],baseline)/amp**2
                pred_state=amp**2*propagated
                obs_raw=amp**2*observed
                targets.append(dict(tolerance=tol,epsilon_relative=eps,amplitude_absolute=amp,
                      delay_xi=lag,read_age_xi=age,
                      observed_quadratic_norm=float(np.linalg.norm(observed)),
                      predicted_quadratic_norm=float(np.linalg.norm(propagated)),
                      relative_vector_error=relative_error(pred_state,obs_raw),
                      absolute_vector_error=float(np.linalg.norm(pred_state-obs_raw)),
                      angle_degrees=angle_deg(propagated,observed),
                      simpson_grid_relative_error=relative_error(k2c,k2),
                      hessian_step_relative_error=relative_error(k2h,k2),
                      source_state_replay_relative_error=state_err,
                      pulse_linear_condition=float(q['pulse_variational_condition'][0])))
        binned.append(pieces);norms.append(k2)
    atomic_save(path,version=np.asarray([VERSION]),tol=np.asarray([tol]),
                source_s=np.array([s0]),source_scale=np.array([slope]),
                source_state_11D=y0,kick_direction_11D=v,
                xi_nodes=xis,source_11D=np.asarray(source_main),
                source_halfstep_11D=np.asarray(source_half),
                source_integrand_11D=g,source_integrand_halfstep_11D=gh,
                source_bin_edges=band_edges,source_bin_contributions_11D=np.asarray(binned),
                predicted_quadratic_at_delays_11D=np.asarray(norms),
                cocycle_condition_nodes=np.asarray(conds),
                source_state_replay_relative_error=np.array([state_err]),
                hessian_step=np.array([h]),cell_xi=np.array([args.cell_xi]))
    # This supplemental diagnostic may be rerun from stored arrays without any new orbit.
    (ROOT/'checkpoint'/f'{tol}_observed_comparison.json').write_text(json.dumps(targets,indent=2))
    (ROOT/'checkpoint'/f'{tol}_provenance.json').write_text(json.dumps(dict(version=VERSION,tolerance=tol,
         original_Q550_source_provenance=prov,mode='OFF_CONSTRAINT_DIAGNOSTIC',
         jacobian='original 11D analytic Jacobian, numerical Hessian of Jacobian',
         source_s=s0,xi_to_s=slope,hessian_step=h,cell_xi=args.cell_xi),indent=2))
    print(f'SAVED {path}',flush=True)

def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-root',type=Path,default=v13.DEFAULT_BASE)
    p.add_argument('--source-bundle-zip',type=Path,default=v13.DEFAULT_SOURCE)
    p.add_argument('--tolerances',nargs='+',choices=list(TOLS),default=list(TOLS))
    p.add_argument('--hessian-step',type=float,default=2e-4)
    p.add_argument('--cell-xi',type=float,default=0.025)
    p.add_argument('--condition-gate',type=float,default=1e8)
    p.add_argument('--source-state-relative-gate',type=float,default=5e-9)
    p.add_argument('--force',action='store_true')
    return p.parse_args()

if __name__=='__main__':
    a=parse_args();history=verify_history()
    if a.hessian_step<=0 or a.cell_xi<=0:raise SystemExit('steps must be positive')
    for tol in a.tolerances:process_tol(tol,a,history)
