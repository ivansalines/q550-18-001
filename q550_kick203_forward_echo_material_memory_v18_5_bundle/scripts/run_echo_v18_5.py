#!/usr/bin/env python3
"""Q550 v18.5 forward write/transport-compensated reverse-kick experiment.

A second *kick* is reversed, NOT the differential equation's time direction.
This is a nonlinear history/echo diagnostic, NOT proof of physical hysteresis.
No user-supplied invariant map exists. All kicks are OFF-CONSTRAINT.
"""
from __future__ import annotations
import argparse, json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
import q550_state_impulse_response_v13 as v13
import history_bridge_v16 as old
from echo_math import AMPS,DELAYS_XI,READ_AGES_XI,TOLS,compensate,transport_from_origin

VERSION='2026-10-16-q550-kick203-v18.5-forward-echo'


def verified_inputs():
    import hashlib
    m=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    for filename,key in [('q550_kick203_five_amplitude_scaling_v18_4_sendback.zip','historical_v18_4_sha256'),
                         ('singular_direction_vectors.csv','singular_direction_sha256')]:
        p=ROOT/'inputs'/filename
        if hashlib.sha256(p.read_bytes()).hexdigest()!=m[key]:
            raise RuntimeError(f'Pinned input differs: {filename}')
    return m


def get_layout(base):
    events=pd.read_csv(base/v13.DEPENDENCIES['events'][0]).set_index('event_n')
    a0=float(events.loc[390,'support_center'])
    a1=float(events.loc[433,'support_center'])-a0
    b0=float(events.loc[472,'support_center'])
    b1=float(events.loc[513,'support_center'])-b0
    if min(a1,b1)<=0:raise RuntimeError('non-positive historical section slope')
    # Original v18.4 write occurs where B catches source kick 203.
    xi_cross=210.8611366031777
    s0=a0+a1*203.
    if abs((b0+b1*xi_cross)-s0)>1e-8:
        raise RuntimeError('v18.4 B crossing / A kick time is not reproduced')
    return s0,a1,xi_cross


def _checkpoint(args,tol,eps,lag):
    return args.out/f'{tol}_eps{eps:g}_delay{lag:g}.npz'


def save_atomic(path, data):
    tmp=path.with_suffix('.incomplete.npz')
    np.savez_compressed(tmp,**data)
    os.replace(tmp,path)


def compute(args):
    verified_inputs()
    if any(v<=0 for v in args.amps) or any(v<=0 for v in args.delays) or any(v<0 for v in args.read_ages):
        raise ValueError('amplitudes/delays must be positive, read ages non-negative')
    if 0 not in args.read_ages:raise ValueError('read ages must include 0 (immediate echo)')
    if sorted(set(args.read_ages))!=list(args.read_ages):raise ValueError('read ages must be sorted and unique')
    if sorted(set(args.delays))!=list(args.delays):raise ValueError('delays must be sorted and unique')
    args.out.mkdir(parents=True,exist_ok=True)
    provenance=v13.validate_inputs(args.base_root,args.source_bundle_zip)
    pa=v13.load_module(args.base_root/v13.DEPENDENCIES['audit'][0],'v185_audit')
    sup=v13.load_module(args.base_root/v13.DEPENDENCIES['support'][0],'v185_support')
    tube=v13.load_module(args.base_root/v13.DEPENDENCIES['tube'][0],'v185_tube')
    params=pa.ModelParameters(gamma_rho=.3)
    k0=(55./6.)*(10./550.)
    def rhs(s,y):return np.asarray(sup.support_rhs_float(np.asarray(y,float),kappa0=k0,params=params,audit=pa),float)
    def jac(s,y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=k0,params=params,audit=pa,base=sup),float)
    s0,a1,xcross=get_layout(args.base_root)
    kick_dir=old.columns_for_vectors(args.singular_vectors,203.,[1])[1]
    lag_xi=np.array(args.delays,float)
    age_xi=np.array(args.read_ages,float)
    lag_s=s0+a1*lag_xi
    finish=s0+a1*(max(lag_xi)+max(age_xi))
    configs=[(tol,eps,lag) for tol in args.tolerances for eps in args.amps for lag in args.delays]
    if not args.force and all(_checkpoint(args,*c).exists() for c in configs):
        print('All requested checkpoints already exist: no Q550 replay required.',flush=True)
        return
    # Costly archived orbit load/tail is shared by all experiments.
    orb=pa.load_orbit(550.,args.base_root/v13.DEPENDENCIES['orbit'][0])
    tail=pa.collect_uniform_tail(orb,tail_periods=70,steps_per_period=160,
                 samples_per_rotation=720,initial_x_perturbation=1e-6,
                 params=params,progress_every=1000000)
    spline=pa.build_state_spline(tail)
    seed=max(float(tail.s[0]),float(tail.s[-1])-.5)
    if seed>=s0:raise RuntimeError('orbit seed exceeds kick source')
    for tol in args.tolerances:
        tasks=[(eps,lag) for eps in args.amps for lag in args.delays
               if args.force or not _checkpoint(args,tol,eps,lag).exists()]
        if not tasks:continue
        rt,at=TOLS[tol]
        baseline=old.dense_ivp(rhs,seed,finish,np.asarray(spline(seed),float),rt,at)
        y0=np.asarray(baseline.sol(s0),float).reshape(11)
        scale=max(1.,np.linalg.norm(y0))
        # Control gate: frozen variational seed, 11x11 transport on baseline.
        def variational(s,p):
            y=np.asarray(baseline.sol(s),float).reshape(11)
            return (jac(s,y)@p.reshape(11,11)).ravel()
        Phi=old.dense_ivp(variational,s0,finish,np.eye(11).ravel(),rt,at)
        def getP(s):return np.asarray(Phi.sol(float(s)),float).reshape(11,11)
        if np.linalg.norm(getP(s0)-np.eye(11))>5e-12:raise RuntimeError('identity source cocycle lost')
        # The exact zero-lag write then anti-write is y0 in floating arithmetic.
        for sign in (-1,1):
            trial=y0+sign*(args.amps[0]*scale)*kick_dir
            z=compensate(trial,y0,kick_dir,args.amps[0]*scale,sign)
            if np.linalg.norm(z-y0)>1e-11:raise RuntimeError('zero-lag cancellation failed')
        fd_h=3e-7*scale
        fd=(rhs(s0,y0+fd_h*kick_dir)-rhs(s0,y0-fd_h*kick_dir))/(2*fd_h)
        jac_err=float(np.linalg.norm(fd-jac(s0,y0)@kick_dir)/max(1e-12,np.linalg.norm(fd)))
        if jac_err>2e-3:raise RuntimeError(f'Q550 Jacobian directional check failed: {jac_err}')
        for eps in args.amps:
            active_lags=[lag for e,lag in tasks if e==eps]
            if not active_lags:continue
            amp=eps*scale
            writes={}
            for sign in (-1,1):
                writes[sign]=old.dense_ivp(rhs,s0,finish,y0+sign*amp*kick_dir,rt,at)
            for delay in active_lags:
                sp=s0+delay*a1
                read_times=sp+age_xi*a1
                base_read=np.asarray(baseline.sol(read_times),float).T
                Pc=getP(sp)
                linear_lag=Pc@kick_dir
                cond=float(np.linalg.cond(Pc))
                if not np.isfinite(cond):raise RuntimeError('nonfinite cocycle conditioning')
                if cond>args.condition_gate:print('WARNING: variational transfer ill-conditioned',tol,eps,delay,cond,flush=True)
                result={
                  'version':np.array([VERSION]),'tol':np.array([tol]),
                  'eps':np.array([eps]),'amp_abs':np.array([amp]),
                  'delay_xi':np.array([delay]),'read_age_xi':age_xi,
                  'source_s':np.array([s0]),'pulse_s':np.array([sp]),'read_s':read_times,
                  'a_section_slope':np.array([a1]),'source_state_11D':y0,
                  'kick_direction_11D':kick_dir,
                  'linear_state_at_pulse_11D':amp*linear_lag,
                  'baseline_pulse_11D':np.asarray(baseline.sol(sp),float).reshape(11),
                  'baseline_read_11D':base_read,
                  'pulse_variational_condition':np.array([cond]),
                  'jacobian_fd_relative_error':np.array([jac_err]),
                }
                for sign in (-1,1):
                    branch=writes[sign]
                    yp=np.asarray(branch.sol(sp),float).reshape(11)
                    echo0=compensate(yp,baseline.sol(sp),linear_lag,amp,sign)
                    echoed=old.dense_ivp(rhs,sp,float(read_times[-1]),echo0,rt,at)
                    yread=np.asarray(echoed.sol(read_times),float).T
                    write_read=np.asarray(branch.sol(read_times),float).T
                    resid0=echo0-result['baseline_pulse_11D']
                    pred=[];transcond=[]
                    for sr in read_times:
                        pr,c=transport_from_origin(getP(sr),Pc,resid0,args.condition_gate)
                        pred.append(pr);transcond.append(c)
                    tag='plus' if sign>0 else 'minus'
                    result[f'{tag}_write_read_11D']=write_read
                    result[f'{tag}_echo_read_11D']=yread
                    result[f'{tag}_linearly_transported_echo_residual_11D']=np.asarray(pred)
                    result[f'{tag}_transport_conditions']=np.array(transcond)
                    # A naive anti-kick at the pulse would subtract the original
                    # vector, NOT its material tangent transport.  Immediate
                    # controls require no additional nonlinear solve.
                    result[f'{tag}_naive_antikick_immediate_11D']=yp-sign*amp*kick_dir
                    result[f'{tag}_matched_antikick_immediate_11D']=echo0
                    closure=np.linalg.norm(yread[0]-echo0)
                    if closure>2e-10:raise RuntimeError(f'dense echo initial-state closure {closure}')
                # All pairs generated, only now write checkpoint atomically.
                path=_checkpoint(args,tol,eps,delay)
                save_atomic(path,result)
                print(f'SAVED {path.name}: delay={delay:g} xi eps={eps:g} '
                      f'||r+ pulse||={np.linalg.norm(result["plus_echo_read_11D"][0]-base_read[0]):.6e} '
                      f'cond(P)={cond:.3g}',flush=True)
    (args.out/'run_provenance.json').write_text(json.dumps({
      'version':VERSION,'configuration':dict(eps=args.amps,delay_xi=args.delays,
         read_ages_xi=args.read_ages,tolerances=args.tolerances),
      'xi_to_support': 'delta_support=(s433-s390)*delta_xi',
      'source_s':s0,'A_slope':a1,'B_crossing_xi':xcross,
      'original_Q550_source_provenance':provenance,
      'constraint_status':'OFF_CONSTRAINT_DIAGNOSTIC',
      'reverse_definition':'TRANSPORT_MATCHED_ANTI_KICK_FORWARD_TIME_NOT_TIME_REVERSAL',
      'material_observable':'Q550 real 11-component ODE state; not dlocal descriptor'
    },indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-root',type=Path,default=v13.DEFAULT_BASE)
    p.add_argument('--source-bundle-zip',type=Path,default=v13.DEFAULT_SOURCE)
    p.add_argument('--singular-vectors',type=Path,default=ROOT/'inputs/singular_direction_vectors.csv')
    p.add_argument('--out',type=Path,default=ROOT/'checkpoint')
    p.add_argument('--amps',nargs='+',type=float,default=list(AMPS))
    p.add_argument('--delays',nargs='+',type=float,default=list(DELAYS_XI))
    p.add_argument('--read-ages',nargs='+',type=float,default=list(READ_AGES_XI))
    p.add_argument('--tolerances',nargs='+',choices=list(TOLS),default=list(TOLS))
    p.add_argument('--condition-gate',type=float,default=1e10)
    p.add_argument('--force',action='store_true')
    compute(p.parse_args())
