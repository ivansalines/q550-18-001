#!/usr/bin/env python3
"""v18.4: nonlinear on-flow Q550 history pairs at 3 new amplitude levels.

The two existing amplitude levels are consumed unchanged from v18.3. Recompute
ONLY +/- histories at three new amplitudes, keeping the original xi, A/B
support maps, kick direction, tolerance, A=hquarter and M construction frozen.
The 11D descriptor response is not a physical material-state coordinate.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import os
import sys
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'), str(ROOT/'scripts')]
import history_bridge_v16 as old
import q550_state_impulse_response_v13 as v13
from run_kick203_event_micro_v17_9 import reconstruct
from run_kick203_second_order_observer_v18_1 import tangent_differences
from bridge_math import factor_secants
from scaling_math import AMPS

VERSION='2026-10-16-Q550-kick203-five-amplitudes-v18.4'
NEW_AMPS=(5e-9,1.5e-8,2e-8)


def read_bridge(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        candidates=[n for n in z.namelist() if n.endswith('/bridge_full11D_arrays.npz')]
        if len(candidates)!=1: raise RuntimeError(f'exactly one v18.3 bridge expected, got {candidates}')
        with np.load(io.BytesIO(z.read(candidates[0])),allow_pickle=False) as data:
            return {k:np.asarray(data[k]) for k in data.files}


def verify_source():
    manifest=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    p=ROOT/'inputs/q550_kick203_material_to_observer_bridge_v18_3_bundle.zip'
    digest=v13.hash_file(p)
    if digest != manifest['pinned_input_sha256']:
        raise RuntimeError(f'v18.3 archive SHA256 differs: {digest}')
    return p


def compute(args):
    verify_source()
    reference=read_bridge(ROOT/'inputs/q550_kick203_material_to_observer_bridge_v18_3_bundle.zip')
    xi=np.asarray(reference['xi'],float)
    sA=np.asarray(reference['sA'],float)
    sB=np.asarray(reference['sB'],float)
    s0=float(reference['s0'][0])
    crossing=float(reference['crossing_xi'][0])
    axis=np.asarray(reference['standard__eps1e-08__frozen_axis_11D'],float)
    axis=axis/np.linalg.norm(axis)
    n=len(xi)
    if n!=47 or not np.all(np.diff(xi)>0):raise RuntimeError('frozen 47-sample grid unexpected')
    if abs(crossing-210.8611366031777)>1e-10:raise RuntimeError('frozen crossing has drifted')
    if abs((s0-sB[np.argmin(abs(xi-crossing))]))>5e-10:raise RuntimeError('support crossing mismatch')
    if abs(np.linalg.norm(axis)-1)>1e-12:raise RuntimeError('frozen axis not unit norm')

    # Check that the ODE code & complete numerical history are the original Q550 sources.
    provenance=v13.validate_inputs(args.base_root,args.source_bundle_zip)
    pa=v13.load_module(args.base_root / v13.DEPENDENCIES['audit'][0], 'v184_audit')
    support=v13.load_module(args.base_root / v13.DEPENDENCIES['support'][0], 'v184_support')
    tube=v13.load_module(args.base_root / v13.DEPENDENCIES['tube'][0], 'v184_tube')
    params=pa.ModelParameters(gamma_rho=.3);k0=(55./6.)*(10./550.)
    def rhs(s,y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=k0,params=params,audit=pa),float)
    def jac(s,y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=k0,params=params,audit=pa,base=support),float)
    events=pd.read_csv(args.base_root/v13.DEPENDENCIES['events'][0]).set_index('event_n')
    a0=float(events.loc[390,'support_center']); a1=float(events.loc[433,'support_center'])-a0
    b0=float(events.loc[472,'support_center']); b1=float(events.loc[513,'support_center'])-b0
    if max(abs(sA-(a0+a1*xi)))>5e-10 or max(abs(sB-(b0+b1*xi)))>5e-10:
        raise RuntimeError('A/B map does not reproduce archived support samples')
    shared=json.loads((args.base_root/v13.DEPENDENCIES['shared'][0]).read_text())
    hquarter=.0625*float(shared['p_line'])/(int(shared['arc_samples'])-1)
    h16=hquarter/4
    stop=float(max(sA.max(),sB.max())+max(.1,4*hquarter))
    direction=old.columns_for_vectors(args.singular_vectors,203.,[1])[1]
    args.out.mkdir(parents=True,exist_ok=True)

    # Resume already complete amplitude checkpoints without repeating expensive ODE integration.
    missing=[(tol,eps) for tol in ('standard','tight') for eps in NEW_AMPS
             if args.force or not (args.out/f'{tol}_eps{eps:g}.npz').exists()]
    if not missing:
        print('All six new (tolerance,amplitude) checkpoints already exist; reusing.',flush=True)
        return
    orb=pa.load_orbit(550.,args.base_root/v13.DEPENDENCIES['orbit'][0])
    tail=pa.collect_uniform_tail(orb,tail_periods=70,steps_per_period=160,
                  samples_per_rotation=720,initial_x_perturbation=1e-6,
                  params=params,progress_every=1000000)
    spline=pa.build_state_spline(tail)
    seed=max(float(tail.s[0]),float(tail.s[-1])-.5)
    if seed>=sB.min():raise RuntimeError('baseline seed beyond original first support section')
    for tol in ('standard','tight'):
        epss=[eps for t,eps in missing if t==tol]
        if not epss:continue
        rt,at={'standard':(1e-12,1e-14),'tight':(3e-13,3e-15)}[tol]
        baseline=old.dense_ivp(rhs,seed,stop,spline(seed),rt,at)
        y0=np.asarray(baseline.sol(s0),float).reshape(11)
        yscale=max(1.,float(np.linalg.norm(y0)))
        print('baseline ready',tol,'norm(y0)=',float(np.linalg.norm(y0)),flush=True)
        for eps in epss:
            pair={}
            for sign in (-1,+1):
                branch=old.dense_ivp(rhs,s0,stop,y0+sign*eps*yscale*direction,rt,at)
                history=old.History(baseline,branch,s0,rhs,jac)
                M,conds,sizes,frames=reconstruct(history,xi,sA,sB,{'hquarter':hquarter},rt,at,args.block_xi)
                BA,BB=frames['hquarter']
                # B non-straddling second-order observer is an independent sensitivity check.
                BB2=np.full((n,11,2),np.nan)
                for j,s in enumerate(sB):
                    if abs(s-s0)<5e-12:continue
                    side='backward2' if s<s0 else 'forward2'
                    BB2[j]=tangent_differences(history,float(s),h16)[side]
                # State differences at the SAME physical support s in each history.
                pair[sign]={'BA':BA,'BB':BB,'M':M,'BB2':BB2,
                            'state_A':np.stack([history.state(float(s)) for s in sA]),
                            'state_B':np.stack([history.state(float(s)) for s in sB]),
                            'max_cond':float(np.max(conds))}
                print('  completed',tol,eps,'history',sign,flush=True)
            scenter=[]; interactions=[]; deltas=[]; onedelta=[]
            for j in range(n):
                Fm=tuple(pair[-1][key][j] for key in ('BA','BB','M'))
                Fp=tuple(pair[+1][key][j] for key in ('BA','BB','M'))
                sec,inter,delta=factor_secants(Fm,Fp,old.descriptor)
                scenter.append(sec);interactions.append(inter);deltas.append(delta)
                if np.isfinite(pair[-1]['BB2'][j]).all() and np.isfinite(pair[+1]['BB2'][j]).all():
                    onedelta.append(old.descriptor(pair[+1]['BA'][j],pair[+1]['BB2'][j],pair[+1]['M'][j])-
                                   old.descriptor(pair[-1]['BA'][j],pair[-1]['BB2'][j],pair[-1]['M'][j]))
                else:onedelta.append(np.full(11,np.nan))
            arr={
                'xi':xi,'eps':np.array([eps]),'axis':axis,
                'centered_delta_11D':np.asarray(deltas),
                'centered_interaction_11D':np.asarray(interactions),
                'centered_factor_secants_3x11':np.asarray(scenter),
                'non_straddling_order2_delta_11D':np.asarray(onedelta),
                'material_difference_A_11D':pair[+1]['state_A']-pair[-1]['state_A'],
                'material_difference_B_11D':pair[+1]['state_B']-pair[-1]['state_B'],
                'source_amp_absolute':np.array([eps*yscale]),
                'max_local_cocycle_condition':np.array([max(pair[s]['max_cond'] for s in (-1,1))]),
            }
            if np.max(np.linalg.norm(arr['centered_delta_11D']-arr['centered_factor_secants_3x11'].sum(axis=1)-arr['centered_interaction_11D'],axis=1))>1e-11:
                raise RuntimeError('exact 11D algebraic reconstruction mismatch')
            if np.max(np.linalg.norm(arr['material_difference_B_11D'][sB<s0],axis=1))>1e-14:
                raise RuntimeError('violated B-before-kick causal zero')
            cp=args.out/f'{tol}_eps{eps:g}.npz'
            tmp=cp.with_name(cp.stem+'.incomplete.npz')
            np.savez_compressed(tmp,**arr)
            os.replace(tmp,cp)
            print('SAVED',cp,flush=True)
    (args.out/'run_provenance.json').write_text(json.dumps({
        'version':VERSION,'new_amplitudes':list(NEW_AMPS),
        'frozen_grid_n':n,'crossing_xi':crossing,'epsilon_unit':1e-8,
        'source_provenance':provenance,'constraint_status':'ODE_FLOW_ONLY_NOT_CONSTRAINT_CERTIFIED',
        'frozen_axis_origin':'v18.3 standard eps=1e-8 frozen A',
        'A_frame':'hquarter centered','M':'historical full11D rebased cocycle',
        'B_frame_centered':'hquarter','B_frame_one_sided':'hsixteenth order2',
    },indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-root',type=Path,default=Path('/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle'))
    p.add_argument('--source-bundle-zip',type=Path,default=Path('/root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip'))
    p.add_argument('--singular-vectors',type=Path,default=ROOT/'inputs/singular_direction_vectors.csv')
    p.add_argument('--out',type=Path,default=ROOT/'checkpoint')
    p.add_argument('--block-xi',type=float,default=.75)
    p.add_argument('--force',action='store_true')
    compute(p.parse_args())
