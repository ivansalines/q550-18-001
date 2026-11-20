#!/usr/bin/env python3
"""Independent nonlinear ODE perturbations: centered Jacobian-column validation."""
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'reference'))
from read_original import load
import q550_state_impulse_response_v13 as source

def main():
    p=argparse.ArgumentParser();p.add_argument('--base-root',type=Path,default=source.DEFAULT_BASE);p.add_argument('--epsilon',type=float,default=2e-4);args=p.parse_args()
    data,history=load();d=data['tight'];s0=d['s0'];slope=d['slope'];y0=d['source_state_11D'];omega=2*np.pi
    for role in ('audit','support','tube'):
        n,correct=source.DEPENDENCIES[role];file=args.base_root/n
        if not file.exists():raise FileNotFoundError(file)
        if source.hash_file(file)!=correct:raise RuntimeError('PINNED ODE SHA mismatch: '+n)
    audit=source.load_module(args.base_root/source.DEPENDENCIES['audit'][0],'v18111_original_audit')
    support=source.load_module(args.base_root/source.DEPENDENCIES['support'][0],'v18111_original_rhs')
    model_params=audit.ModelParameters(gamma_rho=.3);k0=(55./6.)*(10./550.)
    def rhs(y):
        return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=k0,params=model_params,audit=audit),dtype=float)
    # Note: only original RHS is used here, not its Jacobian or tangent propagator.
    # One baseline plus +/- four target columns: nine nonlinear trajectories.
    cols=[10,0,5,9];Y0=[y0.copy()]
    for col in cols:
        for sign in (+1,-1):
            xx=y0.copy();xx[col]+=sign*args.epsilon;Y0.append(xx)
    w0=np.concatenate(Y0);m=len(Y0)
    def f(t,w):
        ww=w.reshape(m,11)
        return np.concatenate([rhs(ww[i]) for i in range(m)])
    check_xi=np.array([0.,1.,2.,4.,6.,8.]);times=s0+slope*check_xi
    sol=solve_ivp(f,(times[0],times[-1]),w0,method='DOP853',rtol=4e-12,atol=5e-13,max_step=slope*.025,t_eval=times)
    if not sol.success:raise RuntimeError(sol.message)
    states=sol.y.T.reshape(len(times),m,11)
    rows=[];max_shape_err=0.;max_clock_err=0.;max_rate_err=0.
    for ii,xi in enumerate(check_xi):
        ind=int(round(xi/.0125));P=d['P_full11D'][ind]
        base=states[ii,0]
        max_rate_err=max(max_rate_err,abs(rhs(base)[10]-omega))
        for k,col in enumerate(cols):
            pos=states[ii,1+2*k];neg=states[ii,2+2*k]
            observed=(pos-neg)/(2*args.epsilon)
            expected=P[:,col]
            shape_err=float(np.linalg.norm(observed[:10]-expected[:10])/max(1.,np.linalg.norm(expected[:10])))
            clock_err=float(abs(observed[10]-(1. if col==10 else 0.)))
            max_shape_err=max(max_shape_err,shape_err);max_clock_err=max(max_clock_err,clock_err)
            rows.append({'xi':float(xi),'input_component_zero_based':col,'relative_material_response_error':shape_err,
                         'absolute_phase_response_error':clock_err,'phase_delta_observed':float(observed[10]),
                         'phase_delta_expected':float(expected[10])})
    # Additional independent *off-trajectory* RHS checks: phase derivative always omega,
    # and changing z,theta at initial probe does not change clock rate.
    rng=np.random.default_rng(18111)
    for _ in range(12):
        pert=y0.copy();pert[:10]+=rng.normal(0,1e-3,10);pert[10]+=rng.normal(0,1e-3)
        max_rate_err=max(max_rate_err,abs(rhs(pert)[10]-omega))
    result={'version':'v18.11.1','kind':'CENTERED_NONLINEAR_RHS_ONLY_SENSITIVITY','original_model_SHA256_verified':True,
         'no_baseline_orbit_regeneration':True,'epsilon':args.epsilon,'xi_probes':check_xi.tolist(),
         'perturbed_initial_coordinates_zero_based':cols,'trajectories':m,
         'max_material_response_relative_error':max_shape_err,'max_phase_response_absolute_error':max_clock_err,
         'max_original_rhs_clock_minus_2pi':max_rate_err,'status':'PASSED' if max_shape_err<2e-3 and max_clock_err<3e-7 and max_rate_err<1e-10 else 'FAILED',
         'not_a_nonlinear_material_constraint_certificate':True}
    out=ROOT/'results';out.mkdir(exist_ok=True)
    with (out/'nonlinear_perturbation_profile.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'nonlinear_perturbation_summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))
    if result['status']!='PASSED':raise RuntimeError('NONLINEAR PHASE PROBE FAILED; inspect results')
if __name__=='__main__':main()
