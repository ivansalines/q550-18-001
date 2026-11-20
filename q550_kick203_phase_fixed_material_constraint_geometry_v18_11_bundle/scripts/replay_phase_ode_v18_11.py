#!/usr/bin/env python3
"""Original Q550 ODE replay: independently integrate exact phase-driven 10D flow.

No baseline 400-orbit solve; uses the archived kick203 state and pinned files.
Jacobian blocks: A=Jzz, b=Jztheta, Phi'=A Phi, q'=Aq+b.
"""
from __future__ import annotations
import argparse,ast,csv,json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
import q550_state_impulse_response_v13 as source
from checkpoint_io_v18_11 import load

def save_csv(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base-root',type=Path,default=source.DEFAULT_BASE)
    ap.add_argument('--max-step-xi',type=float,default=.025);arg=ap.parse_args()
    (ROOT/'results').mkdir(exist_ok=True)
    files={}
    for key in ('audit','support','tube'):
        n,sha=source.DEPENDENCIES[key];p=arg.base_root/n
        if not p.is_file():raise FileNotFoundError('Missing original pinned ODE: '+str(p))
        got=source.hash_file(p)
        if got!=sha:raise RuntimeError(f'Model SHA mismatch for {key}: {got} != {sha}')
        files[key]={'relative_name':n,'sha256':got,'server_path':str(p)}
    audit=source.load_module(arg.base_root/source.DEPENDENCIES['audit'][0],'v1811_pinned_audit')
    support=source.load_module(arg.base_root/source.DEPENDENCIES['support'][0],'v1811_pinned_support')
    tube=source.load_module(arg.base_root/source.DEPENDENCIES['tube'][0],'v1811_pinned_tube')
    params=audit.ModelParameters(gamma_rho=.3);k0=(55./6.)*(10./550.)
    def rhs(y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=k0,params=params,audit=audit),float)
    def jac(y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=k0,params=params,audit=audit,base=support),float)
    data,provenance=load();ref=data['tight'];x=ref['xi_nodes'];s0=float(ref['s0']);slope=float(ref['slope']);y0=ref['source_state_11D'];phase0=float(y0[10]);omega=2*np.pi
    # Full nonlinear 11D state and separate, independently integrated 10D state + Phi + phase impulse q.
    start=np.concatenate([y0,y0[:10],np.eye(10).ravel(),np.zeros(10)])
    def coupled(s,w):
        full=w[:11]; z=w[11:21];phi=w[21:121].reshape(10,10);q=w[121:131]
        theta=phase0+omega*(s-s0)
        red=np.concatenate((z,[theta]))
        f=rhs(red); J=jac(red); A=J[:10,:10]; b=J[:10,10]
        return np.concatenate([rhs(full),f[:10],(A@phi).ravel(),A@q+b])
    end=s0+slope*8
    sol=solve_ivp(coupled,(s0,end),start,method='DOP853',rtol=8e-12,atol=1e-13,
                  max_step=slope*arg.max_step_xi,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    sampled=sol.sol(s0+x*slope).T
    full=sampled[:,:11];reduced=sampled[:,11:21];phi=sampled[:,21:121].reshape(-1,10,10);q=sampled[:,121:131]
    phi_ref=ref['P_full11D'][:,:10,:10];q_ref=ref['P_full11D'][:,:10,10]
    errphi=np.linalg.norm(phi-phi_ref,axis=(1,2))/np.maximum(np.linalg.norm(phi_ref,axis=(1,2)),1e-12)
    errq=np.linalg.norm(q-q_ref,axis=1)/np.maximum(np.linalg.norm(q_ref,axis=1),1e-8)
    state_err=np.linalg.norm(reduced-full[:,:10],axis=1)/np.maximum(np.linalg.norm(full[:,:10],axis=1),1)
    phase_err=np.max(np.abs(full[:,10]-(phase0+omega*x*slope)))
    sample_idx=sorted(set([0,40,80,160,240,320,400,480,560,640]))
    structural=[];rows=[]
    for i in sample_idx:
        yy=full[i];f=rhs(yy);J=jac(yy);a=yy.copy();a[10]+=1e-6
        diff=rhs(a)-f
        struct={'xi':float(x[i]),'rhs_clock_error':float(abs(f[10]-omega)),
          'jac_clock_max':float(np.max(abs(J[10,:]))),
          'clock_rhs_shift_change':float(abs(diff[10])),
          'J_z_phase_column_norm':float(np.linalg.norm(J[:10,10])),
          'phase_shift_linearization_residual':float(np.linalg.norm(diff[:10]-J[:10,10]*1e-6))}
        structural.append(struct)
    for i in range(len(x)):
        rows.append({'xi':float(x[i]),'full_vs_reduced_z_relative_error':float(state_err[i]),
          'phi10_against_archived_relative_error':float(errphi[i]),
          'phase_initial_variation_q_relative_error':float(errq[i]),
          'moving_phase_leaf_residual':float(full[i,10]-phase0-omega*x[i]*slope)})
    save_csv(ROOT/'results/ode_phase_replay_profile.csv',rows)
    save_csv(ROOT/'results/ode_phase_structural_probes.csv',structural)
    outcome={'version':'v18.11','replay_success':True,'original_ode_pinned_sources':files,
      'archived_kick203_state_used':True,'no_400_orbit_reintegration':True,
      'xi_samples':len(x),'source_phase0':phase0,'phase_rate_ode_clock':float(omega),
      'max_full_vs_reduced_10d_state_rel_error':float(np.max(state_err)),
      'max_reduced_10d_phi_rel_error_against_archived_11d':float(np.max(errphi)),
      'median_reduced_10d_phi_rel_error_against_archived_11d':float(np.median(errphi)),
      'max_driven_phase_column_q_rel_error_against_archived_11d':float(np.max(errq)),
      'median_driven_phase_column_q_rel_error_against_archived_11d':float(np.median(errq)),
      'max_moving_phase_clock_residual':float(phase_err),
      'max_original_rhs_phase_deviation_from_2pi':float(max(z['rhs_clock_error'] for z in structural)),
      'max_original_jac_phase_row':float(max(z['jac_clock_max'] for z in structural)),
      'min_original_z_phase_coupling_norm':float(min(z['J_z_phase_column_norm'] for z in structural)),
      'max_original_phase_shift_linearization_residual':float(max(z['phase_shift_linearization_residual'] for z in structural)),
      'phase_fibre_structural_status':'NUMERICAL_ODE_VERIFIED_NOT_MECHANICAL_CONSTRAINT',
      'true_material_constraint_or_metric_claim':False,
      'note':'Varying initial phase drives z and q, so a 10D autonomous vector field cannot be inferred. Phase-clock leaf is time dependent.'}
    # Two independent clocks and cocycle checks must actually agree; fail rather than silently returning success.
    if phase_err>1e-7 or np.max(state_err)>1e-7 or np.max(errphi)>3e-4 or np.max(errq)>3e-4 or outcome['max_original_jac_phase_row']>1e-10:
        raise RuntimeError('Independent ODE reduction regression failed: '+json.dumps(outcome,indent=2))
    (ROOT/'results/ode_phase_replay_summary.json').write_text(json.dumps(outcome,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'replay_success':True,'max_red_state_error':outcome['max_full_vs_reduced_10d_state_rel_error'],
     'max_phi10_error':outcome['max_reduced_10d_phi_rel_error_against_archived_11d'],
     'max_q_error':outcome['max_driven_phase_column_q_rel_error_against_archived_11d']},indent=2))
if __name__=='__main__':main()
