"""Fresh nonlinear raw trajectory test on original pinned Q550 ODE.

Six perturbed trajectories + one baseline at xi4, integrated backward and forward;
same matched hidden seeds as prior v18.20, with independent raw data sampling.
Synthetic noise is added ONLY after ODE integration and normalized paired difference.
"""
from pathlib import Path
import sys,zipfile,json,csv,argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
from checkpoint_io_v18_22 import read
from checkpoint_io_v18_21 import PREFIX as SRC20
from frozen_preflight import observers
from raw_jet_math import basis_times,continuous_observation,jet_compression,direct_estimator,jet_estimator,compressed_observation
from jet_math import derivatives,motion_matrices
import q550_state_impulse_response_v13 as original


def evolve(f,a,b,state,grid):
    sol=solve_ivp(f,(a,b),state,t_eval=grid,method='DOP853',rtol=4e-12,atol=5e-13,max_step=.02)
    if not sol.success:raise RuntimeError('Original nonlinear ODE integration failed: '+sol.message)
    return sol.y.T


def main():
    p=argparse.ArgumentParser();p.add_argument('--base-root',type=Path,default=original.DEFAULT_BASE)
    args=p.parse_args()
    sources={}
    for kind,(filename,digest) in original.DEPENDENCIES.items():
        source=args.base_root/filename
        if not source.is_file():raise FileNotFoundError('Pinned source absent '+str(source))
        if original.hash_file(source)!=digest:raise RuntimeError('Pinned ODE dependency checksum mismatch: '+str(source))
        sources[kind]={'path':str(source),'sha256':digest}
    D,prov,_,source_proof=read();d=D['tight']; Hs=observers(d)
    with zipfile.ZipFile(ROOT/'.runtime/v18_20_original_server_sendback.zip') as z:
        probes=json.loads(z.read(SRC20+'results/matched_hidden_seeds.json'))
    probes=[a for a in probes if a['name'] in ('material3_jet2','raw4_jet1','raw6_jet1')]
    if len(probes)!=3:raise RuntimeError('Seed genealogy mismatch')
    audit=original.load_module(args.base_root/original.DEPENDENCIES['audit'][0],'v1822_orig_audit')
    support=original.load_module(args.base_root/original.DEPENDENCIES['support'][0],'v1822_orig_support')
    params=audit.ModelParameters(gamma_rho=.3); kappa0=(55./6.)*(10./550.)
    slope=float(d['slope']);y0=np.asarray(d['source_state_11D'],float)
    def f(x,y):return slope*np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=kappa0,params=params,audit=audit),float)
    base=evolve(f,0,4,y0,[4.])[-1]
    xi=np.asarray(d['xi_nodes'],float);dx=float(np.median(np.diff(xi)))
    j=int(np.argmin(abs(xi-4.)))
    P=np.asarray(d['P_full11D'][:,:10,:10],float)
    spline=CubicSpline(xi,P,axis=0)
    A=slope*np.asarray(d['J_full11D'][:,:10,:10],float)
    mats=motion_matrices(*derivatives(A,j,dx))
    widths=[.2,.5];epsilon=1e-5;sigma=1e-4;trials=64
    rng=np.random.default_rng(182219)
    # Exact 33 sample times per window; do NOT round/deduplicate to checkpoint nodes.
    times_by_width={width:4.+basis_times(width,33) for width in widths}
    t_all=np.array(sorted(set(float(t) for grid in times_by_width.values() for t in grid)))
    backward=t_all[t_all<=4][::-1]
    forward=t_all[t_all>=4]
    rows=[];recovery=[]
    for seed in probes:
        label=seed['name'];H=Hs[seed['observer']]
        v=np.asarray(seed['seed'],float);v/=np.linalg.norm(v)
        states=[]
        for sgn in [1.,-1.]:
            ini=base+sgn*epsilon*np.r_[v,0.]
            backward_states=evolve(f,4.,float(backward[-1]),ini,backward)
            forward_states=evolve(f,4.,float(forward[-1]),ini,forward)
            dictionary={float(t):y for t,y in zip(backward,backward_states)}
            dictionary.update({float(t):y for t,y in zip(forward,forward_states)})
            states.append(dictionary)
        for width in widths:
            t=times_by_width[width];offset=t-4
            raw=np.array([H@((states[0][float(z)]-states[1][float(z)])[:10]/(2*epsilon)) for z in t])
            phase=max(abs(states[0][float(z)][10]-states[1][float(z)][10]) for z in t)
            O=continuous_observation(spline,P[j],t,H)
            predicted=(O@v).reshape(-1,H.shape[0]);err=np.linalg.norm(raw-predicted,axis=1)
            rows.append(dict(label=label,observer=seed['observer'],half_window=width,epsilon=epsilon,
              sample_count=len(t),max_tangent_error=float(err.max()),
              rms_tangent_error=float(np.sqrt(np.mean(err**2))),phase_difference=float(phase)))
            for degree in [3,5]:
                B4,C,W,Oj,_=jet_compression(offset,width,degree,H.shape[0],[H@q for q in mats])
                pinv_d=direct_estimator(O);pinv_j=jet_estimator(Oj)
                rms_d=[];rms_j=[]
                for _ in range(trials):
                    y=(raw+rng.normal(0,sigma,size=raw.shape)).reshape(-1)
                    xd=pinv_d@y;xj=pinv_j@compressed_observation(y,H.shape[0],B4,W)
                    rms_d.append(float(np.linalg.norm(xd-v)))
                    rms_j.append(float(np.linalg.norm(xj-v)))
                recovery.append(dict(label=label,observer=seed['observer'],half_window=width,
                  poly_degree=degree,noise_sigma=sigma,trials=trials,
                  direct_rms=float(np.sqrt(np.mean(np.square(rms_d)))),
                  jet_rms=float(np.sqrt(np.mean(np.square(rms_j))))))
    with (ROOT/'results/original_nonlinear_raw_validation.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (ROOT/'results/original_nonlinear_noise_reconstruction.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(recovery[0]));w.writeheader();w.writerows(recovery)
    summary={'version':'v18.22','status':'ORIGINAL_Q550_RAW_SAMPLE_NONLINEAR_VALIDATION_COMPLETED',
      'pinned_original_ode_sources':sources,'source_prior_status':source_proof['status'],
      'matched_hidden_probes':len(probes),'nonlinear_paired_trajectories':len(probes),
      'nonlinear_integrations':2*len(probes)*2+1,
      'widths':widths,'epsilon':epsilon,'noise_sigma_normalized':sigma,
      'max_abs_nonlinear_tangent_error':max(a['max_tangent_error'] for a in rows),
      'max_phase_difference':max(a['phase_difference'] for a in rows),
      'limitation':'An independent original ODE integration of synthetic paired initial states, not measurements of physical instrumentation.'}
    (ROOT/'results/original_ode_raw_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if summary['max_abs_nonlinear_tangent_error']>1e-3 or summary['max_phase_difference']>1e-8:
        raise RuntimeError('Original nonlinear raw-sample spot check failed')
    print(json.dumps({'status':summary['status'],'max_tangent_error':summary['max_abs_nonlinear_tangent_error'],
      'max_phase_difference':summary['max_phase_difference']},indent=2))
if __name__=='__main__':main()
