"""Independent original-ODE paired nonlinear validation of v18.20 frozen jet predictions.

Uses already pinned original Q550 sources and original kick203 state; NO 400-orbit replay.
Matching through k means jets up to order k vanish ONLY TO LINEAR ORDER in
material perturbation, not exactly at finite epsilon. This is an off-constraint
observability audit, not proof of physical hereditary memory.
"""
from pathlib import Path
import argparse,json,csv,sys
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
from checkpoint_io_v18_20 import read
from frozen_preflight import observers
import q550_state_impulse_response_v13 as original
from jet_math import derivatives,motion_matrices,scaled_taylor

def integrate(f,start,stop,y0,evalpoints,rtol,atol,max_step):
    sol=solve_ivp(f,(float(start),float(stop)),np.asarray(y0,float),method='DOP853',
      t_eval=np.asarray(evalpoints,float),rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success:raise RuntimeError('Q550 original nonlinear ODE failure: '+sol.message)
    return sol.y.T,sol.nfev

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--base-root',type=Path,default=original.DEFAULT_BASE)
    ap.add_argument('--eps',nargs='+',type=float,default=[3e-6,1e-5])
    ap.add_argument('--max-step-xi',type=float,default=.02)
    ap.add_argument('--rtol',type=float,default=4e-12)
    ap.add_argument('--atol',type=float,default=5e-13)
    args=ap.parse_args()
    if any(e<=0 or e>1e-3 for e in args.eps):raise ValueError('invalid epsilon')
    dct,prov,prior=read();d=dct['tight']
    expected_sources={}
    for kind in ('audit','support','tube'):
        filename,sha=original.DEPENDENCIES[kind];path=args.base_root/filename
        if not path.is_file():raise FileNotFoundError('Pinned original ODE dependency missing: '+str(path))
        digest=original.hash_file(path)
        if digest!=sha:raise RuntimeError('Q550 original model PIN mismatch: '+str(path))
        expected_sources[kind]={'file':str(path),'sha256':digest}
    audit=original.load_module(args.base_root/original.DEPENDENCIES['audit'][0],'v1820_orig_audit')
    support=original.load_module(args.base_root/original.DEPENDENCIES['support'][0],'v1820_orig_support')
    params=audit.ModelParameters(gamma_rho=.3)
    kappa0=(55./6.)*(10./550.)
    slope=float(d['slope']);y0=np.asarray(d['source_state_11D'],float)
    def rhs(y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=kappa0,params=params,audit=audit),float)
    def f(xi,y):return slope*rhs(y)
    base,nfev=integrate(f,0.,4.,y0,[4.],args.rtol,args.atol,args.max_step_xi)
    state=base[-1].copy()
    xi_nodes=d['xi_nodes'];grid_step=float(np.median(np.diff(xi_nodes)));at=int(round(4/grid_step));
    Agrid=slope*d['J_full11D'][:,:10,:10];A,D1,D2=derivatives(Agrid,at,grid_step,1)
    M=motion_matrices(A,D1,D2);Hs=observers(d)
    # Original RHS independent central-difference check of material Jacobian.
    J=np.empty((11,11))
    for j in range(11):
        hh=1e-6*max(1.,abs(state[j]));dv=np.eye(11)[j]*hh
        J[:,j]=(rhs(state+dv)-rhs(state-dv))/(2*hh)
    jac_err=float(np.linalg.norm(J-d['J_full11D'][at])/max(1.,np.linalg.norm(d['J_full11D'][at])))
    if jac_err>2e-4:raise RuntimeError('Pinned original RHS Jacobian fails checkpoint match')
    if np.linalg.norm(J[10])>5e-8:raise RuntimeError('Unexpected phase Jacobian row')
    seeds=json.loads((ROOT/'results/matched_hidden_seeds.json').read_text())
    P=d['P_full11D'][:,:10,:10];Pm_inv=np.linalg.inv(P[at]);rows=[];seed_summaries=[]
    deltas=[-.1,-.05,-.025,0.,.025,.05,.1,.2,.5]
    forward=[4+x for x in deltas if x>=0];backward=[4+x for x in deltas if x<=0][::-1]
    for seed_info in seeds:
        label=seed_info['name'];H=Hs[seed_info['observer']];v=np.asarray(seed_info['seed']);order=seed_info['matched_through_order']
        for eps in args.eps:
            yp=state+eps*np.r_[v,0.];yn=state-eps*np.r_[v,0.]
            plusF,_=integrate(f,4.,forward[-1],yp,forward,args.rtol,args.atol,args.max_step_xi)
            minusF,_=integrate(f,4.,forward[-1],yn,forward,args.rtol,args.atol,args.max_step_xi)
            plusB,_=integrate(f,4.,backward[-1],yp,backward,args.rtol,args.atol,args.max_step_xi)
            minusB,_=integrate(f,4.,backward[-1],yn,backward,args.rtol,args.atol,args.max_step_xi)
            lookup={}
            for ts,pp,nn in [(forward,plusF,minusF),(backward,plusB,minusB)]:
                for x,a,b in zip(ts,pp,nn):lookup[round(x-4.,9)]=(a,b)
            rows_for_seed=[]
            for delta in deltas:
                pp,nn=lookup[round(delta,9)]
                diff=(pp-nn)/(2*eps)
                measured=H@diff[:10]
                j=int(round((4+delta)/grid_step))
                predicted=H@P[j]@Pm_inv@v
                taylor=scaled_taylor(H,M,v,delta)
                tangent_error=float(np.linalg.norm(measured-predicted))
                row=dict(label=label,observer=seed_info['observer'],matched_through_order=order,
                    epsilon=eps,delta_xi=delta,
                    scaled_reduced_separation=float(np.linalg.norm(measured)),
                    predicted_tangent_separation=float(np.linalg.norm(predicted)),
                    predicted_taylor_separation=float(np.linalg.norm(taylor)),
                    abs_nonlinear_tangent_error=tangent_error,
                    abs_tangent_taylor_error=float(np.linalg.norm(predicted-taylor)),
                    phase_gap=float(abs(pp[10]-nn[10])),
                    measured_components=json.dumps(measured.tolist()),
                    taylor_components=json.dumps(taylor.tolist()))
                rows.append(row);rows_for_seed.append(row)
            seed_summaries.append({'label':label,'epsilon':eps,'matched_order':order,
                'maximum_nonlinear_tangent_error':max(q['abs_nonlinear_tangent_error'] for q in rows_for_seed),
                'maximum_phase_gap':max(q['phase_gap'] for q in rows_for_seed),
                'nonlinear_vs_tangent_at_delta_0_1':next(q['abs_nonlinear_tangent_error'] for q in rows_for_seed if q['delta_xi']==.1),
                'scaled_gap_at_0_1':next(q['scaled_reduced_separation'] for q in rows_for_seed if q['delta_xi']==.1)})
    # Deterministic identical-full-state control
    c1,_=integrate(f,4.,4.5,state,[4.5],args.rtol,args.atol,args.max_step_xi)
    c2,_=integrate(f,4.,4.5,state,[4.5],args.rtol,args.atol,args.max_step_xi)
    eq=float(np.max(abs(c1-c2)))
    if eq>1e-12:raise RuntimeError('Deterministic identical-full-state control failed')
    out=ROOT/'results'
    with (out/'nonlinear_jet_validation.csv').open('w',newline='') as f:
        names=[k for k in rows[0] if k not in ('measured_components','taylor_components')]
        wr=csv.DictWriter(f,fieldnames=names);wr.writeheader();wr.writerows([{k:r[k] for k in names} for r in rows])
    report={'version':'v18.20','status':'ORIGINAL_Q550_NONLINEAR_JET_VALIDATION_COMPLETED',
      'pinned_sources':expected_sources,'provenance':prov,'original_jacobian_relative_error':jac_err,
      'original_max_phase_jacobian_row':float(np.max(abs(J[10]))),
      'eps':args.eps,'rtol':args.rtol,'atol':args.atol,'max_step_xi':args.max_step_xi,
      'probe_pair_count':len(seed_summaries),'identical_full_state_control_max_gap':eq,
      'maximum_nonlinear_tangent_abs_error':max(q['maximum_nonlinear_tangent_error'] for q in seed_summaries),
      'max_phase_gap':max(q['maximum_phase_gap'] for q in seed_summaries),
      'matched_probe_summaries':seed_summaries,
      'interpretation_limit':'Off-constraint jets for a single full-state deterministic ODE and diagnostic projections. Full-rank local tangent observability neither certifies nonlinear global injectivity nor physical hereditary memory.'}
    (out/'nonlinear_jet_summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],'pair_count':report['probe_pair_count'],
      'original_jacobian_rel':jac_err,'max_phase_gap':report['max_phase_gap'],
      'max_nonlinear_tangent_abs_err':report['maximum_nonlinear_tangent_abs_error']},indent=2))
if __name__=='__main__':main()
