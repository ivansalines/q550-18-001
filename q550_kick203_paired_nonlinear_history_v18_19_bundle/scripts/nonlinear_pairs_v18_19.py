"""Direct ORIGINAL Q550 ODE: paired histories with matched reduced state/first jet.

Off-constraint test. Equal projected present is constructed, not constitutive hysteresis.
No phase re-registration adjustment is made beyond fixing identical theta at match.
"""
from pathlib import Path
import sys,argparse,json,csv,platform
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
import q550_state_impulse_response_v13 as original
from checkpoint_io_v18_19 import read
from frozen_preflight import observers
from pair_math import chosen_seed


def integrate(f,x0,x1,y0,step,rtol,atol,grid):
    sol=solve_ivp(f,(float(x0),float(x1)),np.asarray(y0,float),method='DOP853',
        rtol=rtol,atol=atol,max_step=step,t_eval=np.asarray(grid,float))
    if not sol.success:raise RuntimeError('Nonlinear ODE integration failed: '+sol.message)
    return sol.y.T,{'nfev':sol.nfev,'message':sol.message}

def write_csv(path,rows):
    if not rows:return
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--base-root',type=Path,default=original.DEFAULT_BASE)
    ap.add_argument('--eps',nargs='+',type=float,default=[1e-6,3e-6,1e-5])
    ap.add_argument('--max-step-xi',type=float,default=.04)
    ap.add_argument('--rtol',type=float,default=4e-12)
    ap.add_argument('--atol',type=float,default=5e-13)
    a=ap.parse_args()
    if any(t<=0 or t>1e-3 for t in a.eps):raise ValueError('Epsilon range must be positive <= 1e-3')
    data,provenance=read(); d=data['tight']
    source_files={}
    for kind in ('audit','support','tube'):
        fname,pin=original.DEPENDENCIES[kind];path=a.base_root/fname
        if not path.is_file():raise FileNotFoundError(f'Original Q550 {kind} missing: {path}')
        got=original.hash_file(path)
        if got!=pin:raise RuntimeError(f'PINNED ORIGINAL ODE SHA MISMATCH: {kind} {got} != {pin}')
        source_files[kind]={'path':str(path),'sha256':got}
    audit=original.load_module(a.base_root/original.DEPENDENCIES['audit'][0],'v1819_pinned_audit')
    support=original.load_module(a.base_root/original.DEPENDENCIES['support'][0],'v1819_pinned_support')
    params=audit.ModelParameters(gamma_rho=.3)
    kappa0=(55./6.)*(10./550.)
    slope=float(d['slope']);s0=float(d['s0']);y0=np.asarray(d['source_state_11D'],float)
    def rhs(y):
        z=np.asarray(y,float)
        return np.asarray(support.support_rhs_float(z,kappa0=kappa0,params=params,audit=audit),float)
    def f(xi,y):return slope*rhs(y)
    xi_m=4.;grid_xi=np.arange(1.,8.00001,.25)
    # Forward baseline is independently integrated from the v18.6 original kick state.
    base_nodes=np.array([0.,1.,2.,3.,4.])
    base,_=integrate(f,0.,xi_m,y0,a.max_step_xi,a.rtol,a.atol,base_nodes)
    ymatch=base[-1].copy()
    i_m=int(round(xi_m/.0125));Jarchive=d['J_full11D'][i_m]
    # Independently differentiated original RHS at the actual matching state.
    J=np.zeros((11,11));diff_sizes=[]
    for i in range(11):
        st=1e-6*max(1.,abs(ymatch[i]))
        ei=np.zeros(11);ei[i]=st
        J[:,i]=(rhs(ymatch+ei)-rhs(ymatch-ei))/(2*st)
        diff_sizes.append(st)
    jac_ref_rel=float(np.linalg.norm(J-Jarchive)/max(1.,np.linalg.norm(Jarchive)))
    if jac_ref_rel>2e-4:raise RuntimeError('Independent ODE Jacobian differs from archived baseline at xi=4: '+str(jac_ref_rel))
    if np.max(np.abs(J[10]))>5e-8:raise RuntimeError('Unexpected phase Jacobian lower-row coupling')
    Hs=observers(d)
    modes=[('material3','jet',False),('material3','jet',True),('material3','state_only',False),
           ('raw4','jet',False),('raw4','state_only',False),('raw6','state_only',False)]
    history_grid=np.arange(4.,.99999,-.25)
    future_grid=np.arange(4.,8.000001,.25)
    Pm=d['P_full11D'][i_m];P=d['P_full11D']; xi_nodes=d['xi_nodes'];
    report=[];rows=[];max_clock=0.;max_state=0.;max_jet=float('nan');max_tangent=0.
    for obs,match,second in modes:
        H=Hs[obs];seed,seed_diag=chosen_seed(H,J[:10,:10],match,second)
        label=obs+'_'+match+('_b' if second else '_a')
        seed11=np.r_[seed,0.]
        predicted_initial=np.linalg.solve(Pm,seed11)
        for eps in a.eps:
            plus=ymatch+eps*seed11;minus=ymatch-eps*seed11
            # Exactly same phase and same *linear* reduced state at meeting time.
            state_match=float(np.linalg.norm(H@(plus[:10]-minus[:10])))
            phase_match=float(abs(plus[10]-minus[10]))
            deriv_match=float(np.linalg.norm(H@(f(xi_m,plus)[:10]-f(xi_m,minus)[:10])))
            if state_match>2e-11 or phase_match>1e-12:raise RuntimeError('Matching initial observations failed')
            # Separate original-RHS solve_ivp calls for each trajectory: no shared
            # adaptive step controller, unlike a concatenated 22D solve.
            bp,bi_p=integrate(f,xi_m,1.,plus,a.max_step_xi,a.rtol,a.atol,history_grid)
            bn,bi_n=integrate(f,xi_m,1.,minus,a.max_step_xi,a.rtol,a.atol,history_grid)
            fp,fi_p=integrate(f,xi_m,8.,plus,a.max_step_xi,a.rtol,a.atol,future_grid)
            fn,fi_n=integrate(f,xi_m,8.,minus,a.max_step_xi,a.rtol,a.atol,future_grid)
            allp=np.concatenate((bp[::-1][:-1],fp),axis=0)
            alln=np.concatenate((bn[::-1][:-1],fn),axis=0)
            both=np.stack((allp,alln),axis=1)
            back_info={'nfev':bi_p['nfev']+bi_n['nfev']}
            fw_info={'nfev':fi_p['nfev']+fi_n['nfev']}
            phase_diffs=both[:,0,10]-both[:,1,10]
            clock_drift=float(np.max(abs(phase_diffs)))
            max_clock=max(max_clock,clock_drift);max_state=max(max_state,state_match)
            traj_err=[];future_dist=[];past_dist=[]
            for j,xi in enumerate(grid_xi):
                yp,yn=both[j]
                diff=yp-yn
                xi_idx=int(round(xi/.0125))
                expected=P[xi_idx]@predicted_initial
                odd=diff/(2*eps)
                tangent_rel=float(np.linalg.norm(odd-expected)/max(1.,np.linalg.norm(expected)))
                max_tangent=max(max_tangent,tangent_rel)
                reduced_sep=float(np.linalg.norm(H@diff[:10]))
                full_sep=float(np.linalg.norm(diff[:10]))
                rows.append({'observer':obs,'match':match,'second_null':int(second),'label':label,
                     'epsilon':eps,'xi':float(xi),'relative_xi':float(xi-xi_m),
                     'reduced_separation':reduced_sep,'full10D_separation':full_sep,
                     'scaled_reduced_separation':reduced_sep/(2*eps),
                     'tangent_response_relative_error':tangent_rel,
                     'phase_separation':float(diff[10])})
                if xi>xi_m+1e-8:future_dist.append((xi,reduced_sep))
                if xi<xi_m-1e-8:past_dist.append((xi,reduced_sep))
            future_at_1=float(np.linalg.norm(H@(both[np.argmin(abs(grid_xi-5.)),0,:10]-both[np.argmin(abs(grid_xi-5.)),1,:10])))
            future_at_4=float(np.linalg.norm(H@(both[-1,0,:10]-both[-1,1,:10])))
            past_at_minus_1=float(np.linalg.norm(H@(both[np.argmin(abs(grid_xi-3.)),0,:10]-both[np.argmin(abs(grid_xi-3.)),1,:10])))
            report.append({'label':label,'observer':obs,'match':match,'second_null':bool(second),
                'epsilon':eps,'hidden_seed':seed.tolist(), 'seed_diagnostics':seed_diag,
                'meeting_reduced_state_gap':state_match,'meeting_reduced_derivative_gap':deriv_match,
                'meeting_reduced_derivative_gap_over_2eps':deriv_match/(2*eps),
                'past_reduced_gap_at_xi3':past_at_minus_1,
                'future_reduced_gap_at_xi5':future_at_1,
                'future_reduced_gap_at_xi8':future_at_4,
                'peak_future_scaled_reduced_gap':max(v for x,v in future_dist)/(2*eps),
                'backward_nfev':back_info['nfev'],'forward_nfev':fw_info['nfev'],
                'max_phase_separation':clock_drift})
    # Same-full-state control: identical deterministic initial state, identical solution.
    control_yp,_=integrate(f,xi_m,8.,ymatch,a.max_step_xi,a.rtol,a.atol,future_grid)
    control_yn,_=integrate(f,xi_m,8.,ymatch,a.max_step_xi,a.rtol,a.atol,future_grid)
    full_state_control_gap=float(np.max(abs(control_yp-control_yn)))
    if full_state_control_gap>1e-12:raise RuntimeError('Deterministic full-state control is not identical')
    # Nonlinear phase clock is independent on each trajectory (checked over all samples).
    if max_clock>3e-8:raise RuntimeError('Unexpected phase mismatch')
    summaries={'version':'v18.19','status':'NONLINEAR_PAIRED_REDUCED_OBSERVER_DIAGNOSTIC_COMPLETED',
      'constraint_status':'OFF_CONSTRAINT_NO_PHYSICAL_INVARIANTS_SUPPLIED',
      'source_files':source_files,'provenance':provenance,
      'ode':{'s0':s0,'ds_dxi':slope,'kappa0':kappa0,'xi_meeting':xi_m,'xi_past_start':1.,'xi_future_end':8.,
             'rtol':a.rtol,'atol':a.atol,'max_step_xi':a.max_step_xi,
             'independent_original_rhs_jacobian_relative_error':jac_ref_rel,
             'original_rhs_max_phase_jacobian_row':float(np.max(abs(J[10])))},
      'pair_count':len(report),'modes':len(modes),'amplitudes':list(a.eps),
      'max_reduced_state_match_error':max_state,'max_phase_difference':max_clock,
      'max_nonlinear_vs_tangent_relative_error':max_tangent,
      'full_state_equal_future_max_gap':full_state_control_gap,
      'pair_summary':report,
      'interpretation_gate':{
        'reduced_nonclosure_can_be_reported_if': 'Same reduced state (+ first jet for jet cases), distinct hidden state, and future reduced responses differ above numerical error',
        'not_established':'Constitutive physical memory, unique material metric, or a torsional flip',
        'why':'For deterministic ODE, hidden present-state components suffice to explain differences. Pairs are created at common xi, then integrated backward/forward; no independent constitutive state law is inferred.',
        'role_of_history':'Past differences corroborate separate histories, but are not a causal independent mechanism for future differences once the full state is specified.'}}
    out=ROOT/'results';out.mkdir(exist_ok=True)
    write_csv(out/'paired_trajectories.csv',rows)
    write_csv(out/'pair_summary.csv',[{k:v for k,v in d.items() if k not in ('hidden_seed','seed_diagnostics')} for d in report])
    (out/'nonlinear_pair_summary.json').write_text(json.dumps(summaries,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':summaries['status'],'pair_count':len(report),
      'max_phase_gap':max_clock,'independent_jacobian_error':jac_ref_rel,
      'max_nonlinear_tangent_error':max_tangent,'control_gap':full_state_control_gap},indent=2))

if __name__=='__main__':main()
