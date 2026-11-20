"""Fresh independently integrated original-ODE check for selected v18.20 hidden seeds.
This validates the response data that motivate the v18.21 frozen noise stress test;
it does NOT turn synthetic measurement noise into experimental instrument noise.
"""
from pathlib import Path
import argparse,sys,json,csv,zipfile,numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
from checkpoint_io_v18_21 import read,PREFIX
from frozen_preflight import observers
import q550_state_impulse_response_v13 as original

def solve(f,lo,hi,y,ts):
    ans=solve_ivp(f,(lo,hi),y,method='DOP853',t_eval=ts,rtol=4e-12,atol=5e-13,max_step=.02)
    if not ans.success:raise RuntimeError('Original-Q550 nonlinear integration failure '+ans.message)
    return ans.y.T

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base-root',type=Path,default=original.DEFAULT_BASE)
    args=ap.parse_args()
    sources={}
    for kind,(filename,sha) in original.DEPENDENCIES.items():
        p=args.base_root/filename
        if not p.is_file():raise FileNotFoundError('Pinned ODE source missing: '+str(p))
        got=original.hash_file(p)
        if got!=sha:raise RuntimeError('Pinned ODE SHA256 mismatch '+str(p))
        sources[kind]=dict(path=str(p),sha256=sha)
    D,prov,prior,_=read();d=D['tight'];Hs=observers(d)
    with zipfile.ZipFile(ROOT/'inputs/v18_20_original_server_sendback.zip') as z:
        seeds=json.loads(z.read(PREFIX+'results/matched_hidden_seeds.json'))
    seeds=[s for s in seeds if s['name'] in ('material3_jet2','raw4_jet1','raw6_jet1')]
    if len(seeds)!=3:raise RuntimeError('Missing original v18.20 hidden probes')
    audit=original.load_module(args.base_root/original.DEPENDENCIES['audit'][0],'v1821_orig_audit')
    support=original.load_module(args.base_root/original.DEPENDENCIES['support'][0],'v1821_orig_support')
    params=audit.ModelParameters(gamma_rho=.3);kappa0=(55./6.)*(10./550.)
    slope=float(d['slope']);y0=np.asarray(d['source_state_11D'],float)
    def rhs(y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=kappa0,params=params,audit=audit),float)
    def f(xi,y):return slope*rhs(y)
    state=solve(f,0,4,y0,[4])[-1]; grid_step=float(np.median(np.diff(d['xi_nodes'])))
    j=int(round(4/grid_step));P=d['P_full11D'][:,:10,:10];pinv=np.linalg.inv(P[j])
    deltas=[0.,.05,.1,.2,.5];times=[4+x for x in deltas];out=[]
    for info in seeds:
        label=info['name'];H=Hs[info['observer']];v=np.asarray(info['seed'],float)
        for eps in (3e-6,1e-5):
            plus=solve(f,4,4.5,state+eps*np.r_[v,0.],times)
            minus=solve(f,4,4.5,state-eps*np.r_[v,0.],times)
            for delta,yp,ym in zip(deltas,plus,minus):
                observed=H@((yp-ym)[:10]/(2*eps))
                idx=int(round((4+delta)/grid_step))
                prediction=H@(P[idx]@pinv@v)
                out.append(dict(label=label,eps=eps,delta_xi=delta,
                    observed_norm=float(np.linalg.norm(observed)),
                    tangent_norm=float(np.linalg.norm(prediction)),
                    absolute_error=float(np.linalg.norm(observed-prediction)),
                    phase_difference=float(abs(yp[10]-ym[10]))))
    worst=max(x['absolute_error'] for x in out)
    if worst>1e-3 or max(x['phase_difference'] for x in out)>1e-8:
        raise RuntimeError('Original ODE nonlinear tangent validation failed')
    with (ROOT/'results/independent_original_ode_check.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(out[0]));wr.writeheader();wr.writerows(out)
    summary=dict(version='v18.21',status='ORIGINAL_Q550_INDEPENDENT_NONLINEAR_SPOT_CHECK_COMPLETED',
        pinned_ode_sources=sources,source_v18_20_status=prior['status'],
        probes=len(seeds),pair_count=6,epsilon=[3e-6,1e-5],samples_per_pair=len(deltas),
        max_abs_nonlinear_tangent_error=worst,
        max_phase_difference=max(x['phase_difference'] for x in out),
        limitation='This checks physical ODE response of diagnostic directions, not physical noise floor, constraints, or global observability.')
    (ROOT/'results/original_ode_check_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(dict(status=summary['status'],pairs=summary['pair_count'],max_abs_error=worst),indent=2))
if __name__=='__main__':main()
