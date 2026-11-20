"""Reproducible frozen real-Q550 Jacobian jet, noise floor, and observer gauge audit."""
from pathlib import Path
import sys,json,csv,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
from checkpoint_io_v18_21 import read
from frozen_preflight import observers
from jet_math import derivatives,motion_matrices,observability,numerical_rank
from noise_gauge_math import oracle_stats,gauge_test

def writecsv(filename,rows):
    if not rows:return
    with (ROOT/'results'/filename).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    D,prov,prior,prior_csv=read();fine=D['tight'];Hs=observers(fine)
    xi=fine['xi_nodes'];dx=float(np.median(np.diff(xi)))
    if max(abs(np.diff(xi)-dx))>1e-10:raise RuntimeError('Nonuniform frozen grid')
    rows=[];grows=[];noise=[];compare=[];windows=[];covariances=[]
    for x in [2.,3.,4.,5.,6.,7.]:
        j=int(round(x/dx))
        mats={}
        for tol in ['standard','tight']:
            d=D[tol];arr=float(d['slope'])*d['J_full11D'][:,:10,:10]
            mats[tol]=motion_matrices(*derivatives(arr,j,dx))
        for name,H in Hs.items():
            for tau in [.03,.05,.1,.2,.5]:
                Ot=observability(H,mats['tight'],tau)[-1]
                Os=observability(H,mats['standard'],tau)[-1]
                s=np.linalg.svd(Ot,compute_uv=False)
                rank,_=numerical_rank(Ot)
                sr=np.linalg.svd(Os,compute_uv=False)
                rows.append(dict(xi=x,observer=name,tau=tau,order=3,rank=rank,sigma_min=float(s[-1]),sigma_max=float(s[0]),
                    condition=float(s[0]/s[-1]),sigma_min_standard=float(sr[-1]),
                    relative_smin_drift=float(abs(sr[-1]-s[-1])/max(s[-1],1e-15))))
                if x==4.:
                    if tau in [.05,.1,.2]:
                        for item in oracle_stats(Ot,[1e-6,1e-4,1e-3],seed=1821+len(noise)):
                            noise.append(dict(observer=name,tau=tau,xi=x,**item))
                    if tau==.1:
                        t=gauge_test(Ot,H.shape[0],3,seed=1821+len(grows))
                        grows.append(dict(observer=name,tau=tau,xi=x,**t))
            if x==4.:
                # Direct time-sample observability. No numerical derivatives are measured.
                P=fine['P_full11D'][:,:10,:10]
                P0inv=np.linalg.inv(P[j]);Pstandard=D['standard']['P_full11D'][:,:10,:10]
                for width in [.05,.1,.2,.5,1.]:
                    offsets=np.linspace(-width,width,9)
                    Ot=np.vstack([H@(P[int(round((4+z)/dx))]@P0inv) for z in offsets])
                    ss=np.linalg.svd(Ot,compute_uv=False)
                    Os=np.vstack([H@(Pstandard[int(round((4+z)/dx))]@
                                   np.linalg.inv(Pstandard[j])) for z in offsets])
                    st=np.linalg.svd(Os,compute_uv=False)
                    windows.append(dict(observer=name,center_xi=x,half_window=width,samples=9,
                       sigma_min=float(ss[-1]),condition=float(ss[0]/ss[-1]),
                       sigma_min_standard=float(st[-1]),
                       weak_unit_state_snr_at_noise_1e4=float(ss[-1]/1e-4)))
                # Scaling jet derivatives AND their uncertainty cannot improve information.
                Ounscaled=observability(H,mats['tight'],tau=1.)[-1]
                m=H.shape[0]
                physical_noise=np.array([1e-4,3e-4,2e-3,1e-2])
                baseline_weights=np.repeat(1./physical_noise,m)
                sing0=np.linalg.svd(baseline_weights[:,None]*Ounscaled,compute_uv=False)
                for t in [.03,.05,.1,.2,.5]:
                    row_scale=np.repeat(t**np.arange(4),m)
                    scaled=row_scale[:,None]*Ounscaled
                    scaled_stdev=row_scale*np.repeat(physical_noise,m)
                    sing=np.linalg.svd(scaled/scaled_stdev[:,None],compute_uv=False)
                    covariances.append(dict(observer=name,tau=t,
                        covariance_aware_sigma_min=float(sing[-1]),
                        covariance_aware_condition=float(sing[0]/sing[-1]),
                        relative_information_singular_error=float(np.linalg.norm(sing-sing0)/np.linalg.norm(sing0))))
            # Verify independent timestep stencils at same xi, tau 0.1
            for jump in [1,2,4]:
                d=D['tight'];arr=float(d['slope'])*d['J_full11D'][:,:10,:10]
                Mm=motion_matrices(*derivatives(arr,j,dx,jump=jump))
                ss=np.linalg.svd(observability(H,Mm,.1)[-1],compute_uv=False)
                compare.append(dict(observer=name,xi=x,stencil_jump=jump,sigma_min=float(ss[-1]),condition=float(ss[0]/ss[-1])))
    writecsv('jet_noise_singular_spectrum.csv',rows)
    writecsv('noise_recovery_trials.csv',noise)
    writecsv('coordinate_gauge_invariance.csv',grows)
    writecsv('jet_stencil_noise_sensitivity.csv',compare)
    writecsv('direct_time_window_observability.csv',windows)
    writecsv('noise_covariance_scale_invariance.csv',covariances)
    # Include sensitivity analyses for real nonlinear checks archived in predecessor.
    (ROOT/'results/nonlinear_v18_20_reference.csv').write_text(prior_csv)
    out={'version':'v18.21','status':'FROZEN_JET_NOISE_AND_GAUGE_PREFLIGHT',
         'source_status':prior['status'],'provenance':prov,'xi_grid':[2,3,4,5,6,7],
         'taus':[.03,.05,.1,.2,.5],
         'noise_model':'Independent isotropic Gaussian perturbations on scaled jet stack, NOT raw instrument samples.',
         'unit_state_reconstruction':'Weakest right-singular unit direction in nominal Euclidean 10D gauge',
         'gauge_model':'Constant invertible material/output transformations with covariance-aware output whitening.',
         'rank_tolerance':1e-8,'trials_per_setting':256,
         'gauge_max_physical_state_error':max(x['physical_state_error'] for x in grows),
         'gauge_max_physical_prediction_error':max(x['physical_prediction_error'] for x in grows),
         'rank_full':all(x['rank']==10 for x in rows),
         'max_covariance_consistent_rescaling_error':max(x['relative_information_singular_error'] for x in covariances),
         'at_xi4_direct_time_window_examples':{str(x['observer'])+'_'+str(x['half_window']):x for x in windows if x['half_window'] in [.05,.5]},
         'most_ill_conditioned':max(rows,key=lambda t:t['condition']),
         'at_xi4_tau0_1':{r['observer']:r for r in rows if r['xi']==4. and r['tau']==.1},
         'claim_limits':['Numerical local tangent rank only; no interval or global certificate.',
                         'Noise model assumes directly measured derivatives with equal noise in scaled jet coordinates.',
                         'Real finite-difference jet noise is correlated and can be amplified; not represented by this optimistic test.',
                         'Condition numbers are gauge- and unit-dependent without physical metrics/covariances.',
                         'No physical material constraint manifold is certified.',
                         'Memory/hysteresis is not diagnosed by rank or conditioning alone.']}
    (ROOT/'results/frozen_noise_gauge_summary.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    if not out['rank_full'] or out['gauge_max_physical_prediction_error']>1e-7 or out['max_covariance_consistent_rescaling_error']>1e-10:
        raise RuntimeError('Unexpected gauge or rank property')
    print(json.dumps({'status':out['status'],'at_xi4_tau0_1':{k:{'smin':v['sigma_min'],'cond':v['condition']} for k,v in out['at_xi4_tau0_1'].items()},
        'max_gauge_physical_state_error':out['gauge_max_physical_state_error'],'server_nonlinear_v1820':out['source_status']},indent=2))
    return out
if __name__=='__main__':main()
