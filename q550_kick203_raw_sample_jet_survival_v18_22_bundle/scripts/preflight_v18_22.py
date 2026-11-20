"""Frozen Q550 tangent raw-observation experiment: matched-noise direct vs reconstructed jets.
All numerical data come from the authenticated real v18.21 original-server sendback.
NO fresh original ODE calculation is claimed in this preflight.
"""
from pathlib import Path
import csv,json,sys
import numpy as np
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
from checkpoint_io_v18_22 import read
from frozen_preflight import observers
from jet_math import derivatives,motion_matrices
from raw_jet_math import *

CONFIG={'center':4.0,'sample_count':33,'halfwidths':[.05,.1,.2,.5,1.0],
        'noise_levels':[1e-5,1e-4,1e-3], 'repeats':128,'poly_degrees':[3,5],
        'source':'original-v18.21 server / v18.20 hidden seeds, original tangent propagator',
        'units':'Noise sigma applied to raw perturbation response normalized by 2*epsilon; no physical units',
        'seed':1822}

def write_csv(fname,rows):
    if not rows:return
    with (ROOT/'results'/fname).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)

def reconstruct():
    data,prov,oldproof,newproof=read()
    d=data['tight'];ds=data['standard'];Hs=observers(d)
    xi=np.asarray(d['xi_nodes'],float)
    dx=float(np.median(np.diff(xi)))
    i=int(np.argmin(np.abs(xi-4.0)))
    if abs(xi[i]-4.0)>1e-9:raise RuntimeError('Reference xi 4 checkpoint unavailable')
    P=np.asarray(d['P_full11D'][:,:10,:10],float)
    Pst=np.asarray(ds['P_full11D'][:,:10,:10],float)
    spline=CubicSpline(xi,P,axis=0)
    spline_standard=CubicSpline(xi,Pst,axis=0)
    A=float(d['slope'])*np.asarray(d['J_full11D'][:,:10,:10],float)
    Mm=motion_matrices(*derivatives(A,i,dx))
    # Matched, originally nonlinear-validated hidden directions from server v18.20.
    from checkpoint_io_v18_21 import PREFIX as SRC20
    import zipfile
    with zipfile.ZipFile(ROOT/'.runtime/v18_20_original_server_sendback.zip') as z:
        seeds=json.loads(z.read(SRC20+'results/matched_hidden_seeds.json'))
    wanted=('material3_jet2','raw4_jet1','raw6_jet1')
    seeds={s['name']:s for s in seeds if s['name'] in wanted}
    if set(seeds)!=set(wanted):raise RuntimeError('Original matched hidden seeds incomplete')
    rng=np.random.default_rng(CONFIG['seed'])
    rows=[];spectra=[];deterministic=[]
    for label in wanted:
        info=seeds[label];key=info['observer'];H=Hs[key]
        x=np.asarray(info['seed'],float);x/=np.linalg.norm(x);m=H.shape[0]
        jet_mats=[H@mat for mat in Mm]
        for width in CONFIG['halfwidths']:
            offsets=basis_times(width,CONFIG['sample_count'])
            actual=offsets
            times=4.+actual
            O=continuous_observation(spline,P[i],times,H)
            Os=continuous_observation(spline_standard,Pst[i],times,H)
            err_tight_standard=float(np.linalg.norm(O-Os)/max(np.linalg.norm(O),1e-14))
            pinv_direct=direct_estimator(O)
            sv=np.linalg.svd(O,compute_uv=False)
            for degree in CONFIG['poly_degrees']:
                B4,C,W,Oj,Hjet=jet_compression(actual,width,degree,m,jet_mats)
                pinv_jet=jet_estimator(Oj)
                check=test_covariance(C,m,W)
                if check>1e-7:raise RuntimeError('Jet whitening/correlation error')
                # The correct regression jet ideal and the polynomial observation differ
                # because the real trajectory contains 4th+ order terms.
                bias=nominal_jet_bias(O,x,B4,W,m,Oj)
                condj=np.linalg.svd(Oj,compute_uv=False)
                spectra.append(dict(label=label,observer=key,window_halfwidth=width,poly_degree=degree,
                  min_singular_direct=float(sv[-1]),min_singular_jet=float(condj[-1]),
                  condition_direct=float(sv[0]/sv[-1]),condition_jet=float(condj[0]/condj[-1]),
                  covariance_identity_error=check,standard_tight_relative=err_tight_standard,
                  jet_model_bias_whitened=float(np.linalg.norm(bias)),
                  noiseless_direct_state_error=float(np.linalg.norm(pinv_direct@(O@x)-x)),
                  noiseless_jet_state_error=float(np.linalg.norm(pinv_jet@compressed_observation(O@x,m,B4,W)-x))))
                raw_clean=(O@x).reshape(-1,m)
                for sigma in CONFIG['noise_levels']:
                    dsq=[];jsq=[];dproj=[];jproj=[]
                    for k in range(CONFIG['repeats']):
                        noise=rng.normal(0,sigma,size=raw_clean.shape)
                        observed=raw_clean+noise
                        yd=observed.reshape(-1)
                        xd=pinv_direct@yd
                        xj=pinv_jet@compressed_observation(yd,m,B4,W)
                        dsq.append(np.linalg.norm(xd-x))
                        jsq.append(np.linalg.norm(xj-x))
                        # direction-specific amplitude recovery, not full-state recovery
                        dproj.append(float(x@xd))
                        jproj.append(float(x@xj))
                    rows.append(dict(label=label,observer=key,window_halfwidth=width,
                      sample_count=CONFIG['sample_count'],degree=degree,noise_sigma=sigma,trials=CONFIG['repeats'],
                      direct_rms=float(np.sqrt(np.mean(np.square(dsq)))),
                      jet_rms=float(np.sqrt(np.mean(np.square(jsq)))),
                      direct_projection_mean=float(np.mean(dproj)),
                      jet_projection_mean=float(np.mean(jproj)),
                      direct_projection_sd=float(np.std(dproj,ddof=1)),
                      jet_projection_sd=float(np.std(jproj,ddof=1)),
                      direct_snr_seed=float(np.linalg.norm(O@x)/sigma),
                      jet_snr_seed=float(np.linalg.norm(Oj@x)/sigma),
                      poly_jet_bias_norm=float(np.linalg.norm(bias)),
                      jet_over_direct_rms=float(np.sqrt(np.mean(np.square(jsq)))/max(np.sqrt(np.mean(np.square(dsq))),1e-14))))
    write_csv('raw_sample_reconstruction_trials.csv',rows)
    write_csv('observation_singular_spectrum.csv',spectra)
    worst=max(z['standard_tight_relative'] for z in spectra)
    output={'version':'v18.22','status':'FROZEN_ORIGINAL_Q550_RAW_SAMPLE_TANGENT_PREFLIGHT',
      'config':CONFIG,'provenance':prov,'predecessor_nonlinear_server_status':newproof['status'],
      'matched_seeds':wanted,'full_rank_all':all(z['min_singular_direct']>1e-10 and z['min_singular_jet']>1e-10 for z in spectra),
      'whiten_max_defect':max(z['covariance_identity_error'] for z in spectra),
      'max_standard_tight_window_drift':worst,
      'limitations':['Raw samples are synthetically noisy observations of the certified tangent propagator, not measured sensor readings.',
        'Noise acts on differences of paired trajectories normalized by their amplitude; physical noise scale is unknown.',
        'Polynomial reconstruction has finite-window truncation bias and correlated coefficient errors.',
        'Jet map uses local third-order Taylor coefficients even for degree-five polynomial regression.',
        'No intrinsic physical metric or full nonlinear history-dependent constitutive law is certified.',
        'Direct and compressed estimators both use model knowledge and a fixed 10D diagnostic gauge.']}
    (ROOT/'results/raw_sample_summary.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':output['status'],'rows':len(rows),'max_standard_tight_window_drift':worst,'all_full_rank':output['full_rank_all'],
      'representative':[z for z in rows if z['window_halfwidth']==.2 and z['noise_sigma']==1e-4 and z['degree']==5]},indent=2))
    if not output['full_rank_all']:raise RuntimeError('Numerical rank deficiency in raw/jet maps')
    return output
if __name__=='__main__':reconstruct()
