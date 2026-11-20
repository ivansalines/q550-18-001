#!/usr/bin/env python3
"""v18.18 — Exact source recurrence closure & observer-model residual audit.

The full source-integrated v18.17 kernel is algebraically a weighted average of
instantaneous pair correlations. In a single trajectory, a zero 'full source'
residual is expected by construction, not physical proof of absent memory.
Nonzero residuals from simpler predictors mean model incompleteness only.
"""
import csv, json, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from checkpoint_io_v18_18 import read
from recurrence_math import paired_decomposition,instantaneous,train_endpoint_calibrator,ridge_lag
TIMES=np.round(np.arange(4.,8.0001,.05),10)
LAGS=np.round(np.arange(1.,2.0001,.025),10)
METRICS={'euclidean':np.ones(10), 'strong_diag':np.geomspace(.5,2.,10)}

def row_stats(arr):
    a=np.asarray(arr,float)
    return {'max_abs':float(np.max(np.abs(a))), 'median_abs':float(np.median(np.abs(a))),
            'rmse':float(np.sqrt(np.mean(a*a))), 'p90_abs':float(np.quantile(np.abs(a),.9))}

def save_csv(p,rows):
    if not rows:return
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in rows for k in row.keys())));w.writeheader();w.writerows(rows)

def archived_surface():
    """Open authenticated original v18.17 comparison table without repacking it."""
    zpath=ROOT/'inputs/v18_17_original_sendback.zip'
    member='q550_kick203_moving_antagonism_ridge_v18_17_bundle/results/moving_ridge_surface.csv'
    result={}
    with zipfile.ZipFile(zpath) as z:
        lines=z.read(member).decode('utf-8').splitlines()
    for row in csv.DictReader(lines):
        if row['method']!='matched_1xi':continue
        key=(row['tolerance'],row['metric'])
        result.setdefault(key,{})[(float(row['xi']),float(row['lag']))]=float(row['kernel'])
    return result

def main():
    RESULT=ROOT/'results';RESULT.mkdir(exist_ok=True)
    DATA,previous,provenance=read()
    previous_surface=archived_surface()
    surfaces=[];timeseries=[];crosschecks=[];events=[]
    tol_cache={}
    with zipfile.ZipFile(ROOT/'inputs/v18_17_original_sendback.zip') as z:
        s17=json.loads(z.read('q550_kick203_moving_antagonism_ridge_v18_17_bundle/results/summary.json'))
    past_events={int(x['event']):float(x['event_xi']) for x in s17['event_comparison'] if x['tolerance']=='standard'}
    for tolerance in ('standard','tight'):
        row=DATA[tolerance]
        cx=PchipInterpolator(row['xi_nodes'],row['integral_source_frame_11D'][:,:10],axis=0)
        source=cx.derivative()
        original=PchipInterpolator(row['xi_nodes'],float(row['slope'])*row['source_integrand_11D'][:,:10],axis=0)
        e=cx(1.)
        for metric_name,w in METRICS.items():
            D=paired_decomposition(source,TIMES,LAGS,w,e,80)
            E=instantaneous(source,TIMES,LAGS,w,e)
            Cal,coef=train_endpoint_calibrator(E,D['exact'],TIMES,6.)
            assert np.max(np.abs(D['exact']-D['closure']))<1e-12
            expected=np.array([[previous_surface[(tolerance,metric_name)][(float(t),float(l))]
                        for l in LAGS] for t in TIMES])
            archive_diff=row_stats(D['exact']-expected)
            if archive_diff['max_abs']>1e-9: raise RuntimeError(f'v18.17 reproduction mismatch: {tolerance} {metric_name}: {archive_diff}')
            indep=paired_decomposition(original,TIMES,LAGS,w,e,80)['exact']
            finer=paired_decomposition(source,TIMES,LAGS,w,e,160)['exact'] if metric_name=='euclidean' else None
            coarser=paired_decomposition(source,TIMES,LAGS,w,e,40)['exact'] if metric_name=='euclidean' else None
            diag={
                'tolerance':tolerance,'metric':metric_name,
                'v18_17_reproduction':archive_diff,
                'exact_source_closure':row_stats(D['exact']-D['closure']),
                'independent_original_integrand':row_stats(D['exact']-indep),
                'unit_direction_only_residual':row_stats(D['exact']-D['unit_direction']),
                'amplitude_gated_unit_residual':row_stats(D['covariance']),
                'instantaneous_endpoint_residual':row_stats(D['exact']-E),
                'scrambled_age_pair_residual':row_stats(D['exact']-D['scrambled']),
                'heldout_xi_gt_6':{
                    'endpoint_raw':row_stats((D['exact']-E)[TIMES>6]),
                    'endpoint_calibrated_train4to6':row_stats((D['exact']-Cal)[TIMES>6]),
                    'unit_only':row_stats((D['exact']-D['unit_direction'])[TIMES>6]),
                    'amplitude_gated_unit':row_stats(D['covariance'][TIMES>6]),
                },
                'endpoint_calibration':{'intercept':float(coef[0]),'slope':float(coef[1])},
                'width_convergence_40_vs_80':None if coarser is None else row_stats(coarser-D['exact']),
                'width_convergence_80_vs_160':None if finer is None else row_stats(finer-D['exact'])
            }
            crosschecks.append(diag)
            # grid-point closure and drift, distinguish heldout actual from fitted.
            for i,z in enumerate(TIMES):
                rids={k:float(ridge_lag(LAGS,arr[i])) for k,arr in
                     {'exact':D['exact'],'unit_only':D['unit_direction'],'amp_gated_unit':D['gate_times_unit'],
                      'instantaneous':E,'endpoint_calibrated':Cal}.items()}
                timeseries.append(dict(tolerance=tolerance,metric=metric_name,xi=float(z),
                    split='train' if z<=6 else 'holdout',
                    **{f'lag_{k}':v for k,v in rids.items()},
                    min_kernel_exact=float(np.min(D['exact'][i])),
                    max_abs_closure=float(np.max(np.abs(D['exact'][i]-D['closure'][i]))),
                    rms_covariance_residual=float(np.sqrt(np.mean(D['covariance'][i]**2))),
                    rms_endpoint_residual=float(np.sqrt(np.mean((D['exact'][i]-E[i])**2))),
                    rms_calibrated_residual=float(np.sqrt(np.mean((D['exact'][i]-Cal[i])**2)))))
                for j,lag in enumerate(LAGS):
                    surfaces.append(dict(tolerance=tolerance,metric=metric_name,xi=float(z),lag=float(lag),
                        exact=float(D['exact'][i,j]),closure=float(D['closure'][i,j]),
                        unit_only=float(D['unit_direction'][i,j]),amp_gate=float(D['gate'][i,j]),
                        amplitude_gated_unit=float(D['gate_times_unit'][i,j]),
                        covariance_residual=float(D['covariance'][i,j]),
                        instantaneous=float(E[i,j]),endpoint_calibrated=float(Cal[i,j]),
                        scrambled=float(D['scrambled'][i,j])))
            if metric_name=='euclidean' and tolerance=='standard':
                for event_no,z in past_events.items():
                    if z < 4:
                        events.append(dict(event=event_no,event_xi=z,coverage='xi_below_profile_range',sample_xi=None))
                        continue
                    i=int(np.argmin(abs(TIMES-z)))
                    events.append(dict(event=event_no,event_xi=z,coverage='nearest_profile_xi',sample_xi=float(TIMES[i]),
                        ridge_exact=float(ridge_lag(LAGS,D['exact'][i])),
                        ridge_instantaneous=float(ridge_lag(LAGS,E[i])),
                        ridge_unit_only=float(ridge_lag(LAGS,D['unit_direction'][i])),
                        rms_amplitude_covariance=float(np.sqrt(np.mean(D['covariance'][i]**2))),
                        rms_endpoint_residual=float(np.sqrt(np.mean((D['exact'][i]-E[i])**2)))))
            tol_cache[(tolerance,metric_name)]=D['exact']
    std_tight={name:row_stats(tol_cache[('standard',name)]-tol_cache[('tight',name)]) for name in METRICS}
    result={'version':'Q550-v18.18-source-recurrence-cancellation-echo-residual',
       'status':'FROZEN_SOURCE_EXACT_CLOSURE_AND_SURROGATE_AUDIT',
       'source':'Authenticated v18.17 server sendback with nested original v18.7 source',
       'profile':{'xi':[4,8,.05],'lag':[1,2,.025],'age_support':[0,1],
                  'quadrature_nodes':80,'training_xi':[4,6],'holdout_xi':'6.05..8'},
       'closure_identity':'K=SUM_age ||q_z|| ||q_old|| cos(q_z,q_old) / sqrt(SUM||q_z||² SUM||q_old||²).',
       'covariance_definition':'K - gate * MEAN_age cos(pair); diagnostic of amplitude/direction variation',
       'source_recurrence_vs_transport':'All arrays derived from one frozen transported source. A source-only identity cannot test an independent material memory.',
       'metrics':crosschecks,'standard_tight':std_tight,'events':events,'provenance':provenance,
       'interpretation_limits':[
         'Full source closure is an algebraic identity: vanishing residual is expected, not proof of no physical memory.',
         'Residual to amplitude-blind or endpoint-only model is model misspecification, not newly discovered constitutive memory.',
         'Train/test split of one autocorrelated trajectory does not provide independent physical-history replication.',
         'Pairwise source correlation and ridge timing depend on chosen observer metric, anchor, projection and window.',
         'A fixed phase 10D linear perturbation subspace is not automatically a material constraint manifold.',
         'Scrambled age pairing is a diagnostic null, not a physically realizable Q550 trajectory.',
         'The lag is ξ separation, not a calibrated physical memory time.']}
    save_csv(RESULT/'source_kernel_decomposition.csv',surfaces)
    save_csv(RESULT/'time_residual_and_ridge.csv',timeseries)
    save_csv(RESULT/'event_probe.csv',events)
    (RESULT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# Q550 v18.18 — Source Recurrence Cancellation & Transported Echo Residual Audit','','## Scientific conclusion','','The exact weighted quadrature of instantaneous projected source pairs reproduces the v18.17 integrated kernel to numerical roundoff. **This is an algebraic identity**: it proves the observable kernel is fully determined by the supplied source on the same frozen trajectory. It neither certifies nor excludes physical hereditary dynamics inside that source.','','Nonzero discrepancies of amplitude-blind, endpoint-only and age-scrambled constructions show limitations of those observer baselines, not constitutive memory. The source-only identity is a high-strength falsification control against overinterpreting an apparent residual.','','## Diagnostic table','','| Gauge | Closure max abs | Independent-g max abs | Amplitude-gated unit residual RMSE | Endpoint holdout RMSE | Endpoint fitted holdout RMSE |','|---|---:|---:|---:|---:|---:|']
    for c in crosschecks:
        lines.append(f"| {c['tolerance']}/{c['metric']} | {c['exact_source_closure']['max_abs']:.3e} | {c['independent_original_integrand']['max_abs']:.3e} | {c['amplitude_gated_unit_residual']['rmse']:.4f} | {c['heldout_xi_gt_6']['endpoint_raw']['rmse']:.4f} | {c['heldout_xi_gt_6']['endpoint_calibrated_train4to6']['rmse']:.4f} |")
    lines+=['','## Definitions','',
      '- Pair direction `c(a)=<q(z-a),q(z-a-L)>_G/(||q(z-a)||_G ||q(z-a-L)||_G)`.',
      '- Exact `K=∫r(a)s(a)c(a) da/sqrt(∫r(a)^2 da ∫s(a)^2 da)`.',
      '- Magnitude gate `H=∫rs da/sqrt(∫r² da ∫s² da)`; directional-only predictor `mean(c)`; amplitude-gated predictor `H*mean(c)`.',
      '- Covariance residual `K-H*mean(c)=H*(weighted_mean(c)-mean(c))` is completely determined by age-varying instantaneous source norms and orientations.',
      '- Endpoint-only prediction is `cos(q(z),q(z-L))`; calibrated endpoint has 2 parameters trained on xi=4..6, tested on xi>6.',
      '- Age scrambled pairs roll the older source across the age grid; not a physical counterfactual.',
      '', '## Checks','','- Authenticated CRC/SHA256 for the original v18.17 server sendback and the nested genealogy through v18.7.',
      '- Recomputed v18.17 kernel from its independent checkpoint and compared against original archived surface values.',
      '- Compared derivative of cumulative A to separately saved original source integrand.',
      '- Standard/tight; alternate diagonal metric; quadrature n=40,80,160; delayed-source age permutation.',
      '', '## Next falsifiable step','','Reconstruct independently perturbed ODE trajectories with matched reduced observables and phase, and explicitly separate linear predicted differences from nonlinear response. Do not claim memory merely from same-trajectory source correlation.','']
    (RESULT/'REPORT.md').write_text('\n'.join(lines))
    show=crosschecks[0]
    print(json.dumps({'status':result['status'],'closure_max':show['exact_source_closure']['max_abs'],
      'source_integrand_max_drift':show['independent_original_integrand']['max_abs'],
      'endpoint_holdout_rmse':show['heldout_xi_gt_6']['endpoint_raw']['rmse'],
      'endpoint_fitted_holdout_rmse':show['heldout_xi_gt_6']['endpoint_calibrated_train4to6']['rmse'],
      'amplitude_gated_residual_rmse':show['amplitude_gated_unit_residual']['rmse'],
      'standard_tight_max':std_tight['euclidean']['max_abs']},indent=2))
if __name__=='__main__':main()
