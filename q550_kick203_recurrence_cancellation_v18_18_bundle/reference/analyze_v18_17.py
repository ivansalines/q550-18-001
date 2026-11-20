#!/usr/bin/env python3
"""Q550 v18.17 — continuous moving antagonism ridge and source recurrence.

Uses exact v18.16 ORIGINAL server sendback and underlying frozen source checkpoint.
Does not integrate Q550 or assert physical memory, physical lag, or material torsion.
"""
import csv,json,argparse,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from checkpoint_io_v18_17 import read
from ridge_math import fixed_support,instantaneous,variable_support,ridge_row,event_nearest

TIMES=np.round(np.arange(4.,8.0001,.05),10)
LAGS=np.round(np.arange(1.,2.00001,.025),10)
METRICS={'euclidean':np.ones(10),
         'moderate_diag':np.linspace(.8,1.2,10),
         'strong_diag':np.geomspace(.5,2.,10)}

def dump_csv(file,rows):
    if not rows:return
    with file.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def summary_stats(rows,selector):
    hits=[r for r in rows if selector(r) and r['lag_quad'] is not None]
    if not hits:return {}
    lag=np.array([r['lag_quad'] for r in hits]);h=np.array([r['kernel_min'] for r in hits]);b=np.array([r['span_002'] for r in hits]);
    return {'n_xi':len(hits),'median_argmin_lag':float(np.median(lag)),
            'lag_range':[float(lag.min()),float(lag.max())],
            'lag_std':float(lag.std()),'median_kernel':float(np.median(h)),
            'negative_fraction':float(np.mean(h<0)),
            'broad_ridge_span_002_median_xi':float(np.median(b)),
            'edge_fraction':float(np.mean([r['edge'] for r in hits]))}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    DATA,prev,provenance=read()
    surfaces=[];ridges=[];check_rows=[];source_checks=[];events=[]
    numerical={}
    for tolerance in ['standard','tight']:
        raw=DATA[tolerance];xx=raw['xi_nodes']
        cumulative=PchipInterpolator(xx,raw['integral_source_frame_11D'][:,:10],axis=0)
        gA=cumulative.derivative()
        gG=PchipInterpolator(xx,float(raw['slope'])*raw['source_integrand_11D'][:,:10],axis=0)
        grid=np.linspace(1.,8.,561)
        mismatch=np.linalg.norm(gA(grid)-gG(grid),axis=1)
        norm=np.linalg.norm(gG(grid),axis=1)
        source_checks.append({'tolerance':tolerance,'max_absolute_source_mismatch':float(np.max(mismatch)),
                              'median_relative_source_mismatch':float(np.median(mismatch/np.maximum(norm,1e-12))),
                              'max_relative_source_mismatch':float(np.max(mismatch/np.maximum(norm,1e-12)))})
        for metric_name,weights in METRICS.items():
            e=cumulative(1.)
            K={}
            for samples in (20,40,80):
                K[samples]=fixed_support(gA,TIMES,LAGS,samples,weights,e)
            max_diff_20_80=float(np.max(np.abs(K[20]-K[80])))
            max_diff_40_80=float(np.max(np.abs(K[40]-K[80])))
            source_independent=fixed_support(gG,TIMES,LAGS,80,weights,e)
            max_source_diff=float(np.max(np.abs(source_independent-K[80])))
            g_now=instantaneous(gA,TIMES,LAGS,weights,e)
            if metric_name in ('euclidean','strong_diag'):
                k_variable=variable_support(gA,TIMES,LAGS,weights,e,step=.025,horizon=3.)
            else:k_variable=None
            near=fixed_support(gA,TIMES,LAGS,80,weights,e,age_lo=0.,age_hi=.5)
            far=fixed_support(gA,TIMES,LAGS,80,weights,e,age_lo=.5,age_hi=1.)
            observations={'matched_1xi':K[80], 'recent_half_xi':near,'older_half_xi':far,
                          'instantaneous_pair':g_now}
            if k_variable is not None:observations['variable_support_3xi']=k_variable
            for kind,grid_values in observations.items():
                for i,t in enumerate(TIMES):
                    ridge=ridge_row(LAGS,grid_values[i]);ridges.append(dict(tolerance=tolerance,metric=metric_name,method=kind,xi=float(t),**ridge))
                    for j,L in enumerate(LAGS):
                        surfaces.append({'tolerance':tolerance,'metric':metric_name,'method':kind,
                                         'xi':float(t),'lag':float(L),'kernel':float(grid_values[i,j])})
            ridge_coarse=np.array([ridge_row(LAGS,r)['lag_quad'] for r in K[20]])
            ridge_fine=np.array([ridge_row(LAGS,r)['lag_quad'] for r in K[80]])
            check_rows.append({'tolerance':tolerance,'metric':metric_name,
                               'max_k_quadrature_20_80':max_diff_20_80,
                               'max_k_quadrature_40_80':max_diff_40_80,
                               'max_ridge_lag_20_80':float(np.max(np.abs(ridge_coarse-ridge_fine))),
                               'max_independent_A_derivative_vs_original_g_K':max_source_diff,
                               'paired_negative_fraction_1xi':float(np.mean(K[80]<0)),
                               'instant_negative_fraction':float(np.mean(g_now<0)),
                               'flattened_instant_vs_matched_pearson':float(np.corrcoef(g_now.ravel(),K[80].ravel())[0,1])})
            if metric_name=='euclidean':
                for j,z in enumerate(prev['events'][tolerance],1):
                    i=event_nearest(TIMES,z)
                    ev={'tolerance':tolerance,'event':j,'event_xi':float(z),
                        'profile_coverage':'insufficient_post_anchor_history' if z<3. else
                        ('boundary_proxy_0.05xi' if z<TIMES[0] else 'nearest_grid_sample'),
                        'nearest_profile_xi':None if z<3. else float(TIMES[i]),
                        'profile_offset_xi':None if z<3. else float(TIMES[i]-z)}
                    for label,array in observations.items():
                        r=ridge_row(LAGS,array[i]) if z>=3. else dict(lag_quad=None,kernel_min=None)
                        ev[f'{label}_ridge_lag']=r['lag_quad'];ev[f'{label}_min_K']=r['kernel_min']
                    ev['K_at_lag_1_5_fixed']=float(K[80][i,np.argmin(np.abs(LAGS-1.5))]) if z>=3. else None
                    events.append(ev)
            numerical[(tolerance,metric_name)]=K[80]
    # Matched source comparisons (similar gradients need not have matching lag argmin).
    max_tol_drift=max(float(np.max(np.abs(numerical['standard',m]-numerical['tight',m]))) for m in METRICS)
    paired=[r for r in ridges if r['tolerance']=='standard' and r['method']=='matched_1xi' and r['metric']=='euclidean']
    lag_euclid=np.array([r['lag_quad'] for r in paired])
    lag_strong=np.array([r['lag_quad'] for r in ridges if r['tolerance']=='standard' and r['method']=='matched_1xi' and r['metric']=='strong_diag'])
    ridge_variation=float(np.max(np.abs(lag_euclid-lag_strong)))
    # Anchor sensitivity (observation only; source remains physically post-anchor u>=1).
    raw=DATA['standard'];cum=PchipInterpolator(raw['xi_nodes'],raw['integral_source_frame_11D'][:,:10],axis=0)
    e_default=cum(1.)
    anchor=[]
    for z_anchor in (.75,1.,1.25):
        e=cum(z_anchor)
        kval=fixed_support(cum.derivative(),TIMES,LAGS,80,np.ones(10),e)
        anchor.append({'reference_xi':z_anchor,
                       'max_abs_K_change_vs_default':float(np.max(np.abs(kval-numerical['standard','euclidean']))),
                       'max_argmin_lag_change_vs_default':float(np.max(np.abs(
                           np.array([ridge_row(LAGS,k)['lag_quad'] for k in kval])-lag_euclid)))})
    dump_csv(args.output/'moving_ridge_surface.csv',surfaces)
    dump_csv(args.output/'ridge_tracking.csv',ridges)
    dump_csv(args.output/'quadrature_and_independent_source_checks.csv',check_rows)
    dump_csv(args.output/'source_integrand_quality.csv',source_checks)
    dump_csv(args.output/'four_event_ridge_comparison.csv',events)
    dump_csv(args.output/'anchor_sensitivity.csv',anchor)
    metrics_comparison={m:summary_stats(ridges,lambda r:r['tolerance']=='standard' and r['metric']==m and r['method']=='matched_1xi') for m in METRICS}
    method_comparison={m:summary_stats(ridges,lambda r:r['tolerance']=='standard' and r['metric']=='euclidean' and r['method']==m) for m in ('matched_1xi','recent_half_xi','older_half_xi','instantaneous_pair','variable_support_3xi')}
    # Invariance of source subspace dimension can be tested, but argmin of metric kernel is not invariant.
    finding={
        'version':'Q550-v18.17-moving-source-antagonism-ridge',
        'status':'FROZEN_SOURCE_RIDGE_OBSERVER_AUDIT',
        'source':'Authenticated original v18.16 sendback and nested v18.7 11D source; phase-fixed 10D projection',
        'xi_profile':{'from':4.,'to':8.,'step':.05,'n':len(TIMES)},
        'lag_scan':{'from':1.,'to':2.,'step':.025,'n':len(LAGS)},
        'reference':'e=A(1), metric orthogonal projection; source from original archived cumulative A',
        'fixed_paired_age_support_xi':[0.,1.],
        'variable_support_horizon_xi':3.,
        'quadrature_nodes_for_convergence':[20,40,80],
        'ridge_metrics':metrics_comparison,'ridge_methods':method_comparison,
        'max_standard_tight_K_drift':max_tol_drift,
        'max_fixed_ridge_metric_location_change_xi':ridge_variation,
        'source_quality':source_checks,'fixed_quadrature_and_independent_source_checks':check_rows,
        'anchor_sensitivity':anchor,
        'event_comparison':events,'provenance':provenance,
        'interpretation_limits':[
            'Correlation ridge is NOT a mechanical torsion, physical lag, phase locking or constitutive memory estimate.',
            'Finite lag search [1,2] censors any true optimum outside that window; edge minima are flagged.',
            'Quadratic refined argmin is a numerical locator, not a confidence interval.',
            'A broad minimum or strongly metric-dependent lag should not be called a unique ridge.',
            'The fixed 1xi paired support compares like-with-like across lags; the variable-support comparator does not.',
            'The source is from a single computed trajectory; source recurrences do not establish propagation, heredity, or independent memory.',
            'All metrics are declared test metrics, not intrinsic material metrics; xi is not calibrated physical time.',
            'Source derivative of PCHIP is checked against archived integrand; differences are interpolation error.',
            'Instantaneous endpoint similarity does not isolate which source contributions survive dynamical transport.'
        ]}
    (args.output/'summary.json').write_text(json.dumps(finding,indent=2,allow_nan=False)+'\n')
    lines=['# Q550 — v18.17 Moving Ridge, Instantaneous-Source Recurrence & Window-Support Audit','',
           'Original archived Q550 kick203 state, no ODE rerun. All interpretations conditional on declared 10D metric/projection.','',
           '## Fixed-support ridge (time >=4, lag 1..2, paired ages 0..1)','',
           '| Metric | Median best lag | Lag range | Median K at minimum | Broad minimum span (K <= min+0.02) | Edge hits |',
           '|---|---:|---|---:|---:|---:|']
    for key,v in metrics_comparison.items():
        lines.append(f"| {key} | {v['median_argmin_lag']:.4f} | {v['lag_range'][0]:.4f}–{v['lag_range'][1]:.4f} | {v['median_kernel']:.4f} | {v['broad_ridge_span_002_median_xi']:.3f} | {v['edge_fraction']:.1%} |")
    lines+=['','## Sensitivity of the minimum to how recurrences are observed','','| Construction | Median lag | Range | Median K |','|---|---:|---|---:|']
    for k,v in method_comparison.items():lines.append(f"| {k} | {v['median_argmin_lag']:.4f} | {v['lag_range'][0]:.3f}–{v['lag_range'][1]:.3f} | {v['median_kernel']:.4f} |")
    lines+=['','## Independent checks','',
            f'- Maximum standard vs tight fixed-support K drift: `{max_tol_drift:.3e}`.',
            f'- Maximum Euclidean vs strong-diagonal ridge location change: `{ridge_variation:.5f} xi`.',
            f"- Maximum derivative-of-A vs original source-g fixed-support K disagreement: `{max(x['max_independent_A_derivative_vs_original_g_K'] for x in check_rows):.3e}`.",
            f"- Maximum 20-to-80 node K quadrature difference: `{max(x['max_k_quadrature_20_80'] for x in check_rows):.3e}`.",
            f"- Maximum 40-to-80 node K quadrature difference: `{max(x['max_k_quadrature_40_80'] for x in check_rows):.3e}`.",
            '', '## Interpretation','',
            'A negative paired-age kernel may be a recurring pattern of a smooth forcing vector, not a delayed response of material degrees of freedom. The constant paired-age support removes one lag-dependent overlap bias from v18.16; the variable-support and instantaneous comparisons quantify other observer dependencies. No causal propagation law is inferred.',
            '', '## Audit trail','','Original original v18.16 server sendback is retained verbatim under inputs/. Every earlier source script, raw output and checkpoint is nested inside. CRC and SHA256 are checked at runtime.',
            '', '## Next falsifiable question','','Repeat the ridge analysis on independently restarted physical trajectories with matched instantaneous phase and state but controlled, different valid histories; only then attempt a constitutive memory claim.','']
    (args.output/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'status':finding['status'], 'ridge_metrics':metrics_comparison,
                      'max_standard_tight_K_drift':max_tol_drift,'max_source_K_difference':max(x['max_independent_A_derivative_vs_original_g_K'] for x in check_rows)},indent=2))
if __name__=='__main__':main()
