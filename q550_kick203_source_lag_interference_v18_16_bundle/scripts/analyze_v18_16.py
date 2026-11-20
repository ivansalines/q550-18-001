#!/usr/bin/env python3
"""Q550 v18.16 source-lag interference / echo persistence diagnostics.

No new ODE replay; no physical source-history removal. The fixed trajectory is
cut into time-of-birth slabs. This audit asks whether an apparent ~1.5 xi lag
is a narrow temporal feature or a selection effect of strongest individual
negative slab pair at four preselected events.
"""
import argparse,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from checkpoint_io_v18_16 import read
from lag_kernel_math import age_slabs,kernel,pairing_identity,effective_lag_minimum,kernel_from_slabs,project
WIDTHS=(.25,.125,.0625,.03125)
LAGS=np.arange(.25,2.25001,.25)
PROFILE_LAGS=np.arange(.25,2.25001,.125)
# Same horizon for matched support; full-history variant available separately.
HORIZON=3.0
MIN_OVERLAP=.5

def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def metrics():
    return {'euclidean':np.ones(10),
            'moderate_diag':np.linspace(.8,1.2,10),
            'strong_diag':np.geomspace(.5,2.0,10)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'results');args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    data,previous,prov=read()
    events=previous['events_xi']
    records=[];events_out=[];lag_profile=[];checks={};independent=[]
    for tol in ('standard','tight'):
        d=data[tol];x=d['xi_nodes'];A=d['integral_source_frame_11D'][:,:10]
        if x.shape!=(641,) or A.shape!=(641,10) or abs(x[80]-1)>1e-12 or np.linalg.norm(A[0])>1e-10:
            raise RuntimeError('Unexpected original archived source or phase grid')
        interp=PchipInterpolator(x,A,axis=0)
        integ=PchipInterpolator(x,d['source_integrand_11D'][:,:10],axis=0).antiderivative()
        slope=float(d['slope'])
        e=interp(1.)
        if np.linalg.norm(e)<1e-7:raise RuntimeError('Degenerate historical anchor')
        for metric_name,weight in metrics().items():
            for j,z in enumerate(events[tol],1):
                for width in WIDTHS:
                    for window in ('post_anchor_full','fixed_3xi'):
                        horizon=None if window=='post_anchor_full' else HORIZON
                        rr,covered=kernel(interp,e,z,width,LAGS,weight,horizon,MIN_OVERLAP)
                        mloc,mval=effective_lag_minimum(rr)
                        for lag,k,overlap,npairs in rr:
                            records.append({'tolerance':tol,'metric':metric_name,'event':j,'event_xi':z,
                                            'event_type':'radial_max' if j%2 else 'radial_min',
                                            'history':window,'bin_width':width,'lag':lag,
                                            'kernel':None if not np.isfinite(k) else float(k),
                                            'support_xi':overlap,'n_pairs':npairs,
                                            'history_covered_xi':covered})
                        if j in (2,4) and window=='post_anchor_full' and metric_name=='euclidean':
                            q,covered=age_slabs(interp,e,z,width,weight,1.,None)
                            edges=np.arange(len(q)+1)*width
                            q_ind=project(e,slope*(integ(z-edges[:-1])-integ(z-edges[1:])),weight)
                            k=int(round(1.5/width))
                            pred,_=kernel_from_slabs(q,k,weight)
                            check,_=kernel_from_slabs(q_ind,k,weight)
                            independent.append({'tolerance':tol,'event':j,'bin_width':width,
                                                'kernel_from_A':float(pred),'kernel_from_integrated_g':float(check),
                                                'abs_disagreement':float(abs(pred-check)),
                                                'relative_slabs_L2':float(np.linalg.norm(q-q_ind)/np.linalg.norm(q))})
                        events_out.append({'tolerance':tol,'metric':metric_name,'event':j,'event_xi':z,
                                           'history':window,'bin_width':width,
                                           'most_negative_lag_1to2':mloc,'most_negative_value_1to2':mval,
                                           'kernel_at_1_5':next((float(r[1]) for r in rr if abs(r[0]-1.5)<1e-9 and np.isfinite(r[1])),None),
                                           'history_covered_xi':covered})
            # Continuous profile evaluated at grid-aligned observed xi; no event-only cherry picking.
            if metric_name=='euclidean':
                for z in x[::4]:  # every 0.05 xi, along whole checked interval
                    if z<1.05:continue
                    for window in ('post_anchor_full','fixed_3xi'):
                        h=None if window=='post_anchor_full' else HORIZON
                        rr,covered=kernel(interp,e,float(z),.0625,PROFILE_LAGS,weight,h,MIN_OVERLAP)
                        mloc,mval=effective_lag_minimum(rr)
                        for L,k,support,npairs in rr:
                            lag_profile.append({'tolerance':tol,'history':window,'xi':float(z),'lag':L,
                                                'kernel':None if not np.isfinite(k) else float(k),
                                                'support_xi':support,'n_pairs':npairs,
                                                'argmin_lag_1to2':mloc,'argmin_value_1to2':mval})
        # precision closure and reference-anchor check independent of kernel
        w=np.ones(10);q,cover=age_slabs(interp,e,6.,.0625,w,anchor=1.)
        residual=np.linalg.norm(q.sum(axis=0)-(interp(6.)-interp(1.) - (np.dot(interp(6.)-interp(1.),e)/np.dot(e,e))*e))
        pbefore=interp(1.)-interp(0.)
        pzero=pbefore - (np.dot(pbefore,e)/np.dot(e,e))*e
        checks[tol]={'slab_closure_residual_xi6':float(residual),'pre_anchor_projected_zero':float(np.linalg.norm(pzero)),
                     'pair_energy_identity_error_xi6':pairing_identity(q,w),
                     'anchor_norm':float(np.linalg.norm(e))}
        if residual>1e-8 or np.linalg.norm(pzero)>1e-8:raise RuntimeError('Post-anchor or slab closure failed')
    # Convert into same metric/lag support region for matched-support persistence.
    stats=[]
    for tol in ('standard','tight'):
        for window in ('post_anchor_full','fixed_3xi'):
            for L in (1.25,1.5,1.75):
                one=[r for r in lag_profile if r['tolerance']==tol and r['history']==window and r['lag']==L and r['kernel'] is not None and r['xi']>=4.0]
                score=np.array([r['kernel'] for r in one]);valid=[r['xi'] for r in one]
                minus=score < 0
                # Run length of consecutive sampled negative-kernel windows, step 0.05 xi
                best_run=cur=0
                for b in minus:
                    cur=cur+1 if b else 0;best_run=max(best_run,cur)
                near2=[r['kernel'] for r in one if abs(r['xi']-events[tol][1])<=.25]
                near4=[r['kernel'] for r in one if abs(r['xi']-events[tol][3])<=.25]
                away=[r['kernel'] for r in one if abs(r['xi']-events[tol][1])>.5 and abs(r['xi']-events[tol][3])>.5]
                stats.append({'tolerance':tol,'history':window,'lag':L,'n_xi':len(score),
                              'negative_fraction':float(minus.mean()) if score.size else None,
                              'longest_negative_run_xi':float(best_run*.05),
                              'median_kernel':float(np.median(score)) if score.size else None,
                              'median_near_min_2':float(np.median(near2)) if near2 else None,
                              'median_near_min_4':float(np.median(near4)) if near4 else None,
                              'median_away_from_minima':float(np.median(away)) if away else None})
    # Replication test standard/tight on events in the same metric/partition.
    std={(r['metric'],r['event'],r['history'],r['bin_width'],r['lag']):r['kernel'] for r in records if r['tolerance']=='standard' and r['kernel'] is not None}
    diffs=[abs(r['kernel']-std[(r['metric'],r['event'],r['history'],r['bin_width'],r['lag'])]) for r in records if r['tolerance']=='tight' and r['kernel'] is not None and (r['metric'],r['event'],r['history'],r['bin_width'],r['lag']) in std]
    drift=float(np.max(abs(np.array(events['standard'])-np.array(events['tight']))))
    if drift>1e-5 or max(diffs)>1e-5:raise RuntimeError('Standard/tight stability regression failure')
    write_csv(args.output/'event_lag_kernel.csv',records)
    write_csv(args.output/'event_lag_minima.csv',events_out)
    write_csv(args.output/'kernel_xi_lag_profile.csv',lag_profile)
    write_csv(args.output/'lag_persistence_stats.csv',stats)
    write_csv(args.output/'independent_source_check.csv',independent)
    chosen=[r for r in events_out if r['tolerance']=='standard' and r['metric']=='euclidean' and r['bin_width']==.0625 and r['history']=='post_anchor_full']
    selected=[r for r in records if r['tolerance']=='standard' and r['metric']=='euclidean' and r['bin_width']==.0625 and r['history']=='post_anchor_full' and r['event'] in (2,4)]
    notes=[]
    for eidx in (2,4):
        rr=[r for r in selected if r['event']==eidx]
        half=next(r for r in rr if abs(r['lag']-1.5)<1e-8)
        b=next(r for r in chosen if r['event']==eidx)
        notes.append({'event':eidx,'xi':float(b['event_xi']),'kernel_lag_1_5':half['kernel'],
                      'strongest_lag_in_1to2':b['most_negative_lag_1to2'],
                      'strongest_lag_kernel':b['most_negative_value_1to2']})
    summary={'version':'Q550-v18.16-source-lag-interference-kernel-echo-persistence',
             'status':'FROZEN_SOURCE_LAG_OBSERVER_AUDIT',
             'model':'Q550 kick203 original archived v18.7 source A in phase-fixed 10D co-moving coordinates',
             'projection':'metric-orthogonal to e=A(1); birth history constrained to u>=1',
             'geometry':'declared diagonal positive-definite weighted metrics, not established physical metric',
             'lag_grid_profile_xi':PROFILE_LAGS.tolist(), 'event_kernel_lags':LAGS.tolist(),
             'age_resolutions_xi':WIDTHS,'fixed_history_horizon_xi':HORIZON,
             'minimum_overlap_xi':MIN_OVERLAP,'events':events,
             'event_euclidean_full_postanchor_bin_00625':chosen,
             'event_2_and_4':notes,
             'persistence_stats':stats,'closure':checks,
             'max_independent_A_vs_g_kernel_disagreement':float(max(r['abs_disagreement'] for r in independent)),
             'standard_tight_max_event_xi_drift':drift,
             'standard_tight_max_kernel_abs_difference':float(max(diffs)),
             'provenance':prov,
             'limitations':[
               'The fixed original trajectory has not been re-integrated or physically history-ablated.',
               'Kernel, sign, strongest lag, angular interpretation and negative-run fraction depend on metric, anchor, bin resolution, lag support, and finite history window.',
               'The paired support shrinks with lag: compare only supported lag/xi combinations; no stationary process or universal delay is presumed.',
               'Persistent correlation of a smooth fixed source is not independent evidence of constitutive material memory or mechanical torsion.',
               'Four earlier event positions are inherited from a metric-dependent radial observable.',
               'This source-pair correlation differs from the maximal opposing single slab-pair diagnostic of v18.15. They need not identify the same lag.',
               'Fraction of samples with negative kernel does not measure fraction of energy or physical material memory.'
             ]}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    lines=['# Q550 v18.16 — Source-Lag Interference Kernel & Echo Persistence Audit','',
           '## Scope','', 'This audit evaluates all xi checkpoints (sampled 0.05 xi) rather than only four selected extrema. It compares the full post-anchor history u>=1 against a fixed recent 3 xi window, then varies slab width and diagonal metric. Frozen original Q550 kick203 v18.7 source inputs remain unchanged.','',
           '## Event summary (Euclidean metric, source from xi>=1, slab width 0.0625)','',
           '| Event | xi | K(lag=1.5) | Most negative lag in [1,2] | K at that lag |','|---|---:|---:|---:|---:|']
    for r in chosen:
        fmt=lambda v: 'n/a' if v is None else f'{v:.4f}'
        lines.append(f"| {r['event']} | {r['event_xi']:.6f} | {fmt(r['kernel_at_1_5'])} | {fmt(r['most_negative_lag_1to2'])} | {fmt(r['most_negative_value_1to2'])} |")
    lines.extend(['','## Persistence (xi>=4, 0.05 xi sample spacing)','',
                  '| Window | Lag xi | Negative fraction | Max negative run xi | Median K |','|---|---:|---:|---:|---:|'])
    for r in stats:
        if r['tolerance']=='standard':lines.append(f"| {r['history']} | {r['lag']:.2f} | {r['negative_fraction']:.3f} | {r['longest_negative_run_xi']:.2f} | {r['median_kernel']:.4f} |")
    lines.extend(['','## Checks','',
        f"Largest standard/tight difference of kernel = {max(diffs):.3e}.",
        f"Maximum standard/tight event drift = {drift:.3e} xi.",
        f"Largest independent source-integrand A-versus-g kernel disagreement = {max(r['abs_disagreement'] for r in independent):.3e}.",
        'Closure checks are in summary.json. Resolution and metric sensitivity are explicit in event_lag_kernel.csv; changes in strongest lag are not a failure but evidence against metric-independent localization.',
        '', '## Interpretation','',
        'The signed K is a normalized sum of dot products of age-offset source slabs in a declared 10D metric. A negative value describes opposing vectors across a time separation, not anti-causality, work, physical reversal or hysteresis.',
        'The previous ~1.5 xi value was selected from the single most opposing slab pair at two preselected extrema; this kernel averages matched pairs over their available overlap and can refute that value as a broadly persistent lag.',
        'For genuine constitutive memory, compare physically realizable trajectories with matched state and phase but differing histories; this audit does not perform that test.',
        '', '## Genealogy','',
        f"Original authentic v18.15 server sendback SHA256: `{prov['v18_15_original_sha256']}`.",
        'Nested original v18.14 -> v18.7/6 data and their source scripts are inside inputs/v18_15_original_sendback.zip. The v18.15 anchor/lag addendum script is retained separately under reference/.',
        'All code and report text are in English. No baseline 400 orbit replay.',''])
    (args.output/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'status':summary['status'],'event_2_4':notes,'max_tight_drift':max(diffs)},indent=2))
if __name__=='__main__':main()
