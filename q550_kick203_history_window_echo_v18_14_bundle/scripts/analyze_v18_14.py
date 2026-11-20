#!/usr/bin/env python3
"""Q550 v18.14 fixed-width source-history ablation, event echoes, and rank stress tests.

Diagnostic only: omitting old source from A is *not* re-solving a modified ODE,
nor an experimental preparation of a second state-matched material history.
"""
from __future__ import annotations
import argparse,csv,json,sys,hashlib
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from checkpoint_io_v18_14 import read
from reorientation_math import quantities,roots
from history_window_math import HORIZONS,orth,eval_at,event_comparison

def csvwrite(file,records):
    if not records:return
    with open(file,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)

def analyze(d,tol):
    x=d['xi_nodes'];A=d['integral_source_frame_11D'][:,:10];g=d['source_integrand_11D'][:,:10];s=float(d['slope']);e=A[80].copy()
    if x.shape!=(641,) or A.shape!=(641,10) or g.shape!=(641,10):raise RuntimeError('Unexpected original checkpoint dimensions')
    if not np.allclose(A[0],0,rtol=0,atol=1e-12):raise RuntimeError('A(0)!=0: source-history reference invalid')
    if not np.allclose(np.diff(x),.0125,atol=1e-10,rtol=0):raise RuntimeError('Unexpected xi grid spacing')
    a=PchipInterpolator(x,A,axis=0);gi=PchipInterpolator(x,g,axis=0)
    mask=np.flatnonzero(x>1+1e-12)
    four=roots(x[mask],quantities(e,A[mask],g[mask],s)['radial_rate'],'pchip')
    if len(four)!=4:raise RuntimeError(f'Full-history event regression: expected 4, got {len(four)}')
    rows=[]; ev=[];errors=[]
    for z in x[mask]:
        full=eval_at(x,A,g,s,e,float(z),ai=a,gi=gi)
        for h in HORIZONS:
            if z-h<0:continue
            r=eval_at(x,A,g,s,e,float(z),h,ai=a,gi=gi)
            old=orth(e,a(z-h));old_excl=old-np.dot(old,r['v'])*r['v']
            signed=float(r['v']@full['v'])
            rows.append({'tolerance':tol,'xi':float(z),'window_xi':h,'signed_cos':signed,
                'line_angle_deg':float(np.degrees(np.arccos(np.clip(abs(signed),0,1)))),
                'old_exclusive_over_full_radius':float(np.linalg.norm(old_excl)/full['rho']),
                'old_radius_over_full_radius':float(np.linalg.norm(old)/full['rho']),
                'recent_radius_over_full_radius':float(r['rho']/full['rho']),
                'recent_turn_speed':r['speed'],'full_turn_speed':full['speed'],
                'projector_difference_fro':float(np.linalg.norm(full['projector']-r['projector'])),
                'recent_radial_rate':r['radial'],'full_radial_rate':full['radial']})
            # Exact additive decomposition of integrals and their moving-boundary derivatives.
            arecent=a(z)-a(z-h);aold=a(z-h)
            frecent=s*(gi(z)-gi(z-h));fold=s*gi(z-h)
            errors.append(max(np.linalg.norm(arecent+aold-a(z)),
                              np.linalg.norm(frecent+fold-s*gi(z))))
    for j,z in enumerate(four,1):
        for h in HORIZONS:
            r=event_comparison(x,A,g,s,e,z,h,ai=a,gi=gi)
            r={'tolerance':tol,'event':j,**r}
            if r['rank3']!=3 or r['relative_rank3_smallest_singular']<=1e-5:
                raise RuntimeError(f'Previous-source rank3 not robust at event {j}, h={h}')
            ev.append(r)
    out={'tolerance':tol,'events_xi':four,'samples':len(mask),'window_lengths_xi':HORIZONS,
        'event_history_metrics':ev,
        'max_additive_history_identity_error':max(errors),
        'min_rank3_smallest_singular_at_events':min(r['relative_rank3_smallest_singular'] for r in ev),
        'window_profiles':{str(h):{
           'n_valid_samples':len([r for r in rows if r['window_xi']==h]),
           'median_line_angle_deg':float(np.median([r['line_angle_deg'] for r in rows if r['window_xi']==h])),
           'median_old_exclusive_over_full_radius':float(np.median([r['old_exclusive_over_full_radius'] for r in rows if r['window_xi']==h]))}
           for h in HORIZONS},
        'scope':'Ablation in archived co-moving source history; not distinct physical trajectories, nor a hysteresis certificate.'}
    return out,rows,ev

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'results');args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    d,prov=read(); reports=[];profile=[];event=[]
    for tol in ('standard','tight'):
        one,rows,ev=analyze(d[tol],tol);reports.append(one);profile+=rows;event+=ev
    aa,bb=reports
    root_diff=float(np.max(abs(np.array(aa['events_xi'])-np.array(bb['events_xi']))))
    event_angle_diff=float(np.max([abs(p['direction_line_angle_deg']-q['direction_line_angle_deg'])
                    for p,q in zip(aa['event_history_metrics'],bb['event_history_metrics'])]))
    if root_diff>1e-5 or event_angle_diff>1e-4:raise RuntimeError('Tolerance pair mismatch')
    summary={'version':'Q550-v18.14-history-window-ablation-transverse-echo',
         'status':'HISTORY_WINDOW_DIRECTIONAL_ECHO_NUMERIC_GAUGE_DECLARED',
         'provenance':prov,'by_tolerance':reports,
         'standard_tight_max_event_xi_drift':root_diff,
         'standard_tight_max_window_event_angle_difference_deg':event_angle_diff,
         'warnings':[
           'This is diagnostic removal of old terms from the original accumulated source, NOT a new nonlinear ODE trajectory.',
           'The old contribution is integrated over [0,xi-h], the recent contribution over [xi-h,xi].',
           'For a sliding window, the differentiated contribution contains BOTH g(xi) and -g(xi-h); neglecting the latter is a wrong kinematic test.',
           'Line angles, norms, speeds and radial events depend on the declared initial Euclidean metric.',
           'Algebraic three-vector independence is invariant under invertible changes of coordinates, when the rank test is numerically well conditioned.',
           'Only four events are sampled; event families do not establish a universal alternation.',
           'No two separately evolved trajectories at matched physical state and phase were compared; physical hysteresis, mechanical constraints and material metric remain uncertified.'
         ]}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    csvwrite(args.output/'history_window_profiles.csv',profile)
    csvwrite(args.output/'event_echo_windows.csv',event)
    report=['# Q550 v18.14 — History-Window Ablation & Transverse Source Echo Audit','',
        '## Premise','',
        'The exact source history A(xi) and integrand g(xi) are retained from the verified Q550 v18.13 genealogy; standard and tight checkpoints verified recursively. No 400-orbit replay.',
        '', 'For h=0.25,0.5,1,2: A_recent=A(xi)-A(xi-h); A_old=A(xi-h). At fixed xi A_total=A_recent+A_old.',
        'Under a moving window A_recent prime = s*(g(xi)-g(xi-h)); A_old prime=s*g(xi-h).',
        'All orientation and scalar comparisons occur in the *declared initial Euclidean 10D frame*.',
        '','## Event-centred echo — angular distance from full history','',
        '| event | xi | h=.25 | h=.50 | h=1.0 | h=2.0 |',
        '|---|---:|---:|---:|---:|---:|']
    for j,z in enumerate(aa['events_xi'],1):
        a=[q for q in aa['event_history_metrics'] if q['event']==j]
        report.append(f"| {j} | {z:.7f} | "+' | '.join(f"{r['direction_line_angle_deg']:.2f} deg" for r in a)+' |')
    report += ['', '## Independent old-history direction outside the recent-history plane','',
        'The vector triple [e, A_recent, A_old] has numerical rank 3 for all four event times and all four windows, using SVD of separately normalized columns. This is a coordinate-invariant *rank* property, not a coordinate-invariant SVD magnitude.',
        f"Smallest event rank-3 singular value in declared gauge: {aa['min_rank3_smallest_singular_at_events']:.6f}.",
        f"Additive history and moving-window derivative identity residual: {aa['max_additive_history_identity_error']:.3g}.",
        f"Max standard/tight event-angle discrepancy: {event_angle_diff:.4g} deg.",
        '', '## Scope and interpretation','',
        'Old-history contributions have distinct orientations at all four events. In this gauge the minima (events 2 and 4) retain a large angle to the complete-history line even with h=2. This is *history support of a cumulative observable*; it is not dynamical hysteresis by itself.',
        'Ablation of A at fixed physical trajectory does NOT model deletion of constitutive memory, because no source perturbation was fed back into the ODE. For physical hysteresis, compare two distinct admissible histories meeting at a common phase-fixed material state and then measure their responses.',
        '', '## Provenance','',
        f"Original v18.13 SHA256: {prov['v18_13_original_sha256']}",
        'All scripts, prior original sendback, summary, profiles, event audit, plots, tests and SHA256 manifest retained.', '']
    (args.output/'REPORT.md').write_text('\n'.join(report))
    print(json.dumps({'status':summary['status'],'events':aa['events_xi'],
         'angles_h2_deg':[q['direction_line_angle_deg'] for q in aa['event_history_metrics'] if q['window_xi']==2.],
         'rank3_min_sigma':aa['min_rank3_smallest_singular_at_events'],
         'maximum_standard_tight_angle_drift_deg':event_angle_diff},indent=2))

if __name__=='__main__':main()
