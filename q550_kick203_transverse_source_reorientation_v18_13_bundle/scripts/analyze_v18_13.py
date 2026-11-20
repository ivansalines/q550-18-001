#!/usr/bin/env python3
"""Q550 v18.13 — residual transverse-source orientation at the four radial turnings.

Reads frozen v18.12 provenance, independently computes geometric quantities, compares
polar source and line-turning angles around events and subjects event detection to
interpolation, tolerance and observer-metric sensitivity checks.
"""
from __future__ import annotations
import sys,argparse,csv,json,hashlib
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from checkpoint_io_v18_13 import read
from reorientation_math import quantities,roots,event_at,root_pair_coherence,span_three_relative_volume
WINDOWS=(0.05,0.10,0.25,0.50)

def dumpcsv(path,rows):
    if not rows:return
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def geom_data(d,T=None):
    x=d['xi_nodes']; inds=np.flatnonzero(x>1.0+1e-12)
    e=d['integral_source_frame_11D'][80,:10].copy()
    A=d['integral_source_frame_11D'][inds,:10].copy()
    g=d['source_integrand_11D'][inds,:10].copy()
    s=float(d['slope']); xx=x[inds]
    if T is not None:
        e=T@e;A=A@T.T;g=g@T.T
    q=quantities(e,A,g,s)
    return xx,e,A,g,s,q

def metrics_audit(d):
    diag=np.linspace(-1.,1.,10)
    gauges={
      'identity':np.eye(10),
      'uniform_scale_1p7':np.eye(10)*1.7,
      'mild_diagonal':np.diag(np.exp(0.10*diag)),
      'strong_diagonal':np.diag(np.exp(0.70*diag))}
    out={}
    for name,T in gauges.items():
        x,e,A,g,s,q=geom_data(d,T)
        out[name]={'turnings_xi':roots(x,q['radial_rate'],'pchip'),
           'turn_rate_min':float(q['speed'].min()),
           'turn_rate_median':float(np.median(q['speed'])),
           'positive_turn_speed_samples':int(np.sum(q['speed']>1e-7)),
           'metric_scope':'Euclidean after the stated fixed GL(10) map, not a derived physical metric'}
    return out

def audit(d,name):
    x,e,A,g,s,q=geom_data(d)
    # Fundamental differentiation identities and consistency with archived residual-plane CSV.
    fd=np.gradient(q['rho'],x,edge_order=2)
    radial_fd=float(np.median(np.abs(fd[2:-2]-q['radial_rate'][2:-2])/
                  np.maximum(np.abs(q['radial_rate'][2:-2]),1.)))
    turn4=roots(x,q['radial_rate'],'linear')
    p4=roots(x,q['radial_rate'],'pchip')
    c4=roots(x,q['radial_rate'],'cubic')
    # Downsampling by 2: phase-locked parity provides an independent discretization stress test.
    even=roots(x[::2],q['radial_rate'][::2],'pchip')
    odd=roots(x[1::2],q['radial_rate'][1::2],'pchip')
    if len(p4)!=4 or len(turn4)!=4 or len(c4)!=4 or len(even)!=4 or len(odd)!=4:
        raise RuntimeError('Four turning events not reproduced across sampling/interpolation')
    root_method_spread=max(abs(np.array(p4)-np.array(turn4)).max(),abs(np.array(p4)-np.array(c4)).max())
    root_decimation_spread=max(abs(np.array(p4)-np.array(even)).max(),abs(np.array(p4)-np.array(odd)).max())
    records=[];wid=[]
    for i,z in enumerate(p4):
        p=event_at(x,e,A,g,s,z)
        if float(p['speed'])<1e-7: raise RuntimeError('Vanishing transverse turn at radial event')
        # Values at roots are interpolant estimates, not certified exact zeros.
        event={'tolerance':name,'event':i+1,'xi_pchip':z,
          'xi_linear':turn4[i],'xi_cubic':c4[i],
          'xi_even_decimated':even[i],'xi_odd_decimated':odd[i],
          'radial_rate_at_root_interp':float(p['radial_rate']),
          'rho_at_root':float(p['rho']),
          'turn_speed_at_root':float(p['speed']),
          'source_outside_fraction_at_root':float(p['fraction_outside']),
          'algebraic_noncollinearity_relative_3volume_gauge':span_three_relative_volume(e,PchipInterpolator(x,A,axis=0)(z),PchipInterpolator(x,g,axis=0)(z))}
        local=[]
        for dt in WINDOWS:
            z0=z-dt;z1=z+dt
            if z0<=x[0] or z1>=x[-1]:raise RuntimeError('Event comparison exceeds data interval')
            vv=root_pair_coherence(x,e,A,g,s,z,dt)
            item={'tolerance':name,'event':i+1,'xi_pchip':z,'halfwidth_xi':dt,**vv}
            wid.append(item);local.append(item)
        event['radial_sign_left_at_0p05']=local[0]['left_radial_rate']
        event['radial_sign_right_at_0p05']=local[0]['right_radial_rate']
        event['coherence_at_0p05']=local[0]['cosine']
        event['coherence_at_0p25']=local[2]['cosine']
        event['coherence_at_0p5']=local[3]['cosine']
        records.append(event)
    out={'tolerance':name,'samples':len(x),'range_xi':[float(x[0]),float(x[-1])],
         'four_pchip_radial_events_xi':p4,
         'max_interpolator_root_difference_xi':float(root_method_spread),
         'max_even_odd_decimation_root_difference_xi':float(root_decimation_spread),
         'radial_derivative_fd_median_relative':radial_fd,
         'speed_positive_count':int(np.sum(q['speed']>1e-7)),
         'max_abs_radial_rate':float(np.max(abs(q['radial_rate']))),
         'min_speed':float(np.min(q['speed'])),
         'median_speed':float(np.median(q['speed'])),
         'min_fraction_source_outside_plane':float(q['fraction_outside'].min()),
         'events':records,'observer_metrics':metrics_audit(d),
         'interpretation':'Radial turning and angle coherence are gauge/metric dependent. Nonzero algebraic wedge is invertible-coordinate invariant. No mechanical torsion, physical metric or hysteresis proven.'}
    profile=[]
    for j,xi in enumerate(x):
        profile.append({'tolerance':name,'xi':float(xi),'radial_rate':float(q['radial_rate'][j]),
           'turn_speed':float(q['speed'][j]),'rho':float(q['rho'][j]),
           'source_outside_fraction':float(q['fraction_outside'][j])})
    return out,records,wid,profile

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default=str(ROOT/'results'));args=parser.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    data,provenance=read()
    reps=[];events=[];windows=[];prof=[]
    for name in ('standard','tight'):
        r,ev,ws,pp=audit(data[name],name);reps.append(r);events+=ev;windows+=ws;prof+=pp
    a,b=reps
    drift=float(np.max(np.abs(np.array(a['four_pchip_radial_events_xi'])-np.array(b['four_pchip_radial_events_xi']))))
    coh_drift=float(np.max(abs(np.array([w['cosine'] for w in windows if w['tolerance']=='standard'])-
                                np.array([w['cosine'] for w in windows if w['tolerance']=='tight']))))
    if drift>1e-5 or coh_drift>1e-5: raise RuntimeError('Standard/tight event or orientation mismatch')
    # The four reorientations need not be identical: coarse span is longer than local window.
    if not all(ev['coherence_at_0p05']>.9 for ev in a['events']):
        raise RuntimeError('Expected local directional continuity absent')
    summary={'version':'Q550-v18.13-transverse-source-reorientation','provenance':provenance,
      'status':'NUMERICAL_LOCAL_ORIENTATION_PERSISTENCE_GAUGE_DECLARED',
      'by_tolerance':reps,
      'standard_tight_max_event_xi_diff':drift,
      'standard_tight_max_direction_coherence_diff':coh_drift,
      'limitations':[
        'Radial turning counts, positions, angles and direction-coherence cosines depend on the chosen metric and initial gauge.',
        'Wedge rank under invertible coordinate changes is invariant; normalized three-volumes and their magnitudes are not.',
        'Interpolation and decimation are numerical stability evidence, not interval certificates between samples.',
        'At radial turning the source is tangential to rho-levels by construction, so source_outside_fraction~=1 is an identity, not an independent discovery.',
        'Nonzero residual turning is not by itself a proof of nonholonomic torsion, material metric, or physical hysteresis.',
        'No two distinct histories brought to the same phase-fixed material state were compared.'
      ]}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    dumpcsv(out/'event_turnings.csv',events)
    dumpcsv(out/'event_orientation_windows.csv',windows)
    dumpcsv(out/'source_components_profile.csv',prof)
    report=['# Q550 v18.13 — Transverse Source Reorientation & Radial Turning Audit','',
      '## Reproduction',
      'Original v18.12 sendback independently checked using CRC and SHA256. Data: standard/tight; 560 nodes each; previous genealogical chain retained inside original nested archive.',
      '', '## Four radial crossings and directional coherence','',
      '| event | xi (standard, PCHIP) | local cos ±0.05 | broad cos ±0.25 | broad cos ±0.50 |',
      '|---|---:|---:|---:|---:|']
    for ev in a['events']:
      report.append(f"| {ev['event']} | {ev['xi_pchip']:.9f} | {ev['coherence_at_0p05']:+.6f} | {ev['coherence_at_0p25']:+.6f} | {ev['coherence_at_0p5']:+.6f} |")
    report += ['',f'- Maximum standard/tight root drift: {drift:.3g} xi.',
      f"- Maximum linear/cubic/PCHIP method drift: {a['max_interpolator_root_difference_xi']:.3g} xi.",
      f"- Maximum even/odd decimation root drift: {a['max_even_odd_decimation_root_difference_xi']:.3g} xi.",
      f'- Maximum standard/tight coherence drift: {coh_drift:.3g}.',
      '- Continuous tangent directions were compared using spherical Levi-Civita parallel transport in a **declared Euclidean initial 10D gauge**.',
      '', '## Metric sensitivity',
      '| fixed observer | number of radial turnings | locations (xi) |','|---|---:|---|']
    for tag,value in a['observer_metrics'].items():
      report.append(f"| {tag} | {len(value['turnings_xi'])} | {', '.join(f'{v:.5f}' for v in value['turnings_xi'])} |")
    report += ['', '## Interpretation and hard limits',
      'At the four radial turnings, rho_prime passes through zero while residual line turning remains strictly positive at sampled/interpolated points. The transverse fraction approaches one at each event algebraically, since radial source projection is zero there; this is **not** separate physical evidence.',
      'Local tangent-direction coherence is strong at ±0.05 xi, while ±0.25 and ±0.50 xi resolve progressive reorientation, including sign reversal for some events at the broadest window.',
      'The effect is a property of the co-moving source-frame residual under the chosen Euclidean metric. True coordinate-free claims are limited to rank/nonvanishing of the appropriate wedge and exact tensorial equations. No distinguished material metric, physical torsion, constraint, or hysteresis loop has been established.',
      '', '## Next decisive experiment',
      'To claim history dependence beyond source accumulation, construct distinct admissible histories sharing the same reduced material state and phase (with verified tolerance) and measure unequal residual response, or derive a material metric from the actual ODE constitutive structure.', '']
    (out/'REPORT.md').write_text('\n'.join(report))
    print(json.dumps({'status':summary['status'],'events':a['four_pchip_radial_events_xi'],
                      'coherence_pm0p05':[ev['coherence_at_0p05'] for ev in a['events']],
                      'coherence_pm0p50':[ev['coherence_at_0p5'] for ev in a['events']],
                      'standard_tight_max_event_drift':drift},indent=2))
if __name__=='__main__':main()
