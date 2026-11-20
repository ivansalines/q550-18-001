#!/usr/bin/env python3
"""Summarize the dense Q550 source-birth sign crossing and plot 11D rates."""
from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'

def main():
    summaries=[];data={}
    for tol in ('standard','tight'):
        path=OUT/f'{tol}_crossing_summary.json';p=OUT/f'{tol}_dense_projection.csv'
        if not path.exists() or not p.exists():raise FileNotFoundError(f'Missing {tol}: run Q550 v18.7 on server')
        summaries.append(json.loads(path.read_text()))
        with p.open() as f:data[tol]=list(csv.DictReader(f))
    a,b=summaries
    if len(a['roots'])!=len(b['roots']):raise RuntimeError('Different number of roots between tolerances')
    root_concordance=[{'index':i+1,'standard_xi':a['roots'][i]['xi_root'],
      'tight_xi':b['roots'][i]['xi_root'],
      'difference_xi':b['roots'][i]['xi_root']-a['roots'][i]['xi_root']} for i in range(len(a['roots']))]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,5))
    for tol in ('standard','tight'):
        rows=data[tol];x=[float(r['xi']) for r in rows];y=[float(r['old_signed_percent']) for r in rows]
        ax.plot(x,y,label=tol,lw=1.4 if tol=='standard' else 1.,linestyle='-' if tol=='standard' else '--')
    ax.axhline(0,color='black',lw=.8)
    for r in a['roots']:ax.axvline(r['xi_root'],lw=.7,color='gray',ls=':',label='Crossing' if r['crossing_index']==1 else None)
    ax.set(xlabel='Delay from kick 203 (xi)',ylabel='Old 0–1 xi signed contribution (%)',title='v18.7: dense full-11D signed source-birth crossing')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(OUT/'source_birth_sign_crossing_dense.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    for tol in ('standard','tight'):
        rows=data[tol];x=[float(r['xi']) for r in rows]
        ax.plot(x,[float(r['dot_rate_metric']) for r in rows],label=f'{tol}: cocycle deformation',lw=1.3)
        ax.plot(x,[float(r['dot_rate_new_source']) for r in rows],label=f'{tol}: new source birth',ls='--',lw=1.3)
    ax.axhline(0,lw=.7,color='black')
    for r in a['roots']:ax.axvline(r['xi_root'],lw=.7,color='gray',ls=':')
    ax.set(xlabel='Delay from kick 203 (xi)',ylabel='d(E·K2)/dξ contribution',title='v18.7: sign-crossing source vs metric-rate mechanism')
    ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(OUT/'crossing_rate_source_vs_transport.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    r=data['tight'];x=[float(t['xi']) for t in r]
    ax.plot(x,[float(t['old_vs_at_one_lab_angle_degrees']) for t in r],label='Old vector: 11D laboratory angle since xi=1')
    ax.plot(x,[float(t['old_vs_at_one_polar_corot_angle_degrees']) for t in r],label='Old vector: polar co-rotating angle since xi=1')
    ax.plot(x,[float(t['angle_old_total_degrees']) for t in r],label='Old vs total (frame-invariant for orthogonal frames)')
    ax.plot(x,[float(t['original_source_frame_early_total_angle_degrees']) for t in r],ls=':',label='Source-frame old vs total angle (not isometric)')
    ax.axhline(90,ls='--',color='gray',lw=.8)
    for c in a['roots']:ax.axvline(c['xi_root'],lw=.7,color='gray',ls=':')
    ax.set(xlabel='Delay from kick 203 (xi)',ylabel='Angle (degrees)',title='v18.7: lab, polar co-moving, and inverse-cocycle angles')
    ax.legend(fontsize=7);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(OUT/'polar_corotation_and_source_frame_angles.png',dpi=180);plt.close(fig)
    status={'version':'2026-10-16-q550-v18.7-material-birth-crossing',
            'root_concordance':root_concordance,
            'num_dense_rows':{t:len(data[t]) for t in data},
            'max_v186_anchor_relative_error':max(a['full_11D_anchor_max_relative_error'],b['full_11D_anchor_max_relative_error']),
            'source':'original frozen v18.6 source and Q550 11D variational replay',
            'mode':'OFF_CONSTRAINT_DIAGNOSTIC',
            'interpretation':'Neither a material rigid-body flip nor constitutive hysteresis is certified.'}
    (OUT/'summary.json').write_text(json.dumps(status,indent=2)+'\n')
    report=['# Q550 v18.7 — Source-Birth Sign Crossing & Co-Moving Transport', '',
      '**Constraint status: OFF_CONSTRAINT_DIAGNOSTIC.** The Euclidean 11D metric is model-coordinate dependent and not a certified physical rigid-body metric.', '',
      '## Definition', '',
      'Early source `E=P(s)*a_early`, total quadratic material source `K2=P(s)*a_total`; roots satisfy `E·K2=0`, with `a_early=int_0^1 P^-1 S2 ds` and `a_total=int_0^xi P^-1 S2 ds`.', '',
      '`d_xi(E·K2) = a_early^T (d_xi G) a_total + a_early^T G (d_xi a_total)`, where `G=P^T P`.', '',
      'The first term is tangent cocycle deformation (including stretch/shear); the second is newly born source. Neither is a phenomenological memory kernel.', '',
      '## Crossings and numerical controls', '']
    for s in summaries:
        report += [f'### {s["tolerance"]}', '',f'- Crossings found: {len(s["roots"])}',f'- Worst v18.6 11D regression error: {s["full_11D_anchor_max_relative_error"]:.3e}',f'- Max cocycle condition: {s["max_cocycle_condition"]:.4f}',f'- Largest old-source inverse-cocycle drift: {s["old_co_moving_inverse_drift_max"]:.3e}', '']
        for i,r in enumerate(s['roots']):
            co=s['root_controls'][i]
            report += [f'- Root {i+1}: xi={r["xi_root"]:.10f}; material source s={r["absolute_s"]:.12g}; angle old/total={r["old_total_angle_deg"]:.6f} deg; signs at ±0.05 xi: {r["old_signed_pct_left_0_05"]:.6g}% / {r["old_signed_pct_right_0_05"]:.6g}%.',
             f'  - d(E·K2)/dxi: cocycle {r["dot_rate_metric"]:+.9g}, new source {r["dot_rate_source"]:+.9g}, total {r["dot_rate_sum"]:+.9g}; finite-difference {r["dot_rate_fd"]:+.9g}.',
             f'  - Root difference with Hessian half step {co["half_hessian_difference_xi"]:+.4e} xi; with coarse Simpson {co["coarse_quadrature_difference_xi"]:+.4e} xi.', '']
    if root_concordance:
        report += ['## Tolerance agreement','']
        for r in root_concordance:report.append(f'- Root {r["index"]}: tight–standard = {r["difference_xi"]:+.6e} xi')
        report.append('')
    report += ['## Interpretation limit', '',
      'The early vector is exactly fixed in the inverse fundamental-matrix coordinates up to floating-point error. The polar co-rotating frame removes one canonical Euclidean rotation of the complete cocycle, **not** a mechanically unique body rotation. The crossing of signed 11D projection is a geometric event of the chosen observable, and does not prove a physical torsional flip or hysteresis.', '',
      'The source data are only sampled every 0.0125 xi. Root interpolation refines the location inside a sampled sign change but does not exclude unresolved pairs of crossings between adjacent nodes. The Simpson and Hessian/tolerance comparisons quantify discretization sensitivity; they are not interval enclosures.', '']
    (OUT/'REPORT.md').write_text('\n'.join(report))
    print(json.dumps(status,indent=2))

if __name__=='__main__':main()
