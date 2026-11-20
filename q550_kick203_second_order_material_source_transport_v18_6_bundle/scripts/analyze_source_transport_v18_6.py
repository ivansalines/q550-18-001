#!/usr/bin/env python3
"""No-solver v18.6 report and figures from checkpoint and v18.5 observations."""
from __future__ import annotations
import io,json,sys,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from source_transport_math import angle_deg,relative_error,signed_even,slope_finite_amplitude

def main():
    p=ROOT/'results';p.mkdir(exist_ok=True)
    files=list((ROOT/'checkpoint').glob('*_source_transport.npz'))
    if not files:raise RuntimeError('Q550 solver checkpoints absent; execute run_source_transport_v18_6.py on server first')
    comp=[];bins=[];sources=[]
    by_tol={}
    for f in files:
        tol=f.name.split('_')[0]
        with np.load(f) as z: d={k:z[k].copy() for k in z.files}
        by_tol[tol]=d
        comp.extend(json.loads((ROOT/'checkpoint'/f'{tol}_observed_comparison.json').read_text()))
        bins11=d['source_bin_contributions_11D']; total=d['predicted_quadratic_at_delays_11D']
        for k,lag in enumerate((.25,.5,1.,2.,4.,8.)):
            for j,(left,right) in enumerate(zip(d['source_bin_edges'][:-1],d['source_bin_edges'][1:])):
                vec=bins11[k,j]
                bins.append(dict(tolerance=tol,delay_xi=lag,source_from_xi=float(left),source_to_xi=float(right),
                    contribution_norm=float(np.linalg.norm(vec)),
                    total_norm=float(np.linalg.norm(total[k])),
                    fraction_of_total_norm=float(np.linalg.norm(vec)/np.linalg.norm(total[k])) if np.linalg.norm(total[k]) else np.nan,
                    signed_along_total=float(vec@total[k]/np.linalg.norm(total[k])) if np.linalg.norm(total[k]) else np.nan,
                    angle_to_total_deg=angle_deg(vec,total[k])))
        for xi,sv,cv in zip(d['xi_nodes'],d['source_11D'],d['cocycle_condition_nodes']):
            sources.append(dict(tolerance=tol,age_xi=float(xi),source_norm=float(np.linalg.norm(sv)),
                source_projection_initial=float(sv@d['kick_direction_11D']),cocycle_condition=float(cv)))
    C=pd.DataFrame(comp).sort_values(['tolerance','epsilon_relative','delay_xi','read_age_xi'])
    B=pd.DataFrame(bins);S=pd.DataFrame(sources)
    C.to_csv(p/'echo_vs_second_variation.csv',index=False)
    B.to_csv(p/'source_age_bins.csv',index=False)
    S.to_csv(p/'source_profile.csv',index=False)
    # Tolerance convergence of the predictions and full 11D vectors (not norms only).
    tols=[]
    if 'standard' in by_tol and 'tight' in by_tol:
        st=by_tol['standard'];ti=by_tol['tight']
        for i,lag in enumerate((.25,.5,1.,2.,4.,8.)):
            a=st['predicted_quadratic_at_delays_11D'][i];b=ti['predicted_quadratic_at_delays_11D'][i]
            tols.append(dict(delay_xi=lag,relative_vector_difference=relative_error(a,b),
                angle_deg=angle_deg(a,b)))
    pd.DataFrame(tols).to_csv(p/'source_tolerance_comparison.csv',index=False)
    # Independently extrapolate v18.5 material echo to zero amplitude at each delay.
    from run_source_transport_v18_6 import AMPS,DELAYS,AGES,read_v185,verify_history
    history=verify_history(); extrap=[]
    for tol in by_tol:
        pred=by_tol[tol]['predicted_quadratic_at_delays_11D']
        for i,lag in enumerate(DELAYS):
            amps=[];vecs=[]
            for eps in AMPS:
                q=read_v185(history,tol,eps,lag)
                a=float(q['amp_abs'][0]);amps.append(a)
                v=signed_even(q['plus_echo_read_11D'][0],q['minus_echo_read_11D'][0],q['baseline_pulse_11D'])/a**2
                vecs.append(v)
            extrapv=slope_finite_amplitude(vecs,amps)
            extrap.append(dict(tolerance=tol,delay_xi=lag,
                extrapolated_quadratic_norm=float(np.linalg.norm(extrapv)),
                predicted_quadratic_norm=float(np.linalg.norm(pred[i])),
                relative_difference=relative_error(pred[i],extrapv),
                angle_deg=angle_deg(pred[i],extrapv)))
    E=pd.DataFrame(extrap);E.to_csv(p/'extrapolated_zero_amplitude.csv',index=False)
    imm=C[C.read_age_xi==0]
    summary=dict(version='2026-10-16-q550-v18.6',status='Q550 source reconstructed',
        counts=dict(physical_solver_tolerances=len(by_tol),echo_comparisons=int(len(C)),
            source_bin_rows=int(len(B)),source_nodes_per_tol=int(len(S)//len(by_tol))),
        median_immediate_relative_vector_error=float(imm.relative_vector_error.median()),
        max_immediate_relative_vector_error=float(imm.relative_vector_error.max()),
        median_all_read_relative_vector_error=float(C.relative_vector_error.median()),
        median_hessian_step_relative_error=float(imm.hessian_step_relative_error.median()),
        median_simpson_relative_error=float(imm.simpson_grid_relative_error.median()),
        median_extrapolated_rel_error=float(E.relative_difference.median()),
        max_source_cocycle_condition=float(S.cocycle_condition.max()),
        reliability_note='All numerical gates descriptive; no invariant-manifold certificate')
    (p/'summary.json').write_text(json.dumps(summary,indent=2))
    report=['# Q550 v18.6 — Material quadratic source transport','',
       'Actual Q550 data. Off-constraint diagnostic. The analytic 11D Jacobian was differenced to obtain the local second directional derivative.','',
       '## Global metrics','',f'- Median immediate relative **11D vector** prediction error: {summary["median_immediate_relative_vector_error"]:.5g}',
       f'- Max immediate relative 11D vector prediction error: {summary["max_immediate_relative_vector_error"]:.5g}',
       f'- Median Simpson coarse/fine error: {summary["median_simpson_relative_error"]:.5g}',
       f'- Median Hessian step/half-step error: {summary["median_hessian_step_relative_error"]:.5g}',
       f'- Median extrapolated zero-amplitude comparison error: {summary["median_extrapolated_rel_error"]:.5g}',
       '', '## Errors by delay and tolerance','',
       '| Tolerance | Delay xi | Median immediate vector error |','|:--|--:|--:|']
    for (tol,lag),sub in imm.groupby(['tolerance','delay_xi']):
       report.append(f'| {tol} | {lag:g} | {sub.relative_vector_error.median():.6g} |')
    report+=['','## Interpretation','',
       '- The source is S2(s)=0.5*D²f(y0(s))[P(s)v,P(s)v]. Its transported age-bin contributions form a signed 11D vector sum, not a sum of positive memory percentages.',
       '- Increasing K2 magnitude with delay does not imply monotonically accumulating source; source contributions can rotate and cancel after propagation.',
       '- All material kicks remain off the certified constraint manifold; no independent hysteresis mechanism is inferred.',
       '- Use the CSV/NPZ for directional sign, source-birth bins, and error checks.']
    (p/'REPORT.md').write_text('\n'.join(report)+'\n')
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8.5,5.0))
    for tol,g in C[(C.read_age_xi==0)&(C.epsilon_relative==3e-8)].groupby('tolerance'):
        ax.plot(g.delay_xi,g.observed_quadratic_norm,marker='o',label=f'{tol}: observed')
        ax.plot(g.delay_xi,g.predicted_quadratic_norm,linestyle='--',marker='x',label=f'{tol}: source integral')
    ax.set(xlabel='Time since first kick (xi units)',ylabel='Quadratic coefficient norm (state units⁻¹)',title='Observed echo versus transported material second variation')
    ax.legend(fontsize=8);ax.grid(alpha=.3);fig.tight_layout();fig.savefig(p/'source_prediction_vs_echo.png',dpi=165);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,4.8))
    for tol,g in imm.groupby('tolerance'):
        e=g.groupby('delay_xi').relative_vector_error.median()
        ax.semilogy(e.index,e.values,marker='o',label=tol)
    ax.set(xlabel='Delay (xi units)',ylabel='Median 11D relative error',title='Second-order material source closure error');ax.legend();ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(p/'closure_error_by_delay.png',dpi=165);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,5.0))
    a=B[(B.tolerance=='tight')&(B.delay_xi==8)].sort_values('source_from_xi')
    if not len(a):a=B[B.delay_xi==8].sort_values('source_from_xi').iloc[:6]
    labels=[f'{l:g}–{r:g}' for l,r in zip(a.source_from_xi,a.source_to_xi)]
    x=np.arange(len(a));ax.bar(x,a.signed_along_total.values)
    ax.set_xticks(x,labels);ax.set(xlabel='Source birth interval (xi units)',ylabel='Projection on final K2 (state units⁻¹)',title='Where the quadratic source was born (transported to xi=8)')
    ax.axhline(0,color='black',lw=.7);fig.tight_layout();fig.savefig(p/'source_birth_contributions.png',dpi=165);plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
