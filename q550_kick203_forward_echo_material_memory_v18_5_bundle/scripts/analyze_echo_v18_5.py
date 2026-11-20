#!/usr/bin/env python3
"""Postprocess Q550 v18.5 material pulse echo without restarting solver."""
from __future__ import annotations
import json, sys, argparse, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from echo_math import AMPS,DELAYS_XI,READ_AGES_XI,TOLS,signed_even_odd,safe_ratio,unit_angle


def summarize(args):
    args.out.mkdir(parents=True,exist_ok=True)
    rows=[]; raw=[]
    for tol in args.tolerances:
        for eps in args.amps:
            for delay in args.delays:
                p=args.checkpoints/f'{tol}_eps{eps:g}_delay{delay:g}.npz'
                if not p.is_file():raise FileNotFoundError(f'Missing complete experiment: {p}')
                with np.load(p,allow_pickle=False) as z:
                    if not (abs(float(z['eps'][0])-eps)<1e-20 and abs(float(z['delay_xi'][0])-delay)<1e-10):
                        raise RuntimeError(f'checkpoint mismatch {p}')
                    amp=float(z['amp_abs'][0]);v=z['kick_direction_11D'];base=z['baseline_read_11D'];
                    ep=z['plus_echo_read_11D']-base
                    em=z['minus_echo_read_11D']-base
                    wp=z['plus_write_read_11D']-base
                    wm=z['minus_write_read_11D']-base
                    even,odd=signed_even_odd(ep,em)
                    linear_write=0.5*(wp-wm)
                    naivep=z['plus_naive_antikick_immediate_11D']-base[0]
                    naivem=z['minus_naive_antikick_immediate_11D']-base[0]
                    pm=z['plus_linearly_transported_echo_residual_11D']
                    mm=z['minus_linearly_transported_echo_residual_11D']
                    echo0=max(0.5*(np.linalg.norm(ep[0])+np.linalg.norm(em[0])),1e-300)
                    ev0=max(np.linalg.norm(even[0]),1e-300)
                    for i,age in enumerate(z['read_age_xi']):
                        evo=float(np.linalg.norm(even[i]));odo=float(np.linalg.norm(odd[i]))
                        avg_echo=float(0.5*(np.linalg.norm(ep[i])+np.linalg.norm(em[i])))
                        ln=float(np.linalg.norm(linear_write[i]))
                        predicted_even=(pm[i]+mm[i])*.5
                        predicted_odd=(pm[i]-mm[i])*.5
                        e_model=np.linalg.norm(even[i]-predicted_even)
                        o_model=np.linalg.norm(odd[i]-predicted_odd)
                        base_kwargs=dict(tolerance=tol,eps_rel=eps,amplitude_absolute=amp,
                            delay_xi=delay,read_age_xi=float(age),
                            pulse_s=float(z['pulse_s'][0]),read_s=float(z['read_s'][i]),
                            pulse_condition=float(z['pulse_variational_condition'][0]))
                        r={**base_kwargs,
                            'echo_even_norm':evo,'echo_odd_norm':odo,'echo_avg_norm':avg_echo,
                            'raw_write_odd_norm':ln,'echo_fraction_of_write':safe_ratio(avg_echo,ln),
                            'echo_even_fraction_of_write':safe_ratio(evo,ln),
                            'echo_norm_relative_to_immediate':safe_ratio(avg_echo,echo0),
                            'even_relative_to_immediate':safe_ratio(evo,ev0),
                            'echo_even_per_amp_sq':evo/amp**2,
                            'echo_odd_per_amp_cube':odo/amp**3,
                            'echo_even_signed_v':float(even[i]@v),
                            'echo_odd_signed_v':float(odd[i]@v),
                            'even_angle_to_pulse_deg':unit_angle(even[0],even[i]),
                            'even_linear_transport_relative_defect':safe_ratio(e_model,np.linalg.norm(predicted_even)),
                            'odd_linear_transport_relative_defect':safe_ratio(o_model,np.linalg.norm(predicted_odd)),
                            'naive_antikick_avg_norm_at_pulse':float(.5*(np.linalg.norm(naivep)+np.linalg.norm(naivem))),
                            'matched_antikick_avg_norm_at_pulse':echo0,
                            'solver_transfer_condition':float(z['plus_transport_conditions'][i]),
                            'jacobian_fd_relative_error':float(z['jacobian_fd_relative_error'][0])}
                        rows.append(r)
                        raw.append((tol,eps,delay,float(age),even[i].copy(),odd[i].copy()))
    df=pd.DataFrame(rows)
    df.to_csv(args.out/'forward_echo_samples.csv',index=False,float_format='%.17g')
    # Compare the *same* physical target (lag/read) across solver tolerances.
    sens=[]
    if {'standard','tight'}.issubset(args.tolerances):
        key={(tol,eps,lag,age):(even,odd) for tol,eps,lag,age,even,odd in raw}
        for eps in args.amps:
            for lag in args.delays:
                for age in args.read_ages:
                    es,os=key['standard',eps,lag,age]
                    et,ot=key['tight',eps,lag,age]
                    sens.append(dict(eps_rel=eps,delay_xi=lag,read_age_xi=age,
                      even_solver_relative=float(np.linalg.norm(es-et)/max(np.linalg.norm(es),np.linalg.norm(et),1e-300)),
                      odd_solver_relative=float(np.linalg.norm(os-ot)/max(np.linalg.norm(os),np.linalg.norm(ot),1e-300)),
                      even_direction_tolerance_deg=unit_angle(es,et)))
    pd.DataFrame(sens).to_csv(args.out/'tolerance_discrepancy.csv',index=False)
    # Causal scaling exponent of echo even/odd by least-squares in logs; 3 amplitudes.
    fits=[]
    for (tol,delay,age),g in df.groupby(['tolerance','delay_xi','read_age_xi']):
        g=g.sort_values('eps_rel')
        valid=(g.echo_even_norm>0)&(g.echo_odd_norm>0)
        if valid.sum()>=3:
            e_exp=float(np.polyfit(np.log(g.eps_rel),np.log(g.echo_even_norm),1)[0])
            o_exp=float(np.polyfit(np.log(g.eps_rel),np.log(g.echo_odd_norm),1)[0])
        else:e_exp=o_exp=float('nan')
        fits.append(dict(tolerance=tol,delay_xi=delay,read_age_xi=age,
                         even_exponent=e_exp,odd_exponent=o_exp,
                         even_amp2_median=float(np.median(g.echo_even_per_amp_sq)),
                         odd_amp3_median=float(np.median(g.echo_odd_per_amp_cube))))
    fitdf=pd.DataFrame(fits)
    fitdf.to_csv(args.out/'echo_amplitude_scaling.csv',index=False)

    figure_specs=[
       ('echo_fraction_vs_delay.png','echo_even_fraction_of_write','Even nonlinear echo / uncorrected write','log'),
       ('even_per_amp2_vs_delay.png','echo_even_per_amp_sq','Echo even / absolute kick amplitude²','log'),
       ('memory_persistence_vs_age.png','even_relative_to_immediate','Even echo persistence at read age','log')
    ]
    for name,field,yl,yscale in figure_specs:
        fig,ax=plt.subplots(figsize=(9,5))
        use=df[(df.tolerance==args.tolerances[0])]
        if name=='memory_persistence_vs_age.png':
            # One curve per delayed reversal, using middle epsilon.
            use=use[use.eps_rel==args.amps[len(args.amps)//2]]
            for lag,g in use.groupby('delay_xi'):
                ax.plot(g.read_age_xi,g[field],marker='o',label=f'lag {lag:g} xi')
            ax.set_xlabel('Age since anti-kick (xi-equivalent)')
        else:
            use=use[use.read_age_xi==0]
            for eps,g in use.groupby('eps_rel'):
                ax.plot(g.delay_xi,g[field],marker='o',label=f'eps {eps:g}')
            ax.set_xlabel('Delay before anti-kick (xi-equivalent)')
        if yscale=='log':ax.set_yscale('log')
        ax.set_ylabel(yl);ax.set_title(yl)
        ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout()
        fig.savefig(args.out/name,dpi=150);plt.close(fig)
    # Dedicated naive vs compensated one-level scale control.
    fig,ax=plt.subplots(figsize=(9,5))
    show=df[(df.tolerance==args.tolerances[0]) & (df.eps_rel==args.amps[0]) &(df.read_age_xi==0)]
    ax.plot(show.delay_xi,show.naive_antikick_avg_norm_at_pulse,'o-',label='naive same-vector reverse kick')
    ax.plot(show.delay_xi,show.matched_antikick_avg_norm_at_pulse,'o-',label='tangent-matched reverse kick')
    ax.set_yscale('log');ax.set_xlabel('Delay (xi-equivalent)');ax.set_ylabel('Absolute 11D state residual')
    ax.set_title('What remains after compensating material transport?')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(args.out/'naive_vs_tangent_matched.png',dpi=150);plt.close(fig)

    top=df[(df.read_age_xi==0)]
    report={
        'version':'2026-10-16-v18.5-forward-echo',
        'status':'FORWARD_ECHO_MATERIAL_MEMORY_DIAGNOSTIC_OFF_CONSTRAINT',
        'samples':len(df), 'checkpoints':len(args.tolerances)*len(args.amps)*len(args.delays),
        'lag_xi':args.delays,'amplitude_relative':args.amps,
        'read_age_xi':args.read_ages,
        'min_max_immediate_echo_fraction_of_write':[
            float(top.echo_fraction_of_write.min()),float(top.echo_fraction_of_write.max())],
        'median_immediate_echo_fraction_by_lag':{
            str(delay):float(top[top.delay_xi==delay].echo_fraction_of_write.median())
            for delay in args.delays},
        'median_even_exponent':float(fitdf.even_exponent.median()),
        'median_odd_exponent':float(fitdf.odd_exponent.median()),
        'max_even_solver_disagreement':float(max((x['even_solver_relative'] for x in sens),default=float('nan'))),
        'interpretation_boundaries':[
          'Echo uses forward time and exact variational anti-kick; NOT time-reversed dynamics.',
          'Response is physical 11D numerical state, NOT 11D dlocal descriptor.',
          'Without invariants plugin no physical constraint certification.',
          'A nonzero echo is nonlinear path residue in state, NOT proof of thermodynamic or rate-independent hysteresis.',
          'The model is deterministic Markovian in full state: identical full states and identical future forcing imply identical futures.',
          'Zero epsilon would have zero echo; small-amp scaling must be checked against standard/tight errors.'
        ]}
    (args.out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    reportlines=[
      '# Q550 kick203 v18.5 forward echo — measured audit',
      '',
      f'Checkpoints: {report["checkpoints"]}; sampled readouts: {len(df)}.',
      f'Lag xi-equivalent: {args.delays}.',
      f'Relative kick amplitudes: {args.amps}.',
      f'Range of mean echo fractions at reversal: {report["min_max_immediate_echo_fraction_of_write"]}.',
      f'Median even scaling exponent (descriptive only): {report["median_even_exponent"]:.6g}.',
      f'Median odd scaling exponent (often close to numerical floor): {report["median_odd_exponent"]:.6g}.',
      '', '## Immediate nonlinear memory left after tangent-matched anti-kick','',
      '| Delay xi | Median echo / uncancelled write |','|---:|---:|',
    ]
    for lag in args.delays:
        reportlines.append(f'| {lag:g} | {report["median_immediate_echo_fraction_by_lag"][str(lag)]:.6g} |')
    reportlines.extend(['','## Scientific limit','','This is an off-constraint nonlinear **echo / residual state** diagnostic; '
      'the reverse operation is a forward-time anti-kick, not backward-time ODE '
      'integration nor a thermodynamic/control-parameter hysteresis loop. '
      'Control for apparent memory: compare full 11D matched states and solver tolerances.','',''])
    (args.out/'REPORT.md').write_text('\n'.join(reportlines))
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoints',type=Path,default=ROOT/'checkpoint')
    p.add_argument('--out',type=Path,default=ROOT/'results')
    p.add_argument('--amps',nargs='+',type=float,default=list(AMPS))
    p.add_argument('--delays',nargs='+',type=float,default=list(DELAYS_XI))
    p.add_argument('--read-ages',nargs='+',type=float,default=list(READ_AGES_XI))
    p.add_argument('--tolerances',nargs='+',default=list(TOLS))
    summarize(p.parse_args())
