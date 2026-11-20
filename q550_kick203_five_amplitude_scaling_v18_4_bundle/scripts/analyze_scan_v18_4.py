#!/usr/bin/env python3
"""v18.4 five-amplitude frozen-11D scaling analysis from archived + new histories.

The full geometric descriptor and its exact factor-secants are evaluated in
one frozen 11D coordinate frame. The material 11D coordinate remains separate.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
from scaling_math import AMPS,EPS_UNIT,fit_odd,safe_median,finite_or_none
from run_scan_v18_4 import read_bridge,verify_source,NEW_AMPS


def load_data(checkpoints):
    verify_source()
    old=read_bridge(ROOT/'inputs/q550_kick203_material_to_observer_bridge_v18_3_bundle.zip')
    result={}
    for tol in ('standard','tight'):
        for eps in AMPS:
            key=f'{tol}__eps{eps:g}'
            if float(eps) in NEW_AMPS:
                path=checkpoints/f'{tol}_eps{eps:g}.npz'
                if not path.exists():
                    raise FileNotFoundError(f'{path} missing; run scripts/run_scan_v18_4.py on the Q550 server first')
                with np.load(path,allow_pickle=False) as z:
                    a={k:np.asarray(z[k]) for k in z.files}
                if not np.array_equal(a['xi'],old['xi']):raise RuntimeError('new xi grid drift')
                if not np.allclose(a['axis'],old['standard__eps1e-08__frozen_axis_11D'],atol=1e-14,rtol=0):
                    raise RuntimeError('frozen axis drift')
                result[key]={
                    'delta':a['centered_delta_11D'],
                    'inter':a['centered_interaction_11D'],
                    'sec':a['centered_factor_secants_3x11'],
                    'order2':a['non_straddling_order2_delta_11D'],
                    'materialA':a['material_difference_A_11D']/(2*float(a['source_amp_absolute'][0])),
                    'materialB':a['material_difference_B_11D']/(2*float(a['source_amp_absolute'][0])),
                    'source':'new Q550 nonlinear ODE +/- histories',
                }
            else:
                result[key]={
                    'delta':old[key+'__centered_delta_11D'],
                    'inter':old[key+'__centered_interaction_11D'],
                    'sec':old[key+'__centered_factor_secants_3x11'],
                    'order2':old[key+'__order2_delta_11D'],
                    'materialA':old[key+'__material_uA_11D'],
                    'materialB':old[key+'__material_uB_11D'],
                    'source':'pinned v18.3 historical measured level',
                }
    return old,result


def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    old, exp=load_data(args.checkpoints)
    xi=old['xi'];cross=float(old['crossing_xi'][0]);sB=old['sB'];s0=float(old['s0'][0])
    clean=old['centered_stencil_clean'].astype(bool)
    event=np.abs(sB-s0)<5e-12
    post=clean&(sB>s0);pre=clean&(sB<s0)
    if (len(xi),int(post.sum()),int(pre.sum()))!=(47,19,17):
        raise RuntimeError(f'Expected frozen 47/19/17 geometry, got {len(xi)}/{post.sum()}/{pre.sum()}')
    axis=np.asarray(old['standard__eps1e-08__frozen_axis_11D'],float)
    axis/=np.linalg.norm(axis)
    samp=[];fits=[];summary={
        'version':'2026-10-16-q550-kick203-v18.4-five-amplitude-scaling',
        'amplitudes_relative':AMPS.tolist(), 'reference_amplitude':EPS_UNIT,
        'tolerances':['standard','tight'], 'xi_B_crossing':cross,
        'clean_after_B_count':int(post.sum()),'clean_before_B_count':int(pre.sum()),
        'common_axis_source':'v18.3 standard epsilon 1e-8 frozen A; fixed for all amplitudes/tolerances',
        'constraint_certificate':'NOT_CERTIFIED; histories solve the ODE only',
        'interpretation':'I_ABC is a descriptor-factor mixed secant, not independently measured physical hysteresis',
    }
    max_close=0.;min_preB_material=0.
    for tol in ('standard','tight'):
        delta=np.stack([exp[f'{tol}__eps{eps:g}']['delta'] for eps in AMPS])
        inter=np.stack([exp[f'{tol}__eps{eps:g}']['inter'] for eps in AMPS])
        sec=np.stack([exp[f'{tol}__eps{eps:g}']['sec'] for eps in AMPS])
        ord2=np.stack([exp[f'{tol}__eps{eps:g}']['order2'] for eps in AMPS])
        massA=np.stack([exp[f'{tol}__eps{eps:g}']['materialA'] for eps in AMPS])
        massB=np.stack([exp[f'{tol}__eps{eps:g}']['materialB'] for eps in AMPS])
        shape=(len(AMPS),len(xi),11)
        if any(a.shape!=shape for a in (delta,inter,ord2,massA,massB)) or sec.shape!=(5,47,3,11):
            raise RuntimeError('experiment shape mismatch')
        max_close=max(max_close,float(np.max(np.linalg.norm(delta-(sec.sum(axis=2)+inter),axis=2))))
        min_preB_material=max(min_preB_material,float(np.max(np.linalg.norm(massB[:,sB<s0,:],axis=2))))
        fit_d=fit_odd(AMPS,delta,(1,3,5))
        fit_i3=fit_odd(AMPS,inter,(3,))
        fit_i35=fit_odd(AMPS,inter,(3,5))
        # The one-sided second-order observer is independently fit for all
        # non-event samples; centered clean post-B nodes agree on physical xi.
        orderfit=fit_odd(AMPS,ord2[:,~event],(1,3,5))
        s=delta@axis;i=inter@axis;a=sec.sum(axis=2)@axis
        ratio_vec=inter/((AMPS/EPS_UNIT)[:,None,None]**3)
        ray_norm=np.linalg.norm(ratio_vec,axis=2)
        for k,eps in enumerate(AMPS):
            for j,x in enumerate(xi):
                samp.append(dict(tolerance=tol,eps_rel=float(eps),eps_over_1e8=float(eps/EPS_UNIT),
                  data_origin=exp[f'{tol}__eps{eps:g}']['source'],xi=float(x),xi_minus_B=float(x-cross),
                  clean_before_B=bool(pre[j]),clean_after_B=bool(post[j]),event=bool(event[j]),
                  delta_norm=float(np.linalg.norm(delta[k,j])),interaction_norm=float(np.linalg.norm(inter[k,j])),
                  interaction_over_eps_cube_norm=float(ray_norm[k,j]),
                  signed_delta=float(s[k,j]),signed_interaction=float(i[k,j]),
                  signed_single_factor_sum=float(a[k,j]),
                  signed_order2=float(ord2[k,j]@axis) if not event[j] else np.nan,
                  material_A_sensitivity_norm=float(np.linalg.norm(massA[k,j])),
                  material_B_sensitivity_norm=float(np.linalg.norm(massB[k,j])),
                  closure_norm=float(np.linalg.norm(delta[k,j]-(sec[k,j].sum(axis=0)+inter[k,j])))))
        for j,x in enumerate(xi):
            # Relative tests of functional forms at fixed physical support.
            ino=float(np.linalg.norm(inter[-1,j])); deno=float(np.linalg.norm(inter[1,j]))
            eff=np.log(ino/deno)/np.log(3.) if deno>1e-24 and ino>0 else np.nan
            # Fifth-vs-third order effect at the largest amplitude: use the
            # two-term fit; interpret only if solver resolution supports it.
            i3=fit_i35['coefficients'][0,j];i5=fit_i35['coefficients'][1,j]
            cubic_norm=float(np.linalg.norm(i3)*(3.**3)); quintic_norm=float(np.linalg.norm(i5)*(3.**5))
            r=ray_norm[:,j]; rel_range=(float(np.ptp(r))/max(float(np.median(r)),1e-300))
            fits.append(dict(tolerance=tol,xi=float(x),xi_minus_B=float(x-cross),
                clean_after_B=bool(post[j]),clean_before_B=bool(pre[j]),
                power_1e8_to_3e8=float(eff),
                scaled_cubic_coefficient_relative_range=rel_range,
                cubic_coefficient_norm=float(np.linalg.norm(i3)),
                quintic_to_cubic_at_3e8=float(quintic_norm/max(cubic_norm,1e-300)),
                interaction_fit_cubic_only_max_abs_resid=float(np.max(np.linalg.norm(fit_i3['residual'][:,j],axis=1))),
                interaction_fit_cubic_quintic_max_abs_resid=float(np.max(np.linalg.norm(fit_i35['residual'][:,j],axis=1))),
                interaction_fit_cubic_quintic_max_LOO_error=float(np.max(np.linalg.norm(fit_i35['loo_residual'][:,j],axis=1))),
                descriptor_fit_linear_cubic_quintic_max_LOO_error=float(np.max(np.linalg.norm(fit_d['loo_residual'][:,j],axis=1))),
                signed_C1=float(fit_d['coefficients'][0,j]@axis),
                signed_C3=float(fit_d['coefficients'][1,j]@axis),
                signed_C5=float(fit_d['coefficients'][2,j]@axis)))
        local=pd.DataFrame([r for r in fits if r['tolerance']==tol])
        postfit=local[local.clean_after_B]
        summary[tol]={
            'interaction_effective_p_1e8_to_3e8_postB_median':safe_median(postfit.power_1e8_to_3e8),
            'interaction_cubic_normalized_vector_norm_rel_range_postB_median':safe_median(postfit.scaled_cubic_coefficient_relative_range),
            'quintic_to_cubic_at_largest_epsilon_postB_median':safe_median(postfit.quintic_to_cubic_at_3e8),
            'centered_delta_odd_fit_design_condition':fit_d['design_condition'],
            'interaction_odd_fit_design_condition':fit_i35['design_condition'],
            'order2_odd_fit_design_condition':orderfit['design_condition'],
            'postB_all_interaction_signed_positive':bool(np.all(i[:,post]>0)),
            'postB_all_residual_signed_negative':bool(np.all(s[:,post]<0)),
            'postB_signed_linear_residual_normalized_median':safe_median(s[:,post]/(AMPS[:,None]/EPS_UNIT)),
            'postB_max_abs_interaction_over_delta_norm':float(np.max(np.linalg.norm(inter[:,post,:],axis=2)/np.maximum(np.linalg.norm(delta[:,post,:],axis=2),1e-300))),
            'preB_material_zero_max':float(np.max(np.linalg.norm(massB[:,sB<s0,:],axis=2))),
        }
    samp=pd.DataFrame(samp);fitdf=pd.DataFrame(fits)
    tolerance_rows=[]
    for eps in AMPS:
        key_s=f'standard__eps{eps:g}'
        key_t=f'tight__eps{eps:g}'
        for j,x in enumerate(xi):
            d0,d1=exp[key_s]['delta'][j],exp[key_t]['delta'][j]
            i0,i1=exp[key_s]['inter'][j],exp[key_t]['inter'][j]
            e0,e1=exp[key_s]['order2'][j],exp[key_t]['order2'][j]
            denom_i=max(float(np.linalg.norm(i0)),float(np.linalg.norm(i1)),1e-300)
            denom_d=max(float(np.linalg.norm(d0)),float(np.linalg.norm(d1)),1e-300)
            tolerance_rows.append(dict(eps_rel=float(eps),xi=float(x),clean_after_B=bool(post[j]),
               interaction_abs_discrepancy=float(np.linalg.norm(i0-i1)),
               interaction_relative_discrepancy=float(np.linalg.norm(i0-i1)/denom_i),
               delta_abs_discrepancy=float(np.linalg.norm(d0-d1)),
               delta_relative_discrepancy=float(np.linalg.norm(d0-d1)/denom_d),
               order2_abs_discrepancy=float(np.linalg.norm(e0-e1)) if np.isfinite(e0).all() and np.isfinite(e1).all() else np.nan))
    tol_df=pd.DataFrame(tolerance_rows)
    for eps in AMPS:
        seg=samp[samp.eps_rel==eps]
        aa=seg[seg.tolerance=='standard'].sort_values('xi')
        bb=seg[seg.tolerance=='tight'].sort_values('xi')
        if not np.array_equal(aa.xi.to_numpy(),bb.xi.to_numpy()):raise RuntimeError('tolerance xi grid mismatch')
    if max_close>1e-10:raise RuntimeError(f'descriptor factor closure failed {max_close}')
    if min_preB_material>1e-12:raise RuntimeError(f'causality B failed: {min_preB_material}')
    summary['max_factor_identity_closure']=max_close
    summary['max_before_B_material_sensitivity_norm']=min_preB_material
    summary['all_amplitudes_validated']=len(samp)==2*5*47
    # Keep threshold and physical interpretation visible; do NOT automatically
    # certify a cubic asymptote from a polynomial fit or exact algebraic closure.
    summary['conclusion_gate']='DIAGNOSTIC_ONLY_PENDING_ASYMPTOTIC_AND_CONSTRAINT_AND_REVERSE_PATH'
    tolerance_summaries=[]
    for eps in AMPS:
        sub=tol_df[(tol_df.eps_rel==eps)&(tol_df.clean_after_B)]
        tolerance_summaries.append(dict(eps_rel=float(eps),
             median_postB_interaction_relative_solver_difference=safe_median(sub.interaction_relative_discrepancy),
             maximum_postB_interaction_relative_solver_difference=float(sub.interaction_relative_discrepancy.max()),
             median_postB_delta_relative_solver_difference=safe_median(sub.delta_relative_discrepancy)))
    summary['cross_tolerance_postB']=tolerance_summaries
    # A tolerance disagreement indicates limited numerical resolution, not a
    # physical origin for the interaction. Equal tolerances alone cannot
    # exclude correlated systematic errors.
    summary['smallest_amplitude_numerical_floor_flag']=bool(
        tolerance_summaries[0]['median_postB_interaction_relative_solver_difference']>0.05)
    samp.to_csv(args.out/'five_amplitude_samples.csv',index=False,float_format='%.16g')
    tol_df.to_csv(args.out/'tolerance_discrepancy.csv',index=False,float_format='%.16g')
    fitdf.to_csv(args.out/'samplewise_odd_fits.csv',index=False,float_format='%.16g')
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (args.out/'REPORT.md').write_text('# Q550 v18.4 — five-amplitude scaling results\n\n'
         f'Fixed A axis: `v18.3 standard eps=1e-8`; B crossing xi={cross:.15f}.\n\n'
         'Amplitudes: '+', '.join(f'{e:g}' for e in AMPS)+'; tolerances standard and tight.\n\n'
         f'Exact 11D factor-decomposition max closure: `{max_close:.5g}`.\n\n'
         'No constraint certificate, no forward/reverse memory loop, and no material-flip certificate.\n\n'
         'The cubic term in the geometric descriptor is algebraically compatible with symmetric +/- secants.\n'
         'Asymptotic claims must account for solver noise and five-point odd-order fit stability.\n\n'
         f'Standard: median interaction effective exponent over 19 clean post-B samples = '
         f'`{summary["standard"]["interaction_effective_p_1e8_to_3e8_postB_median"]}`.\n'
         f'Tight: median interaction effective exponent = '
         f'`{summary["tight"]["interaction_effective_p_1e8_to_3e8_postB_median"]}`.\n\n'
         'See `five_amplitude_samples.csv`, `samplewise_odd_fits.csv`, and diagnostic figures.\n')
    for name,field,ylabel,logy in [
         ('interaction_scaled','interaction_over_eps_cube_norm','||I(eps)||(1e-8/eps)^3',False),
         ('interaction_raw','interaction_norm','||I_ABC(eps)||',True),
         ('signed_cancellation','signed_delta','signed Delta along frozen A',False)]:
        fig,axes=plt.subplots(1,2,figsize=(12,4.2),sharey=False)
        for a,tol in zip(axes,('standard','tight')):
            sub=samp[(samp.tolerance==tol)&(samp.clean_after_B)]
            y=[float(np.median(sub[sub.eps_rel==eps][field])) for eps in AMPS]
            a.plot(AMPS/EPS_UNIT,y,'o-',label=field)
            if name=='interaction_raw':
                base=y[1];a.plot(AMPS/EPS_UNIT,[base*(eps/AMPS[1])**3 for eps in AMPS],'--',label='cubic through eps=1e-8')
                a.set_yscale('log')
            if name=='signed_cancellation':
                alt=[float(np.median(sub[sub.eps_rel==eps]['signed_interaction'])) for eps in AMPS]
                other=[float(np.median(sub[sub.eps_rel==eps]['signed_single_factor_sum'])) for eps in AMPS]
                a.plot(AMPS/EPS_UNIT,alt,'s-',label='signed interaction')
                a.plot(AMPS/EPS_UNIT,other,'^-',label='single-factor secants')
                a.axhline(0,color='gray',linewidth=.6)
            a.set(title=tol,xlabel='epsilon / 1e-8',ylabel=ylabel)
            a.grid(alpha=.2);a.legend(fontsize=8)
        fig.tight_layout();fig.savefig(args.out/f'{name}.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.5,4.5))
    grouped=[tol_df[(tol_df.eps_rel==eps)&(tol_df.clean_after_B)] for eps in AMPS]
    y=[float(np.median(g.interaction_relative_discrepancy)) for g in grouped]
    y2=[float(np.median(g.delta_relative_discrepancy)) for g in grouped]
    ax.plot(AMPS/EPS_UNIT,y,'o-',label='interaction standard vs tight')
    ax.plot(AMPS/EPS_UNIT,y2,'s-',label='residual standard vs tight')
    ax.axhline(0.05,color='gray',linewidth=.7,linestyle='--',label='5% resolution warning')
    ax.set(xlabel='epsilon / 1e-8',ylabel='median relative disagreement (19 clean post-B samples)',
           title='Tolerance sensitivity across all five amplitudes')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout()
    fig.savefig(args.out/'tolerance_sensitivity.png',dpi=170);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoints',type=Path,default=ROOT/'checkpoint')
    p.add_argument('--out',type=Path,default=ROOT/'results')
    run(p.parse_args())
