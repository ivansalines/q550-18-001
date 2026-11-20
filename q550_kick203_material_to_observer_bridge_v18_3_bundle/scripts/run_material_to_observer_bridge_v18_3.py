#!/usr/bin/env python3
"""Q550 v18.3: reproducible MATERIAL->(BA,BB,M)->dlocal observer bridge.

Offline analysis of the actual archived v17.9, v18.1 and v18.2 results.
No new orbit or artificial full Jacobian is generated. The descriptor 11D
and state 11D are NOT identified as the same physical observable.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
from history_bridge_v16 import descriptor
from frame_geometry import shapley_three_factors
from bridge_math import factor_secants, match_grid, unit, angle_deg

VERSION='2026-10-16-Q550-kick203-material-to-observer-full11D-bridge-v18.3'


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
    return h.hexdigest()


def exactly_one(z, suffix):
    names=[name for name in z.namelist() if name==suffix or name.endswith('/'+suffix)]
    if len(names)!=1:raise ValueError(f'expected exactly one {suffix}, got {names}')
    return z.read(names[0])


def npz_bytes(blob):
    with np.load(io.BytesIO(blob),allow_pickle=False) as f:
        return {key:np.asarray(f[key]) for key in f.files}


def load_sources(path):
    with zipfile.ZipFile(path) as z2:
        material=npz_bytes(exactly_one(z2,'full11D_material_arrays.npz'))
        summary=json.loads(exactly_one(z2,'results/kick203_v18_2/summary.json'))
        b1=exactly_one(z2,'q550_kick203_second_order_observer_v18_1_sendback.zip')
    with zipfile.ZipFile(io.BytesIO(b1)) as z1:
        observer=npz_bytes(exactly_one(z1,'second_order_full11D_arrays.npz'))
        b179=exactly_one(z1,'q550_kick203_event_micro_v17_9_sendback.zip')
    with zipfile.ZipFile(io.BytesIO(b179)) as z179:
        frame=npz_bytes(exactly_one(z179,'full_frames_cocycle_arrays.npz'))
        shap=npz_bytes(exactly_one(z179,'oriented_shapley_arrays.npz'))
    return material,summary,observer,frame,shap,{'v18_1_zip_sha256':hashlib.sha256(b1).hexdigest(),
                                             'v17_9_zip_sha256':hashlib.sha256(b179).hexdigest()}


def verify_archived_observer(xi, selected, sh_old, sh_new, tol, eps, amp, axis):
    prefix=f'{tol}__kick203__eps{eps:g}'
    centered=sh_old[prefix+'__hquarter__dlocal_delta'][selected]
    centered_alloc=sh_old[prefix+'__hquarter__shapley_ABC'][selected]
    primary=f'{prefix}__A_hquarter'
    second=sh_new[primary+'__B_hsixteenth__non_straddling2__dlocal_delta']
    second_alloc=sh_new[primary+'__B_hsixteenth__non_straddling2__shapley_ABC']
    return centered,centered_alloc,second,second_alloc


def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    bundle_hash=sha256(args.sendback)
    if args.strict_manifest:
        provenance=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
        if provenance['v18_2_sendback_sha256'] != bundle_hash:
            raise RuntimeError('v18.2 input sendback hash mismatch')
    material,old_summary,observer,frames,hshap,nested_hash=load_sources(args.sendback)
    xi=np.asarray(material['xi'],float); n=len(xi)
    if not np.array_equal(xi,np.asarray(material['v18_1_xi'],float)) or not np.array_equal(xi,observer['xi']):
        raise RuntimeError('v18.1 material and observer grids fail exact alignment')
    full_xi=frames['standard__kick203__eps0__sign+0__xi']
    idx=match_grid(full_xi,xi)
    sA=np.asarray(material['sA'],float); sB=np.asarray(material['sB'],float)
    s0=float(material['kick_support_s'][0]); crossing=float(material['xi_B_crossing'][0])
    if not np.isclose(observer['crossing_xi'][0],crossing,rtol=0,atol=1e-12):
        raise RuntimeError('crossing mismatch')
    hhalf=float(observer['B_support_heights'][0]); hquarter=float(observer['B_support_heights'][1]); h16=float(observer['B_support_heights'][3])
    # Use physical-support distances, not xi distances, for causal/stencil gates.
    centered_clean=np.abs(sB-s0)>hquarter + 5e-14
    event=np.abs(sB-s0)<=5e-12
    left_clean=centered_clean&(sB<s0)
    right_clean=centered_clean&(sB>s0)
    near=~centered_clean
    if not (np.any(left_clean) and np.any(right_clean) and np.any(near)):
        raise RuntimeError('insufficient clean left/right and contaminated event samples')
    rows=[]; signrows=[]; arrays={'xi':xi,'sA':sA,'sB':sB,
                                 's0':np.array([s0]),'crossing_xi':np.array([crossing]),
                                 'centered_stencil_clean':centered_clean,
                                 'centered_stencil_straddling':near}
    max_errors={'source_frame_descriptor':0.,'archived_shapley':0.,'factor_split':0.,
                'observer_shapley':0.,'observer_identity':0.,'material_causality_before_B':0.}
    experiments=[]
    for tol in ('standard','tight'):
        basezero=f'{tol}__kick203__eps0__sign+0'
        for eps in (1e-8,3e-8):
            history=f'{tol}__kick203__eps{eps:g}'
            key=f'{tol}__eps{eps:g}'
            primary=f'{history}__A_hquarter'
            axis=unit(observer[primary+'__frozen_A_11D'])
            a_plus=history+'__sign+1';a_minus=history+'__sign-1'
            cached={}; source_delta=[]; secants_all=[]; interaction_all=[]; sh_exact=[]
            for j,k in enumerate(idx):
                F0=tuple(frames[basezero+'__'+suffix][k] for suffix in ('hquarter__BA','hquarter__BB','M'))
                Fm=tuple(frames[a_minus+'__'+suffix][k] for suffix in ('hquarter__BA','hquarter__BB','M'))
                Fp=tuple(frames[a_plus+'__'+suffix][k] for suffix in ('hquarter__BA','hquarter__BB','M'))
                origin=descriptor(*F0)
                plus=descriptor(*Fp); minus=descriptor(*Fm)
                max_errors['source_frame_descriptor']=max(max_errors['source_frame_descriptor'],
                    np.linalg.norm(plus-frames[a_plus+'__hquarter__dlocal'][k]),
                    np.linalg.norm(minus-frames[a_minus+'__hquarter__dlocal'][k]),
                    np.linalg.norm(origin-frames[basezero+'__hquarter__dlocal'][k]))
                secs,interaction,delta=factor_secants(Fm,Fp,descriptor)
                allocated,alldelta,err=shapley_three_factors(Fm,Fp)
                max_errors['factor_split']=max(max_errors['factor_split'],
                    np.linalg.norm(delta-np.sum(secs,axis=0)-interaction))
                max_errors['archived_shapley']=max(max_errors['archived_shapley'],err,
                     np.linalg.norm(delta-alldelta))
                secants_all.append(secs);interaction_all.append(interaction)
                source_delta.append(delta);sh_exact.append(allocated)
                # Source factors and material states are from the SAME +- histories,
                # but section A and section B represent DIFFERENT support times.
                if args.save_factors:
                    for nm,vec in [('dBA',Fp[0]-Fm[0]),('dBB',Fp[1]-Fm[1]),('dM',Fp[2]-Fm[2])]:
                        cached.setdefault(nm,[]).append(vec)
            source_delta=np.asarray(source_delta)
            secants_all=np.asarray(secants_all)
            interaction_all=np.asarray(interaction_all)
            sh_exact=np.asarray(sh_exact)
            centered,archived_alloc,second,second_alloc=verify_archived_observer(xi,idx,hshap,observer,tol,eps,0.,axis)
            max_errors['archived_shapley']=max(max_errors['archived_shapley'],
                     np.nanmax(np.linalg.norm(source_delta-centered,axis=1)),
                     np.nanmax(np.linalg.norm(sh_exact-archived_alloc,axis=(1,2))))
            valid=np.isfinite(second).all(axis=1)&(~event)
            if not np.all(valid[~event]):raise RuntimeError('invalid one-sided second-order sample outside event')
            max_errors['observer_shapley']=max(max_errors['observer_shapley'],
                     np.nanmax(np.linalg.norm(second-np.sum(second_alloc,axis=1),axis=1)))
            correction=second-centered
            ua=material[key+'__uA'];ub=material[key+'__uB']
            la=material[key+'__linear_A'];lb=material[key+'__linear_B']
            B_norm=np.linalg.norm(ub,axis=1)
            max_errors['material_causality_before_B']=max(max_errors['material_causality_before_B'],
                float(np.max(B_norm[sB<s0])))
            if ua.shape!=(n,11) or ub.shape!=(n,11):raise RuntimeError('unexpected material response shapes')
            subtag=f'{tol}__eps{eps:g}'
            # The material amplitude is dimensionally distinct from dlocal.
            matches=[m for m in old_summary['experiments'] if m['tolerance']==tol and np.isclose(m['eps_rel'],eps,rtol=1e-12)]
            if len(matches)!=1:raise RuntimeError('unmatched v18.2 source amplitude')
            amp=float(matches[0]['amp_absolute'])
            corrected_predicted=source_delta+correction
            max_errors['observer_identity']=max(max_errors['observer_identity'],
                    np.nanmax(np.linalg.norm(corrected_predicted-second,axis=1)))
            scalar=np.einsum('ij,j->i',second,axis)
            signed_centered=source_delta@axis
            signed_A=second_alloc[:,0,:]@axis
            signed_B=second_alloc[:,1,:]@axis
            signed_M=second_alloc[:,2,:]@axis
            arrays.update({
                key+'__material_uA_11D':ua, key+'__material_uB_11D':ub,
                key+'__centered_delta_11D':source_delta,
                key+'__order2_delta_11D':second,
                key+'__B_stencil_correction_11D':correction,
                key+'__centered_factor_secants_3x11':secants_all,
                key+'__centered_interaction_11D':interaction_all,
                key+'__centered_shapley_3x11':sh_exact,
                key+'__order2_shapley_3x11':second_alloc,
                key+'__frozen_axis_11D':axis,
                key+'__material_bodyA_rank1_11D':material[key+'__body_A_rank1_11D'],
                key+'__material_bodyA_rank2_11D':material[key+'__body_A_rank2_11D'],
            })
            if args.save_factors:
                for nm,vecs in cached.items():arrays[key+'__'+nm]=np.asarray(vecs)
            prefix=f'{history}__A_hquarter'
            # Pure observer sensitivity, when same ODE histories and frozen A/M are used.
            d1=observer[prefix+'__B_hsixteenth__non_straddling1__dlocal_delta']
            for j in range(n):
                def nrm(x):return float(np.linalg.norm(x))
                rows.append(dict(tolerance=tol,eps_rel=eps,amp_absolute=amp,xi=float(xi[j]),
                   xi_minus_B=float(xi[j]-crossing),sA=float(sA[j]),sB=float(sB[j]),
                   before_B_restart=bool(sB[j]<s0),centered_stencil_clean=bool(centered_clean[j]),
                   centered_stencil_straddling=bool(near[j]),event_sample_censored=bool(event[j]),
                   material_A_norm=nrm(ua[j]),material_B_norm=nrm(ub[j]),
                   material_A_predictor_error=nrm(ua[j]-la[j]),
                   material_B_predictor_error=nrm(ub[j]-lb[j]),
                   material_A_signed_rank1=float(material[key+'__body_A_rank1_11D'][j]@material[key+'__direction_initial']),
                   material_A_signed_rank2=float(material[key+'__body_A_rank2_11D'][j]@material[key+'__direction_initial']),
                   centered_delta_norm=nrm(source_delta[j]),
                   centered_signed=float(signed_centered[j]),
                   order2_delta_norm=nrm(second[j]) if valid[j] else np.nan,
                   order2_signed=float(scalar[j]) if valid[j] else np.nan,
                   stencil_correction_norm=nrm(correction[j]) if valid[j] else np.nan,
                   stencil_correction_signed=float(correction[j]@axis) if valid[j] else np.nan,
                   observer_order1_vs_2_norm=nrm(d1[j]-second[j]) if valid[j] else np.nan,
                   shapley_A_signed=float(signed_A[j]) if valid[j] else np.nan,
                   shapley_B_signed=float(signed_B[j]) if valid[j] else np.nan,
                   shapley_M_signed=float(signed_M[j]) if valid[j] else np.nan,
                   shapley_A_norm=nrm(second_alloc[j,0]) if valid[j] else np.nan,
                   shapley_B_norm=nrm(second_alloc[j,1]) if valid[j] else np.nan,
                   shapley_M_norm=nrm(second_alloc[j,2]) if valid[j] else np.nan,
                   centered_secant_A_signed=float(secants_all[j,0]@axis),
                   centered_secant_B_signed=float(secants_all[j,1]@axis),
                   centered_secant_M_signed=float(secants_all[j,2]@axis),
                   centered_interaction_signed=float(interaction_all[j]@axis),
                   centered_interaction_norm=nrm(interaction_all[j]),
                   centered_interaction_fraction=float(nrm(interaction_all[j]) / max(nrm(source_delta[j]),1e-30)),
                   normalized_material_A_norm=nrm(ua[j]),
                   normalized_descriptor_response_norm=nrm(second[j])/(2*amp) if valid[j] else np.nan))
            post_sign=scalar[right_clean];pre_sign=scalar[left_clean]
            if not np.all(np.isfinite(post_sign)) or not np.all(np.isfinite(pre_sign)):
                raise RuntimeError('missing clean signed observations')
            signrows.append(dict(tolerance=tol,eps_rel=eps,n_pre=int(left_clean.sum()),
                n_post=int(right_clean.sum()),pre_positive=int(np.sum(pre_sign>0)),
                post_negative=int(np.sum(post_sign<0)),median_pre_signed=float(np.median(pre_sign)),
                median_post_signed=float(np.median(post_sign)),
                median_post_stencil_correction_over_delta=float(np.median(
                    np.linalg.norm(correction[right_clean],axis=1)/
                    np.maximum(np.linalg.norm(second[right_clean],axis=1),1e-30))),
                median_post_interaction_over_centered=float(np.median(
                    np.linalg.norm(interaction_all[right_clean],axis=1)/
                    np.maximum(np.linalg.norm(source_delta[right_clean],axis=1),1e-30)))))
            experiments.append(dict(tolerance=tol,eps_rel=eps,absolute_kick_amplitude=amp,
                    clean_before_n=int(left_clean.sum()),clean_after_n=int(right_clean.sum()),
                    sign_before_positive=bool(np.all(pre_sign>0)),sign_after_negative=bool(np.all(post_sign<0)),
                    post_signed_median=float(np.median(post_sign)),
                    material_A_at_B_crossing_minus_plus_angle_deg=float(angle_deg(
                        ua[np.where(xi<crossing)[0][-1]],ua[np.where(xi>crossing)[0][0]]))))
    for key,err in max_errors.items():
        if err>args.closure_gate:raise RuntimeError(f'{key}: closure or causality gate failed: {err} > {args.closure_gate}')
    df=pd.DataFrame(rows);df.to_csv(args.out/'bridge_samples.csv',index=False,float_format='%.16g')
    sdf=pd.DataFrame(signrows);sdf.to_csv(args.out/'robust_sign_audit.csv',index=False,float_format='%.16g')
    np.savez_compressed(args.out/'bridge_full11D_arrays.npz',**arrays)
    # Figures strictly from archived samples; never connect across censored event.
    primary=df[(df.tolerance=='standard')&(df.eps_rel==1e-8)].copy()
    fig,axes=plt.subplots(1,2,figsize=(13,4.7))
    for ax in axes:
        for col,label in [('centered_signed','archived centered v17.9'),
                          ('order2_signed','one-sided second order v18.1'),
                          ('stencil_correction_signed','observer correction')]:
            ax.plot(primary.xi_minus_B,primary[col],marker='.',markersize=3,linewidth=1.2,label=label)
        radius=hquarter/((sB[-1]-sB[0])/(xi[-1]-xi[0]))
        ax.axvspan(-radius,radius,color='grey',alpha=.13)
        ax.axvline(0,color='black',linestyle=':',linewidth=.8)
        ax.grid(alpha=.22)
        ax.set_xlabel('xi - xi_B')
    axes[0].set(xlim=(-.041,.041),ylabel='signed dlocal delta on frozen A',
                title='Full local transition')
    axes[1].set(xlim=(0,.041),ylim=(-.00029,.00012),title='Post-kick residual detail')
    axes[1].legend(fontsize=7,loc='lower right')
    fig.suptitle('v18.3 — source-history response and B observer correction')
    fig.tight_layout();fig.savefig(args.out/'signed_history_vs_observer.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5))
    for col,label in [('shapley_A_signed','Shapley A'),('shapley_B_signed','Shapley B'),('shapley_M_signed','Shapley M'),
                      ('order2_signed','total one-sided order2')]:
        ax.plot(primary.xi_minus_B,primary[col],marker='.',markersize=3,linewidth=1.1,label=label)
    ax.axvline(0,color='black',linestyle=':',linewidth=.8)
    ax.set(xlabel='xi - xi_B',ylabel='signed full11D descriptor contribution',title='v18.3 — exact 3-factor Shapley observer attribution')
    ax.set_xlim(-.041,.041)
    ax.grid(alpha=.22);ax.legend(fontsize=8);fig.tight_layout()
    fig.savefig(args.out/'oriented_A_B_M_shapley.png',dpi=170);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(10,8),sharex=True)
    for col,label in [('material_A_norm','state sensitivity ||u_A||'),('material_B_norm','state sensitivity ||u_B||')]:
        axs[0].plot(primary.xi_minus_B,primary[col],marker='.',markersize=2,label=label)
    axs[0].set_ylabel('state sensitivity per unit kick amplitude');axs[0].legend();axs[0].grid(alpha=.22)
    for col,label in [('centered_delta_norm','centered descriptor ||Delta||'),('order2_delta_norm','order2 descriptor ||Delta||'),
                      ('stencil_correction_norm','observer correction ||Delta||')]:
        axs[1].plot(primary.xi_minus_B,primary[col],marker='.',markersize=2,label=label)
    axs[1].set(xlabel='xi - xi_B',ylabel='raw descriptor norm',title='Different observables — do NOT identify vertical scales')
    axs[1].legend();axs[1].grid(alpha=.22)
    for ax in axs:ax.axvline(0,color='black',linestyle=':',linewidth=.8);ax.set_xlim(-.041,.041)
    fig.tight_layout();fig.savefig(args.out/'material_causality_vs_observer.png',dpi=170);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,4.4))
    for col,label in [('centered_secant_A_signed','conditional A'),('centered_secant_B_signed','conditional B'),
                      ('centered_secant_M_signed','conditional M'),('centered_interaction_signed','interaction')]:
        axs[0].plot(primary.xi_minus_B,primary[col],linewidth=1.1,label=label)
    axs[0].legend(fontsize=7);axs[0].set(ylabel='signed secant on frozen A',xlabel='xi - xi_B',title='Centered factor secants + exact interaction')
    axs[1].plot(primary.xi_minus_B,primary.centered_interaction_fraction,marker='.',markersize=2)
    axs[1].set(xlabel='xi - xi_B',ylabel='||interaction|| / ||centered delta||',title='Nonadditivity of separate factor responses')
    for ax in axs:ax.grid(alpha=.22);ax.axvline(0,color='black',linestyle=':',linewidth=.8);ax.set_xlim(-.041,.041)
    fig.tight_layout();fig.savefig(args.out/'conditional_factor_secants.png',dpi=170);plt.close(fig)
    summary=dict(version=VERSION,outer_source_sha256=bundle_hash,nested_source_sha256=nested_hash,
        n_xi=n,xi_B_crossing=crossing,s0=s0,hquarter_B_support=hquarter,
        h16_B_support=h16,centered_event_stencil_halfwidth_xi=float(hquarter/((sB[-1]-sB[0])/(xi[-1]-xi[0]))),
        n_clean_before=int(left_clean.sum()),n_clean_after=int(right_clean.sum()),
        n_centered_stencil_straddling=int(near.sum()),event_exact_point_censored=int(event.sum()),
        max_identity_errors=max_errors,experiments=experiments,
        physical_constraint_status=old_summary['status'],
        interpretation=('Exact secants split the archival centered physical-history dlocal response into '
                        'factor-conditional A/B/M terms plus a measured nonlinear interaction. The difference '
                        'between archived centered dlocal response and v18.1 B one-sided second-order response '
                        'is an exact observer change for the SAME +/- histories. Material 11D states are '
                        'tracked as causal controls, not equated to geometric descriptor 11D. The B observer '
                        'response has an impulse-entry discontinuity, which is NOT a material torsional flip. '
                        'No certified conserved physical invariant or globally identified full 11x11 pushforward.'))
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# v18.3 — Material-to-Observer full11D Bridge','','## Source alignment',
      f'- {n} identical xi samples; physical B crossing xi={crossing:.15f}.',
      f'- B support kick s0={s0:.15f}; hquarter stencil half-width={summary["centered_event_stencil_halfwidth_xi"]:.9g} xi.',
      f'- Clean left/right: {int(left_clean.sum())}/{int(right_clean.sum())}; centered straddling samples: {int(near.sum())}.','',
      '## Exact finite-history identity','',
      '`Delta_dlocal_cent = S_A + S_B + S_M + Interaction` (conditional midpoint factor secants).',
      '`Delta_dlocal_order2 = Delta_dlocal_cent + ObserverStencilCorrection`.','',
      'These are exact vector identities, not causal identification of an intrinsic torsion.','',
      '## Signed residual away from centered-stencil contamination','']
    for v in experiments:
      lines.append(f'- {v["tolerance"]}, eps={v["eps_rel"]:g}: pre positive={v["sign_before_positive"]}; '
         f'post negative={v["sign_after_negative"]}; post median={v["post_signed_median"]:.8e}.')
    lines+=['','## Controls','']
    for k,v in max_errors.items():lines.append(f'- {k}: {v:.4e}')
    lines+=['','## Cautions','',
         'The directional `+/-` material impulse is an on-flow perturbation; its invariant manifold is not independently certified.',
         'The section-A and section-B states are sampled at different physical support coordinates.',
         'The descriptor and state sensitivities are not interchanged, subtracted, or implicitly dot-producted.',
         'The factor midpoint construction does not identify the complete Jacobian from material states to dlocal.',
         'Centered stencils near B crossing mix both sides of the restart; the one-sided sample exactly at the crossing is censored.',
         'An observable residual sign flip is not automatically a material torsional flip.','']
    (args.out/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'status':'OK', 'xi_samples':n,'clean_before':int(left_clean.sum()),
        'clean_after':int(right_clean.sum()),'straddling':int(near.sum()),'errors':max_errors,
        'signs':signrows},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sendback',type=Path,default=ROOT/'inputs/q550_kick203_on_shell_full11d_material_transport_v18_2_sendback.zip')
    parser.add_argument('--out',type=Path,default=ROOT/'results/kick203_v18_3')
    parser.add_argument('--closure-gate',type=float,default=1e-9)
    parser.add_argument('--no-strict-manifest',dest='strict_manifest',action='store_false',default=True)
    parser.add_argument('--save-factors',action='store_true',default=True)
    run(parser.parse_args())
