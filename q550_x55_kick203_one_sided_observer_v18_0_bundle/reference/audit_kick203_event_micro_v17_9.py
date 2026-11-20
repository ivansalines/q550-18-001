#!/usr/bin/env python3
"""Event-aligned signed 11D Shapley one-sided lens; no physical flip claims.

The frozen 11D reference axis is the A direction at the first point at least
0.025 xi after B crossing. Diagnostics use the original source Shapley arrays.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def unit(x):
    n=float(np.linalg.norm(x))
    return np.asarray(x,float)/n if n>1e-30 else np.full_like(x,np.nan)


def analyze(args):
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(args.allocation_csv)
    src=np.load(args.shapley_npz)
    refs=np.load(args.reference_shapley_npz) if args.reference_shapley_npz else None
    old_df=pd.read_csv(args.reference_allocation_csv) if args.reference_allocation_csv else None
    xb=float((df.xi-df.B_relative).median())
    records=[]; regressions=[]; per_variant=[]; probes=[]
    for key in sorted(k for k in src.files if k.endswith('__shapley_ABC')):
        tag=key.removesuffix('__shapley_ABC')
        parts=tag.split('__')
        tol=parts[0];kick=float(parts[1].removeprefix('kick'))
        er=float(parts[2].removeprefix('eps'))
        stencil=parts[3]
        part=df[(df.tolerance==tol)&np.isclose(df.kick_xi,kick)&np.isclose(df.eps_rel,er)&(df.stencil==stencil)].sort_values('xi').reset_index(drop=True)
        X=part.xi.to_numpy(float)
        if len(part)!=len(src[key]):raise RuntimeError(f'{tag}: inconsistent arrays/csv {len(part)} vs {len(src[key])}')
        abc=np.asarray(src[key],float)
        delta=np.asarray(src[tag+'__dlocal_delta'],float)
        err=np.linalg.norm(np.sum(abc,axis=1)-delta,axis=1)
        if not (np.nanmax(err)<1e-9):raise RuntimeError('algebraic closure broken '+tag)
        kref=int(np.argmin(abs(X-(xb+.025))))
        e1=unit(abc[kref,0]);S=abc[:,0]+abc[:,1]
        e2=unit(S[kref]-np.dot(S[kref],e1)*e1)
        if not np.isfinite(e2).all():raise RuntimeError('Cannot define frozen transverse axis')
        for i,x in enumerate(X):
            A,B,M=abc[i]
            a,b,m,d,s=[float(np.linalg.norm(v)) for v in (A,B,M,delta[i],S[i])]
            cosAB=float(np.dot(A,B)/a/b) if a*b>0 else np.nan
            cosSM=float(np.dot(S[i],M)/s/m) if s*m>0 else np.nan
            eA=unit(A)
            row=dict(variant=tag,xi=float(x),distance_from_crossing=float(x-xb),side='left' if x<xb-1e-11 else 'right' if x>xb+1e-11 else 'at_cross',
              A_norm=a,B_norm=b,M_norm=m,S_norm=s,delta_norm=d,
              AB_angle_deg=float(np.degrees(np.arccos(np.clip(cosAB,-1,1)))) if np.isfinite(cosAB) else np.nan,
              SM_angle_deg=float(np.degrees(np.arccos(np.clip(cosSM,-1,1)))) if np.isfinite(cosSM) else np.nan,
              A_minus_B_norm=a-b, B_over_A=b/a if a else np.nan,
              S_over_AplusB=s/(a+b) if (a+b)>0 else np.nan,
              S_local_parallel=float(np.dot(S[i],eA)) if np.isfinite(eA).all() else np.nan,
              M_local_parallel=float(np.dot(M,eA)) if np.isfinite(eA).all() else np.nan,
              delta_local_parallel=float(np.dot(delta[i],eA)) if np.isfinite(eA).all() else np.nan,
              S_frozen_parallel=float(np.dot(S[i],e1)),
              M_frozen_parallel=float(np.dot(M,e1)),
              delta_frozen_parallel=float(np.dot(delta[i],e1)),
              S_frozen_transverse=float(np.dot(S[i],e2)),
              M_frozen_transverse=float(np.dot(M,e2)),
              delta_frozen_transverse=float(np.dot(delta[i],e2)),
              S_transverse_11D=float(np.linalg.norm(S[i]-np.dot(S[i],eA)*eA)) if np.isfinite(eA).all() else np.nan,
              closure_error=float(err[i]),
              frozen_axis_xi=float(X[kref]))
            records.append(row)
        # Need regression on exact 0.025 old samples; never interpolate.
        if refs is not None and key in refs.files:
            refabc=np.asarray(refs[key],float)
            refpart=old_df[(old_df.tolerance==tol)&np.isclose(old_df.kick_xi,kick)&np.isclose(old_df.eps_rel,er)&(old_df.stencil==stencil)].sort_values('xi').reset_index(drop=True)
            oldx=refpart.xi.to_numpy(float)
            diffs=[];norms=[];nmatched=0
            for j,x in enumerate(oldx):
                i=int(np.argmin(abs(X-x)))
                if abs(X[i]-x)<1e-8:
                    diffs.append(float(np.linalg.norm(abc[i]-refabc[j])))
                    norms.append(float(np.linalg.norm(refabc[j])))
                    nmatched+=1
            regressions.append(dict(variant=tag,matched=nmatched,
                  abs_RMS=float(np.sqrt(np.mean(np.square(diffs)))) if diffs else np.nan,
                  relative_RMS=float(np.sqrt(np.mean(np.square(diffs)))/max(np.sqrt(np.mean(np.square(norms))),1e-300)) if diffs else np.nan,
                  max_abs=float(np.max(diffs)) if diffs else np.nan))
        xi_b=xb
        below=np.where(X<xi_b)[0];above=np.where(X>xi_b)[0]
        if len(below)==0 or len(above)==0:raise RuntimeError('no two sides of event')
        left=int(below[-1]);right=int(above[0]);
        minleft=float(xi_b-X[left]);minright=float(X[right]-xi_b)
        per_variant.append(dict(variant=tag,crossing=xi_b,nearest_left=float(X[left]),nearest_right=float(X[right]),
           left_gap=minleft,right_gap=minright,
           left_B_norm=float(np.linalg.norm(abc[left,1])),right_B_norm=float(np.linalg.norm(abc[right,1])),
           left_delta_norm=float(np.linalg.norm(delta[left])),right_delta_norm=float(np.linalg.norm(delta[right])),
           right_AB_angle_deg=float(np.degrees(np.arccos(np.clip(np.dot(abc[right,0],abc[right,1])/max(np.linalg.norm(abc[right,0])*np.linalg.norm(abc[right,1]),1e-300),-1,1))))))
        for offset in [.025,.01,.005,.001,.00025,.00005,.00001]:
            row=dict(variant=tag,offset_xi=offset)
            for sign,side in [(-1,'left'),(+1,'right')]:
                target=xi_b+sign*offset
                j=int(np.argmin(abs(X-target)))
                row[side+'_actual_offset']=float(X[j]-xi_b)
                row[side+'_B_norm']=float(np.linalg.norm(abc[j,1]))
                row[side+'_delta_norm']=float(np.linalg.norm(delta[j]))
                row[side+'_S_norm']=float(np.linalg.norm(S[j]))
            probes.append(row)
    res=pd.DataFrame(records);res.to_csv(out/'all_variants_signed_vectors.csv',index=False)
    pd.DataFrame(per_variant).to_csv(out/'event_two_sided_summary.csv',index=False)
    pd.DataFrame(probes).to_csv(out/'event_symmetric_offsets.csv',index=False)
    pd.DataFrame(regressions).to_csv(out/'reference_step0025_regression.csv',index=False)
    primary=res[res.variant=='standard__kick203__eps1e-08__hhalf'].copy()
    window=primary[(primary.xi>=xb-.15)&(primary.xi<=xb+.35)].copy()
    # Common two-panel view, separate baseline on each side; no line across event.
    fig,ax=plt.subplots(2,1,figsize=(9,7),sharex=True)
    for side,group in window.groupby('side',sort=False):
        if side=='at_cross':continue
        for col,label in [('A_norm','||A||'),('B_norm','||B||'),('S_norm','||A+B||'),('delta_norm','||delta||')]:
            ax[0].plot(group.distance_from_crossing,group[col],label=f'{label} [{side}]')
        for col,label in [('S_frozen_parallel','S e1'),('M_frozen_parallel','M e1'),('delta_frozen_parallel','delta e1'),
                          ('delta_frozen_transverse','delta e2')]:
            ax[1].plot(group.distance_from_crossing,group[col],label=f'{label} [{side}]')
    for a in ax:a.axvline(0,ls=':',color='black',lw=1);a.legend(fontsize=7,ncol=2)
    ax[0].set_ylabel('full11D norm')
    ax[1].set_ylabel('frozen signed 11D projection')
    ax[1].set_xlabel('xi minus B crossing')
    fig.suptitle('Q550 kick203 event-micro lens — Shapley diagnostics, not physical impulse')
    fig.tight_layout();fig.savefig(out/'event_micro_primary_two_sides.png',dpi=180);plt.close(fig)
    # Symmetric paired event-limit scan on logarithmic delta.
    pairs=pd.DataFrame(probes)
    sub=pairs[pairs.variant=='standard__kick203__eps1e-08__hhalf'].sort_values('offset_xi')
    fig,ax=plt.subplots(figsize=(8,5))
    for col,label in [('left_B_norm','B norm left'),('right_B_norm','B norm right'),('left_delta_norm','delta left'),('right_delta_norm','delta right')]:
        ax.semilogx(sub.offset_xi,sub[col],marker='o',label=label)
    ax.set_xlabel('distance to crossing |delta xi|')
    ax.set_ylabel('observed norm')
    ax.set_title('One-sided probes: allocation discontinuity vs sample refinement')
    ax.legend();fig.tight_layout();fig.savefig(out/'one_sided_event_limits.png',dpi=180);plt.close(fig)
    summary={
       'version':'2026-10-16-Q550-kick203-event-micro-lens-v17.9',
       'source':'v17.5 event microgrid, unchanged 11D history/frame/Shapley',
       'B_crossing_xi':xb,
       'variants':per_variant,
       'prior_grid_matched_overlaps':regressions,
       'classification':'EVENT-ALIGNED OFF-SHELL SHAPLEY ALLOCATION; NOT PHYSICAL KICK IMPULSE OR MATERIAL FLIP',
       'caveat':'B=0 left of crossing may be built into common restart and allocation. The event has a source-defined sidedness. Denser samples do not prove a physical impulse or a continuous birth.',
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    with open(out/'REPORT.md','w') as f:
       f.write('# Q550 kick203 event-micro lens v17.9\n\n')
       f.write('Original v17.5 full11D cocycle and exact symmetric Shapley definition are preserved.\n\n')
       f.write(f'B crossing xi = `{xb:.15f}`\n\n')
       f.write('## Nearest left/right probes\n\n')
       f.write('```text\n'+pd.DataFrame(per_variant).to_string(index=False)+'\n```')
       f.write('\n\n## Exact overlap with v17.5 xi step 0.025\n\n')
       f.write('```text\n'+pd.DataFrame(regressions).to_string(index=False)+'\n```' if regressions else 'No prior reference provided')
       f.write('\n\n## Interpretation boundary\n\n')
       f.write('The Shapley hybrids are off-shell mathematical counterfactuals, not physical trajectories.\n')
       f.write('At the B crossing, plus/minus histories change section sampling side relative to the common restart; a one-sided allocation jump must not be interpreted as a physical shock.\n')
       f.write('The signed frozen reference is fixed within each numerical variant and is not a gauge-independent physical torsion axis.\n')
       f.write('Event crossing is explicitly sampled for diagnostics but exact-side branch selection at floating-point equality may be convention-dependent. Prefer one-sided offsets.\n')
    print('B crossing xi',xb)
    print('Nearest left/right by variant')
    print(pd.DataFrame(per_variant).to_string(index=False))
    print('Matched archival overlap summary')
    print(pd.DataFrame(regressions).to_string(index=False) if regressions else 'not provided')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--allocation-csv',type=Path,required=True)
    p.add_argument('--shapley-npz',type=Path,required=True)
    p.add_argument('--reference-allocation-csv',type=Path)
    p.add_argument('--reference-shapley-npz',type=Path)
    p.add_argument('--out',type=Path,required=True)
    analyze(p.parse_args())
