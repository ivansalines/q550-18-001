#!/usr/bin/env python3
"""Frozen v18.7 Q550 11D metric-sensitivity audit; optional EXTERNAL constraint Jacobians.

Uses original archived 11D source and cocycle.  No new Q550 ODE solves.
No physically justified state metric/constraint functions are inferred from arrays.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,sys,zipfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from metric_math_v18_8 import (VERSION, metric_dot, projectors_from_jacobians,
  fixed_metric_dot_and_rate,source_metric_dot_and_rate,roots_from_samples,
  fixed_w_from_reference,mixed_metric,alpha_crossing_extinction)
ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=ROOT/'inputs/v18_7_original_sendback.zip'
RES=ROOT/'results'
TOLS=['standard','tight'];REFS=[1.,2.,4.,8.]
ALPHAS=[0,.05,.10,.15,.20,.25,.30,.40,.50,.75,1.]
CONTRASTS=[1.1,1.5,2.,4.,10.]
DRAWS=200


def sha(b):return hashlib.sha256(b).hexdigest()

def wr_csv(path,rows):
    if not rows:return
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def load_sources():
    if not ORIGINAL.exists():raise FileNotFoundError(ORIGINAL)
    with zipfile.ZipFile(ORIGINAL) as z:
        bad=z.testzip()
        if bad:raise RuntimeError(f'original zip is corrupt: {bad}')
        manifest=json.loads(z.read('SHA256SUMS_OUTPUT.json'))['files']
        for file,expected in manifest.items():
            if sha(z.read(file))!=expected:raise RuntimeError(f'original hash mismatch {file}')
        d={}
        with zipfile.ZipFile(io.BytesIO(z.read('inputs/v18_6_original_sendback.zip'))) as s:
            for tol in TOLS:
                with np.load(io.BytesIO(s.read(f'checkpoint/{tol}_source_transport.npz')),allow_pickle=False) as a:
                    slope=float(a['source_scale'][0])
                with np.load(io.BytesIO(z.read(f'checkpoint/{tol}_sign_crossing.npz')),allow_pickle=False) as a:
                    d[tol]={key:a[key].copy() for key in a.files}
                d[tol]['slope']=slope
                d[tol]['original_roots']=np.asarray(json.loads(z.read(f'results/{tol}_crossing_summary.json'))['roots'],object)
    return d,{'input_sha256':sha(ORIGINAL.read_bytes()),'original_member_hashes_verified':len(manifest)}


def ensemble_series(E,K,J,P,g,x,slope,lim,draws,seed):
    rng=np.random.default_rng(seed)
    ep=slope*np.einsum('nij,nj->ni',J,E)
    kp=slope*(np.einsum('nij,nj->ni',J,K)+np.einsum('nij,nj->ni',P,g))
    rows=[]
    for ratio in lim:
        for i in range(draws):
            w=np.exp(0.5*np.log(ratio)*rng.uniform(-1,1,11))
            n=np.einsum('ni,i,ni->n',E,w,K)
            nd=np.einsum('ni,i,ni->n',ep,w,K)+np.einsum('ni,i,ni->n',E,w,kp)
            roots=roots_from_samples(x,n,nd)
            rows.append({'contrast_limit':ratio,'draw_index':i,'n_roots':len(roots),
                         'roots_xi_json':json.dumps([round(v,9) for v in roots]),
                         'actual_weight_ratio':float(max(w)/min(w)),
                         'weights_json':json.dumps([float(v) for v in w])})
    return rows


def constraints_run(input_path,provenance_path,primary):
    if not input_path: return None,[]
    if not provenance_path:raise ValueError('--constraint-provenance required with --constraint-npz')
    prov=json.loads(Path(provenance_path).read_text())
    for need in ['constraint_name','constraint_formula','coordinate_mapping','jacobian_construction','source_or_reference']:
        if not prov.get(need):raise ValueError(f'constraint provenance requires nonempty {need}')
    with np.load(input_path,allow_pickle=False) as d:
        xi=d['xi_nodes'];JC=d['constraint_jacobians']
    x=primary['xi_nodes']
    if xi.shape!=x.shape or np.max(abs(xi-x))>1e-10:raise ValueError('constraint xi_nodes mismatch')
    projector,rank,gaps=projectors_from_jacobians(JC)
    ident=np.eye(11)
    rows=[];summary=[]
    for tol in TOLS:
        data=primary if tol=='standard' else constraints_run.tol_data
        E=data['old_source_material_11D'];K=data['quadratic_total_material_11D']
        J=data['J_full11D'];P=data['P_full11D'];g=data['source_integrand_11D'];slope=data['slope']
        leakageE=[];leakageK=[]
        for i in range(80,len(x)):
            leakE=np.linalg.norm((ident-projector[i])@E[i])/max(np.linalg.norm(E[i]),1e-20)
            leakK=np.linalg.norm((ident-projector[i])@K[i])/max(np.linalg.norm(K[i]),1e-20)
            leakageE.append(leakE);leakageK.append(leakK)
        for eta in [0.,.01,.1,1.]:
            W=projector+eta*(ident-projector)
            values=metric_dot(E,K,W)
            # Time-varying constraints metric: derivative evaluated numerically
            # from W'(xi) plus analytic E',K', not pretending W is constant.
            Wprime=np.gradient(W,x,axis=0,edge_order=2)
            ep=slope*np.einsum('nij,nj->ni',J,E)
            kp=slope*(np.einsum('nij,nj->ni',J,K)+np.einsum('nij,nj->ni',P,g))
            dv=metric_dot(ep,K,W)+metric_dot(E,kp,W)+metric_dot(E,K,Wprime)
            roots=roots_from_samples(x,values,dv)
            rows.append({'tolerance':tol,'metric':f'projected_eta_{eta:g}', 'n_roots':len(roots),
                         'roots_xi_json':json.dumps([round(v,9) for v in roots]),
                         'max_early_normal_leakage':float(max(leakageE)),
                         'max_total_normal_leakage':float(max(leakageK)),
                         'median_early_normal_leakage':float(np.median(leakageE)),
                         'median_total_normal_leakage':float(np.median(leakageK))})
        summary.append({'tolerance':tol,'constraint_jacobian_rank':rank,
                        'max_early_normal_leakage':float(max(leakageE)),
                        'max_total_normal_leakage':float(max(leakageK))})
    result={'status':'EXTERNAL_UNVERIFIED_CONSTRAINT_JACOBIANS','input_sha256':sha(Path(input_path).read_bytes()),
            'provenance_sha256':sha(Path(provenance_path).read_bytes()),'supplied_provenance':prov,
            'rank':rank,'diagnostics':summary,'warning':'This is a first-order projection of OFF-CONSTRAINT second-order source data; NOT an on-constraint second variation. Nonzero normal components can be demanded by constraint curvature.'}
    return result,rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--constraint-npz',type=Path)
    ap.add_argument('--constraint-provenance',type=Path)
    ap.add_argument('--draws',type=int,default=DRAWS)
    args=ap.parse_args()
    RES.mkdir(parents=True,exist_ok=True)
    d,integrity=load_sources()
    all_rows=[];ens_rows=[];q={}
    figures={}
    for tol in TOLS:
        a=d[tol];xi=a['xi_nodes'];P=a['P_full11D'];J=a['J_full11D'];g=a['source_integrand_11D'];
        E=a['old_source_material_11D'];K=a['quadratic_total_material_11D'];
        T=a['integral_source_frame_11D'];slope=a['slope'];early=T[80]
        assert xi.shape==(641,) and P.shape==(641,11,11)
        if np.max(abs(E-np.einsum('nij,j->ni',P,early)))>2e-10:raise RuntimeError('old source full11D fails reconstruction')
        if np.max(abs(K-np.einsum('nij,nj->ni',P,T)))>2e-10:raise RuntimeError('total source full11D fails reconstruction')
        lab,ld=fixed_metric_dot_and_rate(E,K,J,P,g,slope,np.eye(11))
        src,sd=source_metric_dot_and_rate(early,T,g,slope)
        old_roots=[float(r['xi_root']) for r in a['original_roots']]
        check_roots=roots_from_samples(xi,lab,ld)
        if len(check_roots)!=3 or max(abs(u-v) for u,v in zip(check_roots,old_roots))>5e-8:raise RuntimeError('v18.7 roots not reproduced')
        if np.min(src[xi>=1])<=0:raise RuntimeError('source-frame dot sign was expected to be positive')
        rootsets={}; lobs=alpha_crossing_extinction(xi,lab,src)
        for alpha in ALPHAS:
            N,D=mixed_metric(lab,ld,src,sd,alpha)
            roots=roots_from_samples(xi,N,D)
            name=f'mixed_pullback_alpha_{alpha:g}'
            rootsets[name]=roots
            all_rows.append({'tolerance':tol,'metric_family':'dynamic_pullback_mix','metric':name,
                             'n_roots':len(roots),'roots_xi_json':json.dumps([round(v,11) for v in roots]),
                             'root_fidelity':'Cubic-Hermite reconstruction on grid, not interval-certified'})
        for ref in REFS:
            W=fixed_w_from_reference(P[int(round(ref/.0125))]);N,D=fixed_metric_dot_and_rate(E,K,J,P,g,slope,W)
            roots=roots_from_samples(xi,N,D)
            name=f'fixed_reference_whitening_xi_{ref:g}'
            rootsets[name]=roots
            all_rows.append({'tolerance':tol,'metric_family':'constant_whitened','metric':name,
                             'n_roots':len(roots),'roots_xi_json':json.dumps([round(v,11) for v in roots]),
                             'root_fidelity':'Cubic-Hermite reconstruction on grid, not interval-certified'})
        # Coordinate-suppression diagnostics are NOT constraints. They localize
        # sensitivity to each recorded channel without claiming a physical normal.
        for channel in range(11):
            for eta in [.01,.1]:
                w=np.ones(11);w[channel]=eta
                N,D=fixed_metric_dot_and_rate(E,K,J,P,g,slope,np.diag(w))
                roots=roots_from_samples(xi,N,D)
                all_rows.append({'tolerance':tol,'metric_family':'coordinate_ablation_not_constraint',
                  'metric':f'channel_{channel:02d}_weight_{eta:g}',
                  'n_roots':len(roots),'roots_xi_json':json.dumps([round(v,11) for v in roots]),
                  'root_fidelity':'Cubic-Hermite reconstruction on grid, not interval-certified'})
        en=ensemble_series(E,K,J,P,g,xi,slope,CONTRASTS,args.draws,184)
        for item in en:ens_rows.append({'tolerance':tol,**item})
        count={str(r):{str(c):sum(rr['contrast_limit']==r and rr['n_roots']==c for rr in en)
                        for c in sorted({rr['n_roots'] for rr in en if rr['contrast_limit']==r})} for r in CONTRASTS}
        q[tol]={'original_roots_xi':old_roots,'reconstructed_roots_xi':check_roots,
                'root_max_absolute_replay_discrepancy_xi':float(max(abs(u-v) for u,v in zip(check_roots,old_roots))),
                'inverse_cocycle_positive_on_all_561_nodes':bool(np.all(src[xi>=1]>0)),
                'pullback_root_extinction_thresholds_approx':lobs,
                'ensemble_distributions':count,
                'fixed_reference_counts':{name:len(roots) for name,roots in rootsets.items() if name.startswith('fixed_reference')},
                'condition_P_max':float(max(np.linalg.cond(p) for p in P[::10]))}
        figures[tol]={'xi':xi,'lab':lab,'pullback':src,
                       'mix10':mixed_metric(lab,ld,src,sd,.1)[0],
                       'mix25':mixed_metric(lab,ld,src,sd,.25)[0],
                       'mix50':mixed_metric(lab,ld,src,sd,.5)[0],
                       'roots':old_roots}
    constraints_run.tol_data=d['tight']
    cs,csrows=constraints_run(args.constraint_npz,args.constraint_provenance,d['standard'])
    wr_csv(RES/'metric_root_sensitivity.csv',all_rows)
    wr_csv(RES/'diagonal_metric_ensemble.csv',ens_rows)
    if csrows:wr_csv(RES/'constraint_projection_diagnostic.csv',csrows)
    summary={'version':VERSION,'constraint_status':('UNVERIFIED_EXTERNAL_JACOBIAN_PROJECTION' if cs else 'OFF_CONSTRAINT_DIAGNOSTIC'),
             'original_integrity':integrity,'by_tolerance':q,'constraint_projection':cs,
             'interpretation_limits':['Metric contrast changes the 11D observer: no physical 11D metric is known.',
                'The dynamical pullback metric uses the cocycle and is NOT a fixed mechanical/energy metric.',
                'Diagonal reweighting experiments are coordinate-unit sensitivities, not physical perturbations.',
                'The constrained protocol requires genuine constraint equations plus on-manifold 1st/2nd variation and independent validation.',
                'Grid-resolved crossings are not an interval-certified count of all roots.']}
    (RES/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    fig,ax=plt.subplots(figsize=(9,5)); f=figures['standard'];sel=f['xi']>=1
    for field,label in [('lab','Euclidean'),('mix10','10% pulled-back'),('mix25','25% pulled-back'),('mix50','50% pulled-back')]:
        vals=f[field][sel]; ax.plot(f['xi'][sel],vals/max(float(np.max(abs(f['lab'][sel]))),1e-12),label=label)
    ax.axhline(0,lw=1,linestyle=':'); ax.set(xlabel='Delay from kick 203 (xi)',ylabel='Old × total, normalized to Euclidean peak',title='Q550 v18.8: sign crossings depend on metric');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(RES/'mixed_metric_signs.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4.5));
    for i,r in enumerate(CONTRASTS):
        counts=[sum(row['tolerance']=='standard' and row['contrast_limit']==r and row['n_roots']==k for row in ens_rows) for k in [0,1,2,3,4,5,6,7,8]]
        # Total share of metrics retaining exactly three roots.
        ax.bar(i,float(counts[3])/args.draws,width=.65)
    ax.set_xticks(range(len(CONTRASTS)),[str(v) for v in CONTRASTS]);ax.set_ylim(0,1.04)
    ax.set(xlabel='Maximum possible ratio of coordinate weights',ylabel='Fraction with exactly 3 resolved crossings',title=f'Deterministic diagonal-metric scan (n={args.draws} per scale)');ax.grid(alpha=.2,axis='y');fig.tight_layout();fig.savefig(RES/'diagonal_metric_robustness.png',dpi=160);plt.close(fig)
    lines=['# Q550 v18.8 — Material Metric & Constraint Sensitivity Audit','',
      'Measured from the user-provided v18.7 full11D checkpoints; no new Q550 ODE solve.','',
      '## Falsifiable result','',
      '- Euclidean (11D state-coordinate) metric retains exactly three sign-crossing roots.',
      '- Continuous positive-definite pullback mixture removes the first pair, then the final root, as alpha increases.',
      '- A constant positive-definite metric obtained by whitening at a different reference xi may yield a different number of roots.',
      '- All these experiments measure observer/metric dependence, **not** a physical material metric.',
      '', '## Data', '', '| Metric | standard roots (grid-resolved) |', '|---|---|']
    for row in all_rows:
        if row['tolerance']=='standard' and (row['metric'].startswith('mixed_pullback_alpha_') or row['metric'].startswith('fixed_reference')):
            lines.append(f"| {row['metric']} | {row['n_roots']} |")
    lines+=['','## Diagonal scale sensitivity','',
      f"Deterministic seed 184, {args.draws} draws per maximum weight ratio. No ensemble constitutes a probability distribution over plausible physical metrics.", '']
    for r,dist in q['standard']['ensemble_distributions'].items(): lines.append(f'- Maximum ratio {r}: '+', '.join(f'{num} crossings={cnt}' for num,cnt in dist.items()))
    lines+=['','## What is and is not constrained','',
      f"Constraint input: {'externally supplied but NOT physically certified' if cs else 'NOT SUPPLIED (off-constraint diagnostic)' }.",
      'Optional user-supplied rank-stable constraint Jacobians produce tangent/normal leakage and projected-metric diagnostics. A projector onto ker(DC) is NOT a second-order on-manifold construction: the true second variation can have normal curvature prescribed by D²C.',
      '', 'The original v18.7 source data came from unconstrained kicks. Even a genuine DC cannot retroactively make those trajectories constrained.',
      '', '## Reproduction','', 'Run `bash scripts/run_on_server.sh` from the bundle root (no server-specific Q550 files needed).',
      'For optional constraint arrays see `templates/CONSTRAINT_INPUT.md`. The original v18.7 archive retains scripts, inputs, data, manifests and complete v18.6/v18.5 history.',
      '', '## Numeric qualifications','', 'Cubic Hermite reconstruction resolves roots bracketed by the 0.0125-xi grid. Closely spaced, unbracketed pairs may be missed. Mixed-metric annihilation thresholds are interpolated diagnostics, not rigorous bifurcation certificates.']
    (RES/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'version':VERSION,'summary':str(RES/'summary.json'), 'euclidean_roots_standard':q['standard']['reconstructed_roots_xi'],
        'metric_extinction_thresholds':q['standard']['pullback_root_extinction_thresholds_approx'],
        'constraint_status':summary['constraint_status']},indent=2))

if __name__=='__main__':main()
