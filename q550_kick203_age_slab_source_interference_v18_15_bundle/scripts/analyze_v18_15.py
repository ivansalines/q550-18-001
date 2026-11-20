#!/usr/bin/env python3
"""Q550 18.15: age-slab interference and source echo localization; no ODE restart.

The source history is the original cumulative *observer* quantity. Splitting
its age support changes no physical trajectory, and is not a hysteresis test.
"""
import argparse,csv,json,hashlib,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from checkpoint_io_v18_15 import read
from age_slab_math import metrics,orth
from reorientation_math import quantities,roots
WIDTHS=(.5,.25,.125,.0625)

def csv_out(p,rows):
    if not rows:return
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def one(d,tol,baseline):
    x=d['xi_nodes'];A=d['integral_source_frame_11D'][:,:10];g=d['source_integrand_11D'][:,:10];s=float(d['slope']);e=A[80].copy()
    if x.shape!=(641,) or A.shape!=(641,10) or g.shape!=(641,10):
        raise RuntimeError('Q550 source shape incompatible with v18.14 checkpoints')
    if abs(x[80]-1.)>1e-10 or np.linalg.norm(A[0])>1e-10:raise RuntimeError('Invalid reference origin or anchor')
    ai=PchipInterpolator(x,A,axis=0);gi=PchipInterpolator(x,g,axis=0)
    mask=x>1+1e-12
    ev=roots(x[mask],quantities(e,A[mask],g[mask],s)['radial_rate'],'pchip')
    if len(ev)!=4:raise RuntimeError('Expected four v18.14 radial events')
    if max(abs(a-b) for a,b in zip(ev,baseline))>1e-6:raise RuntimeError('Original v18.14 event regression mismatch')
    slabs_out=[];compare=[];pair_out=[];events=[]
    for j,z in enumerate(ev,1):
        variants=[]
        for width in WIDTHS:
            stats,rows,pairs=metrics(ai,e,z,width,gi,s)
            if stats['closure_error']>1e-8 or stats['axial_closure_error']>1e-8 or stats['pair_identity_error']>1e-6:
                raise RuntimeError(f'Age-slab algebraic closure fails at {tol}/{j}/{width}')
            variants.append(stats)
            compare.append({'tolerance':tol,'event':j,**stats})
            for row in rows:slabs_out.append({'tolerance':tol,'event':j,'event_xi':z,'width':width,**row})
            if width==.25:
                pair_out.extend([{'tolerance':tol,'event':j,**pair} for pair in pairs[:15]])
        principal=variants[1] # 0.25 xi, directly comparable to prior ablation
        events.append({'tolerance':tol,'event':j,'xi':float(z),'type':'radial_max' if j%2 else 'radial_min',
          'rho':principal['rho'], 'coarse_negative_mass':variants[0]['negative_axial'],
          'neg_axial_width_025':principal['negative_axial'],
          'neg_axial_width_00625':variants[-1]['negative_axial'],
          'cancellation_width_025':principal['destructive_fraction'],
          'cancellation_width_00625':variants[-1]['destructive_fraction'],
          'old_recent_cos_h2':principal['old_recent_cos'],
          'old_recent_dot_h2':principal['old_recent_dot'],
          'stable_age_horizon_10deg_025':principal['age_horizon_10deg'],
          'stable_age_horizon_10deg_00625':variants[-1]['age_horizon_10deg'],
          'negative_age_centroid_00625':variants[-1]['negative_age_centroid'],
          'old2_norm':principal['old2_norm'],'recent2_norm':principal['recent2_norm'],
          'quad_independent_total_error':principal['independent_source_quadrature_relative_error']})
    return ev,events,compare,slabs_out,pair_out

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'results');args=p.parse_args()
    args.output.mkdir(exist_ok=True,parents=True)
    dataset,prov=read();ref=prov['v18_14_events'];rows=[];cmp=[];slabs=[];pairs=[];roots_all={}
    for tol in ('standard','tight'):
        event,sr,cr,sl,pa=one(dataset[tol],tol,ref)
        roots_all[tol]=event;rows+=sr;cmp+=cr;slabs+=sl;pairs+=pa
    drift=float(np.max(abs(np.asarray(roots_all['standard'])-np.asarray(roots_all['tight']))))
    sens={}
    for e in range(1,5):
        x=[r for r in rows if r['tolerance']=='standard' and r['event']==e][0]
        y=[r for r in rows if r['tolerance']=='tight' and r['event']==e][0]
        sens[str(e)]={'negative_mass_difference_width_00625':abs(x['neg_axial_width_00625']-y['neg_axial_width_00625']),
                      'old_recent_cos_difference':abs(x['old_recent_cos_h2']-y['old_recent_cos_h2'])}
    if drift>1e-5:raise RuntimeError('Standard/tight event mismatch')
    summary={'version':'Q550-v18.15-age-slab-source-interference-echo-localization',
      'status':'SOURCE_AGE_INTERFERENCE_NUMERIC_GAUGE_DECLARED',
      'reference':'frozen v18.14 sendback / v18.13 ODE source genealogy',
      'metric':'initial co-moving Euclidean 10D (not independently physical)',
      'age_partition_widths':WIDTHS,'events_xi':roots_all,
      'by_event_standard':[x for x in rows if x['tolerance']=='standard'],
      'standard_tight_max_root_drift':drift,
      'standard_tight_metrics_differences':sens,'provenance':prov,
      'warnings':[
       'Every age slab is an additive partition of one fixed, already observed co-moving source history; the ODE was not altered or re-integrated.',
       'Norms, dot products, signed reinforcement/opposition, cancellation fraction, angular stabilization horizons and ages of echo depend on declared Euclidean metric.',
       'Age-bin refinement convergence is checked against changes in partition-dependent interference metrics and independently integrated source g; additive closure alone is an algebraic identity, not a convergence proof.',
       'The PCHIP of A and the PCHIP of g are separately archived interpolants: independent quadrature differences measure consistency but need not vanish identically.',
       'A negative signed slab does not mean negative physical energy. No material metric, phase-fixed physical constraint manifold or constitutive hysteresis is certified.',
       'Only four event anchors from the earlier declared metric; no claim of universal max/min family behavior.'
      ]}
    csv_out(args.output/'event_summary.csv',rows)
    csv_out(args.output/'resolution_convergence.csv',cmp)
    csv_out(args.output/'event_age_slabs.csv',slabs)
    csv_out(args.output/'strongest_negative_interference_pairs.csv',pairs)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    report=['# Q550 v18.15 — Age-Slab Source Interference & Echo Localization Audit','',
      '## Test definition','',
      'At each pre-existing four radial event times, the unchanged cumulative source A(xi) in the initial co-moving 10D representation is decomposed by **age** a=xi-u, most recent first:',
      '', 'S_i=A(xi-a_i)-A(xi-a_{i+1}); w_i=Q_e S_i; Q_e=I-ee^T/(e.e), e=A(1).',
      'The sum of all w_i equals w=Q_e A(xi). For v=w/||w||, the signed axial contribution is L_i=v.w_i, the lateral contribution is w_i-L_i*v, and the pairwise interaction is w_i.w_j.',
      'By construction sum L_i=||w|| and sum lateral=0; sign, angles and interference magnitude are Euclidean-gauge dependent.',
      '', '## Event summary (standard tolerance)','',
      '| Event | xi | min/max | old versus recent cosine (h=2) | cancellation (h=.25) | cancellation (h=.0625) | 10-degree stable age horizon (h=.0625) |',
      '|---|---:|---|---:|---:|---:|---:|']
    for r in rows[:4]:
        report.append(f"| {r['event']} | {r['xi']:.6f} | {r['type']} | {r['old_recent_cos_h2']:.3f} | {r['cancellation_width_025']:.3f} | {r['cancellation_width_00625']:.3f} | {r['stable_age_horizon_10deg_00625']:.3f} |")
    report.extend(['','## Controls','',f'Max standard/tight event drift: {drift:.3e} xi.',
      'The CSV `resolution_convergence.csv` measures signed-negative mass, cancellation, angular stabilization and independent source-integrand quadrature for slab sizes 0.5, 0.25, 0.125, 0.0625 xi.',
      'The CSV `strongest_negative_interference_pairs.csv` identifies pairs with negative dot products. These are gauge-specific interference diagnostics, not forces or material energy flows.',
      '', '## Interpretation limits','',
      'This is **not** an ablation of physical memory in the ODE. A different source-history decomposition from the same fixed trajectory is a diagnostic of the observer. Proof of physical hysteresis requires distinct physically realizable histories meeting at matched state and phase.',
      'Age axes refer to the original co-moving xi parameter, not directly to laboratory time.',
      '', '## Source genealogy','',f"Embedded, authenticated original v18.14 sendback SHA256: `{prov['v18_14_original_sha256']}`.",
      'All original v18.14 and nested v18.13 -> v18.6 generation scripts remain inside the preserved original sendback archive.',
      'No 400-orbit recomputation.',''])
    (args.output/'REPORT.md').write_text('\n'.join(report))
    print(json.dumps({'status':summary['status'],'events':[{'i':r['event'],'negative_mass_025':r['neg_axial_width_025'],'negative_mass_00625':r['neg_axial_width_00625'],'old_recent_cos_h2':r['old_recent_cos_h2'],'h10_00625':r['stable_age_horizon_10deg_00625']} for r in rows[:4]],'max_tolerance_event_drift':drift},indent=2))
if __name__=='__main__':main()
