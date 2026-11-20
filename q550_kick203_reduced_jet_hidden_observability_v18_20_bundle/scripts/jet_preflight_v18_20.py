"""Frozen original-server-Jacobian jet observability audit: NO fresh nonlinear Q550 solves."""
from pathlib import Path
import sys,json,csv,hashlib,os
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from checkpoint_io_v18_20 import read
from frozen_preflight import observers
from jet_math import derivatives,motion_matrices,observability,numerical_rank,null_vectors,canonical,scaled_taylor

def savecsv(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def analyze():
    D,provenance,prior=read();d=D['tight'];dx=np.diff(d['xi_nodes']);h=float(np.median(dx))
    if np.max(abs(dx-h))>1e-10:raise RuntimeError('nonuniform frozen xi grid')
    Hs=observers(d);xi_grid=[2.,3.,4.,5.,6.,7.];ranks=[];comparison=[];checks=[];seed_details=[];null_states={}
    for xi in xi_grid:
        j=int(round(xi/h))
        arrays={k:float(v['slope'])*v['J_full11D'][:,:10,:10] for k,v in D.items()}
        M_by_tol={}
        for tol,Agrid in arrays.items():
            a,a1,a2=derivatives(Agrid,j,h,1);M_by_tol[tol]=motion_matrices(a,a1,a2)
        for k in range(4):
            num=float(np.linalg.norm(M_by_tol['standard'][k]-M_by_tol['tight'][k]))
            comparison.append(dict(xi=xi,jet_order=k,absolute_difference=num,
                relative_difference=num/max(1.,float(np.linalg.norm(M_by_tol['tight'][k])))))
        fine=M_by_tol['tight']
        coarse=[]
        for jump in [1,2,4]:
            Agrid=arrays['tight'];A,A1,A2=derivatives(Agrid,j,h,jump)
            coarse.append(motion_matrices(A,A1,A2))
        for order in (2,3):
            checks.append({'xi':xi,'order':order,
               'step1_step2_rel':float(np.linalg.norm(coarse[0][order]-coarse[1][order])/max(1,np.linalg.norm(coarse[0][order]))),
               'step2_step4_rel':float(np.linalg.norm(coarse[1][order]-coarse[2][order])/max(1,np.linalg.norm(coarse[1][order])))})
        for name,H in Hs.items():
            for tau in [0.05,0.1,0.2]:
                Olist=observability(H,fine,tau=tau)
                for order,O in enumerate(Olist):
                    rk,s=numerical_rank(O)
                    rk2,_=numerical_rank(O,rtol=1e-7)
                    rk3,_=numerical_rank(O,rtol=1e-9)
                    ranks.append(dict(xi=xi,observer=name,tau=tau,jet_order=order,rank=rk,
                      rank_tol_1e7=rk2,rank_tol_1e9=rk3,unseen=10-rk,
                      smax=float(s[0]),smin_nonzero=float(s[rk-1]),
                      smallest_to_largest_if_full=float(s[-1]/s[0]) if rk==10 else 0.,
                      cond_observable=float(s[0]/s[rk-1])))
            if xi==4.:
                Olist=observability(H,fine,tau=0.1)
                for order in range(3):
                    n=null_vectors(Olist[order])
                    if n.shape[1]==0:continue
                    Dnext=H@fine[order+1]@n
                    u,s,vt=np.linalg.svd(Dnext,full_matrices=False)
                    v=canonical(n@vt[0]);dNext=np.linalg.norm(H@fine[order+1]@v)
                    labels={'material3':[0,1,2],'raw4':[0,1],'raw6':[0,1]}
                    if order not in labels.get(name,[]):continue
                    if dNext<1e-7:continue
                    key=name+'_jet'+str(order)
                    defects=[float(np.linalg.norm(H@fine[k]@v)) for k in range(order+1)]
                    seed_details.append({'name':key,'observer':name,'matched_through_order':order,
                      'hidden_dim':int(n.shape[1]),'seed':v.tolist(),
                      'matched_jet_defects':defects,'next_jet_gain':float(dNext)})
                    null_states[key]=(H,v,fine,order)
    # Compare actual ARCHIVED cocycle with Taylor prediction near matching time
    P=d['P_full11D'][:,:10,:10]; i_m=int(round(4/h));Pinv=np.linalg.inv(P[i_m]);trace=[]
    for label,(H,v,M,order) in null_states.items():
        for delta in [-.1,-.05,-.025,.025,.05,.1]:
            j=int(round((4+delta)/h));actual=H@(P[j]@Pinv@v)
            pred=scaled_taylor(H,M,v,delta)
            trace.append({'label':label,'matched_through':order,'delta_xi':delta,
                          'actual_norm':float(np.linalg.norm(actual)),
                          'taylor_3_norm':float(np.linalg.norm(pred)),
                          'abs_taylor_error':float(np.linalg.norm(actual-pred)),
                          'delta_observable_error':float(np.linalg.norm(actual-pred))})
    out=ROOT/'results';out.mkdir(exist_ok=True)
    savecsv(out/'jet_rank_hierarchy.csv',ranks)
    savecsv(out/'jet_tolerance_drift.csv',comparison)
    savecsv(out/'jet_stencil_sensitivity.csv',checks)
    savecsv(out/'jet_taylor_cocycle_check.csv',trace)
    (out/'matched_hidden_seeds.json').write_text(json.dumps(seed_details,indent=2)+'\n')
    at4={r['observer']:[int(next(a['rank'] for a in ranks if a['xi']==4. and a['observer']==r['observer'] and a['tau']==.1 and a['jet_order']==k)) for k in range(4)] for r in ranks if r['xi']==4.}
    if at4!={'material3':[3,6,9,10],'raw4':[4,7,10,10],'raw6':[6,9,10,10]}:
        raise RuntimeError('Unexpected Q550 jet rank sequence: '+str(at4))
    if any(row['rank']!=row['rank_tol_1e7'] or row['rank']!=row['rank_tol_1e9'] for row in ranks):
        raise RuntimeError('Jet numerical ranks unstable across threshold')
    summary={'version':'v18.20','status':'FROZEN_ORIGINAL_SERVER_JACOBIAN_JET_PREFLIGHT_NO_FRESH_NONLINEAR_SOLVES',
        'data_provenance':provenance,'original_v18_19_pair_count':prior['pair_count'],
        'xi_points':xi_grid,'time_scale_tau':[.05,.1,.2],
        'phase_fixed_material_dimension':10,'at_xi4_ranks':at4,
        'max_jet_tol_relative_drift':max(z['relative_difference'] for z in comparison),
        'max_second_jet_step_sensitivity':max(c['step1_step2_rel'] for c in checks if c['order']==2),
        'max_third_jet_step_sensitivity':max(c['step1_step2_rel'] for c in checks if c['order']==3),
        'probe_count':len(seed_details),
        'claim_limits':['Exact-rank statements are numerical local-linear diagnostics, not interval-certified rank.',
                        'Full rank does not imply well-conditioned recovery under measurement noise.',
                        'The 10D material block is phase-driven and is not an autonomous subsystem.',
                        'The measured H matrices are diagnostic projections, not material observables certified on physical constraints.',
                        'No independent hereditary constitutive memory is established.']}
    (out/'frozen_jet_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'status':summary['status'],'ranks':at4,'max_jet_standard_tight_rel':summary['max_jet_tol_relative_drift'],
         'max_third_jet_stencil_rel':summary['max_third_jet_step_sensitivity'],'probe_count':len(seed_details)},indent=2))
    return summary
if __name__=='__main__':analyze()
