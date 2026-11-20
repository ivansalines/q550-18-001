#!/usr/bin/env python3
"""Optional original-ODE replay and externally documented structure candidate audit.

Baseline is restarted at the actual v18.6 kick source state and support time,
not reconstructed from a new 400-orbit Q550 integration. This independent replay
requires the original pinned Q550 ODE source files on the server.
"""
from __future__ import annotations
import argparse,hashlib,importlib.util,json,sys,csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(HERE/'reference'),str(HERE/'scripts')]
import q550_state_impulse_response_v13 as old
from analyze_cached_v18_9 import get_data
from structural_math_v18_9 import roots_bracketed


def import_python(path,name):return old.load_module(path,name)


def diff_jacobian(fun,y):
    v=np.asarray(fun(y),float).reshape(-1);n=len(y);out=np.zeros((len(v),n))
    for j in range(n):
        h=1e-6*max(1.,abs(y[j]));q=np.eye(1,n,j).flatten()*h
        out[:,j]=(np.asarray(fun(y+q),float).reshape(-1)-np.asarray(fun(y-q),float).reshape(-1))/(2*h)
    return out


def diff_hessian_on_u(fun,y,u):
    norm=np.linalg.norm(u)
    if norm<1e-100:return np.zeros(len(np.asarray(fun(y)).reshape(-1)))
    h=1e-4*max(1.,np.linalg.norm(y))/norm
    return (np.asarray(fun(y+h*u),float).reshape(-1)-2*np.asarray(fun(y),float).reshape(-1)
            +np.asarray(fun(y-h*u),float).reshape(-1))/(h*h)


def read_provenance(path):
    z=json.loads(Path(path).read_text())
    for key in ('name','physical_reference','coordinate_mapping','metric_or_constraint_formula','derivation_method'):
        if not isinstance(z.get(key),str) or not z[key].strip():raise ValueError('missing structure provenance '+key)
    return z


def save_rows(path,rows):
    if not rows:return
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-root',type=Path,default=old.DEFAULT_BASE)
    p.add_argument('--structure-module',type=Path,default=None)
    p.add_argument('--structure-provenance',type=Path,default=None)
    p.add_argument('--out',type=Path,default=HERE/'results')
    p.add_argument('--max-step-xi',type=float,default=.05)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    if a.structure_module is not None and a.structure_provenance is None:p.error('must provide --structure-provenance')
    data,provenance=get_data();d=data['standard'];xi=d['xi_nodes'];s0=d['s0'];slope=d['slope']
    # Must validate actual model source identity before comparing coordinates.
    names=['audit','support','tube'];reused={}
    for key in names:
        name,h=old.DEPENDENCIES[key];path=a.base_root/name
        if not path.is_file():raise FileNotFoundError(f'original pinned Q550 ODE module unavailable: {path}')
        seen=old.hash_file(path)
        if seen!=h:raise RuntimeError(f'Original {key} source hash mismatch; expected {h} observed {seen}')
        reused[key]={'path':str(path),'sha256':seen}
    audit=import_python(a.base_root/old.DEPENDENCIES['audit'][0],'v189_model_audit')
    support=import_python(a.base_root/old.DEPENDENCIES['support'][0],'v189_model_support')
    tube=import_python(a.base_root/old.DEPENDENCIES['tube'][0],'v189_model_tube')
    params=audit.ModelParameters(gamma_rho=.3);kappa0=(55./6.)*(10./550.)
    def rhs(s,y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=kappa0,params=params,audit=audit),float)
    def jac(s,y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=kappa0,params=params,audit=audit,base=support),float)
    y0=d['source_state_11D'];end=s0+8.*slope
    sol=solve_ivp(rhs,(s0,end),y0,method='DOP853',rtol=1e-12,atol=1e-14,
       max_step=max(1e-12,slope*a.max_step_xi),dense_output=True)
    if not sol.success:raise RuntimeError('actual Q550 ODE solve failed: '+sol.message)
    states=np.asarray(sol.sol(s0+xi*slope),float).T
    ind=np.arange(0,len(xi),8)
    if ind[-1]!=len(xi)-1:ind=np.r_[ind,len(xi)-1]
    errors=[]
    for i in ind:
        s=s0+xi[i]*slope; diff=jac(s,states[i])-d['J_full11D'][i]
        errors.append(float(np.linalg.norm(diff)/max(1e-300,np.linalg.norm(d['J_full11D'][i]))))
    maxJ=max(errors)
    if maxJ>2e-4:raise RuntimeError(f'ODE/J archive mismatch: {maxJ:.3g}; refusing candidate interpretations')
    summary={'version':'v18.9','ode_module_sha256':reused,'used_archived_s0':float(s0),
       'used_archived_source_state':True,'ode_replay_success':True,
       'max_rel_jacobian_archive_regression':maxJ,
       'jacobian_comparison_nodes':len(ind),
       'max_step_xi':a.max_step_xi,
       'source_history':'v18.6 full11D original source state; no 400-orbit replay',
       'structure_input_status':'ABSENT','physical_certification':'NOT_ESTABLISHED'}
    if a.structure_module:
        structure=import_python(a.structure_module,'v189_external_structure')
        prov=read_provenance(a.structure_provenance)
        summary['external_provenance']=prov
        summary['structure_input_status']='EXTERNAL_CANDIDATE_UNVERIFIED'
        rows=[];metric_vals=[];constraint_vals=[];jac_ranks=[];lin_leak=[];quad_leak=[];flow_leak=[]
        for i,(x,y) in enumerate(zip(xi,states)):
            s=s0+x*slope;field=rhs(s,y)
            U=d['P_full11D'][i]@d['kick_direction_11D']
            K=d['quadratic_total_material_11D'][i]
            r={'xi':float(x)}
            if hasattr(structure,'metric'):
                G=np.asarray(structure.metric(y),float)
                if G.shape!=(11,11) or not np.isfinite(G).all() or np.max(abs(G-G.T))>1e-7:
                    raise ValueError('metric must be finite symmetric 11x11')
                lam=np.linalg.eigvalsh(G)
                if lam[0]<=0:raise ValueError('supplied metric is not SPD throughout replay')
                metric_vals.append(G)
                r['candidate_metric_min_eigen']=float(lam[0])
                r['candidate_metric_cond']=float(lam[-1]/lam[0])
                E=d['old_source_material_11D'][i]
                r['candidate_metric_dot_old_total']=float(E@G@K)
            if hasattr(structure,'constraints'):
                Cfun=structure.constraints;C=np.asarray(Cfun(y),float).reshape(-1)
                if C.size==0:raise ValueError('empty constraints array')
                JC=(np.asarray(structure.constraint_jacobian(y),float)
                    if hasattr(structure,'constraint_jacobian') else diff_jacobian(Cfun,y))
                if JC.shape!=(len(C),11):raise ValueError('wrong JC dimension')
                singular=np.linalg.svd(JC,compute_uv=False)
                rank=int(np.sum(singular>1e-9*max(singular[0],1e-30)))
                jac_ranks.append(rank);constraint_vals.append(C)
                L=JC@U;S=JC@K+.5*diff_hessian_on_u(Cfun,y,U)
                flow=JC@field
                lin_leak.append(np.linalg.norm(L));quad_leak.append(np.linalg.norm(S));flow_leak.append(np.linalg.norm(flow))
                r['candidate_constraint_first_variation_residual_norm']=float(np.linalg.norm(L))
                r['candidate_constraint_second_variation_residual_norm']=float(np.linalg.norm(S))
                r['candidate_constraint_rhs_transport_residual_norm']=float(np.linalg.norm(flow))
                r['candidate_constraint_rank']=rank
            rows.append(r)
        if metric_vals:
            N=np.asarray([r['candidate_metric_dot_old_total'] for r in rows])
            Nd=np.gradient(N,xi,edge_order=2)
            summary['candidate_metric_euclidean_unit_aware_roots_xi']=roots_bracketed(xi,N,Nd)
            summary['candidate_metric_condition_max']=float(max(r['candidate_metric_cond'] for r in rows))
        if constraint_vals:
            CV=np.asarray(constraint_vals)
            conservation=float(np.max(np.linalg.norm(CV-CV[0],axis=1)))
            summary['candidate_constraint_conservation_drift_norm']=conservation
            summary['candidate_constraint_jacobian_ranks_unique']=sorted(set(jac_ranks))
            summary['candidate_constraint_first_variation_residual_norm_max']=float(max(lin_leak))
            summary['candidate_constraint_second_variation_residual_norm_max']=float(max(quad_leak))
            summary['candidate_constraint_rhs_transport_residual_norm_max']=float(max(flow_leak))
            summary['warning']='Archived kicks are off-constraint. The correct order-2 constraint residual includes 0.5 D2C[u,u]. Results are NOT constrained physical laws.'
        save_rows(a.out/'candidate_structure_diagnostics.csv',rows)
        summary['physical_certification']='EXTERNAL_STRUCTURE_CANDIDATE_DIAGNOSTICS_ONLY'
    (a.out/'model_replay_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'ode_replay_success':True,'max_rel_jacobian_archive_regression':maxJ,
       'structure_input_status':summary['structure_input_status']},indent=2))

if __name__=='__main__':main()
