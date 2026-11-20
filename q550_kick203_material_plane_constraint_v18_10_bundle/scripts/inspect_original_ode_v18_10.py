#!/usr/bin/env python3
"""Pinned-original Q550 ODE: dynamical coordinate inventory, not a physical certificate.

Uses v18.6's archived actual kick state, does not regenerate baseline orbit.
No C(y), physical G(y) or material projection is guessed.
"""
from __future__ import annotations
import argparse,ast,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
import q550_state_impulse_response_v13 as src
from checkpoint_io_v18_10 import load


def defs(path,targets):
    """AST inventory of actual declarations, not semantic equation extraction."""
    parsed=ast.parse(path.read_text()); found=[]
    for n in ast.walk(parsed):
        if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name in targets:
            found.append({'name':n.name,'line':n.lineno,'end_line':n.end_lineno,
              'arguments':[a.arg for a in n.args.args],
              'lines_of_source':n.end_lineno-n.lineno+1})
    return found


def main():
    p=argparse.ArgumentParser();p.add_argument('--base-root',type=Path,default=src.DEFAULT_BASE)
    p.add_argument('--output',type=Path,default=ROOT/'results')
    p.add_argument('--max-step-xi',type=float,default=.05)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    loaded,history=load();base=loaded['standard'];x=base['xi_nodes'];s0=base['s0'];slope=base['slope'];
    pinned={};target={'audit':['ModelParameters'],'support':['support_rhs_float'],'tube':['float_support_jacobian']}
    for k in ('audit','support','tube'):
        name,sha=src.DEPENDENCIES[k];path=a.base_root/name
        if not path.is_file():raise FileNotFoundError('Missing original pinned Q550 source '+str(path))
        got=src.hash_file(path)
        if got!=sha:raise RuntimeError(f'Original Q550 ODE drift detected for {path}: expected {sha}, got {got}')
        pinned[k]={'path':str(path),'sha256':sha,'function_inventory':defs(path,target[k])}
    audit=src.load_module(a.base_root/src.DEPENDENCIES['audit'][0],'v1810_ode_audit')
    support=src.load_module(a.base_root/src.DEPENDENCIES['support'][0],'v1810_ode_support')
    tube=src.load_module(a.base_root/src.DEPENDENCIES['tube'][0],'v1810_ode_tube')
    params=audit.ModelParameters(gamma_rho=.3);kappa0=(55./6.)*(10./550.)
    def rhs(s,y):return np.asarray(support.support_rhs_float(np.asarray(y,float),kappa0=kappa0,params=params,audit=audit),float)
    def jac(s,y):return np.asarray(tube.float_support_jacobian(np.asarray(y,float),kappa0=kappa0,params=params,audit=audit,base=support),float)
    sol=solve_ivp(rhs,(s0,s0+8*slope),base['source_state_11D'],method='DOP853',
      rtol=1e-12,atol=1e-14,max_step=slope*a.max_step_xi,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    at=s0+x*slope; states=sol.sol(at).T
    ind=np.r_[np.arange(0,len(x),8),len(x)-1];ind=np.unique(ind)
    fs=[];jd=[];errors=[];residual_rows=[]
    for i in ind:
        y=states[i];f=rhs(at[i],y);j=jac(at[i],y)
        fs.append(f);jd.append(j)
        delta=j-base['J_full11D'][i]
        errors.append(np.linalg.norm(delta)/max(1e-300,np.linalg.norm(base['J_full11D'][i])))
        residual_rows.append({'xi':float(x[i]),'jacobian_rel_archive':float(errors[-1]),
          'state_norm':float(np.linalg.norm(y)),'rhs_norm':float(np.linalg.norm(f))})
    if max(errors)>2e-4:raise RuntimeError(f'Archived ODE mismatch {max(errors)}; refuse analysis')
    fs=np.asarray(fs);jd=np.asarray(jd)
    row_norm=np.max(np.linalg.norm(jd,axis=2),axis=0)
    col_norm=np.max(np.linalg.norm(jd,axis=1),axis=0)
    rows=[]
    for j in range(11):
        variation=float(np.max(fs[:,j])-np.min(fs[:,j]))
        rows.append({'state_index_zero_based':j,'max_jacobian_row_norm':float(row_norm[j]),
          'max_jacobian_column_norm':float(col_norm[j]),
          'rhs_min':float(np.min(fs[:,j])),'rhs_max':float(np.max(fs[:,j])),
          'rhs_variation':variation,'identically_zero_jacobian_row_numeric_test':bool(row_norm[j]<1e-11),
          'rhs_identically_zero_on_sample_numeric_test':bool(np.max(abs(fs[:,j]))<1e-11),
          'interpretation':'DIAGNOSTIC_ONLY_NO_AUTOMATIC_CONSTRAINT'})
    with (a.output/'ode_coordinate_inventory.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (a.output/'ode_jacobian_replay.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(residual_rows[0]));w.writeheader();w.writerows(residual_rows)
    result={'version':'v18.10','pinned_ode_source':pinned,'replay_success':True,
     'uses_archived_kick203_state':True,'jacobian_comparison_nodes':len(ind),
     'max_rel_archive_jacobian_difference':float(max(errors)),
     'row_norm_near_zero_indices':[int(i) for i in np.where(row_norm<1e-11)[0]],
     'rhs_near_zero_indices':[int(i) for i in np.where(np.max(abs(fs),axis=0)<1e-11)[0]],
     'ode_definition_only_has_no_verified_material_metric':True,
     'verified_nontrivial_constraint_function':False,
     'coordinate_mapping_verified_from_model_definitions':False,
     'constraint_projected_source_generated':False,
     'warning':'State-coordinate inventory from Jacobian and RHS samples is NOT proof of a physical constraint, first integral, energy, mechanical metric, or true body plane.'}
    (a.output/'ode_structure_summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'ode_replay_success':True,'max_rel_jacobian':max(errors),
       'jacobian_zero_row_indices':result['row_norm_near_zero_indices'],
       'physical_metric_certified':False},indent=2))

if __name__=='__main__':main()
