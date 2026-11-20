#!/usr/bin/env python3
"""Independently propagate 11D fundamental matrix using archived J only."""
from __future__ import annotations
import csv,json,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from read_original import load

def block_defects(J,P):
    return {'jac_last_row_max':float(np.max(np.abs(J[:,10,:]))),
            'cocycle_last_row_first10_max':float(np.max(np.abs(P[:,10,:10]))),
            'cocycle_last_corner_minus_one_max':float(np.max(np.abs(P[:,10,10]-1.))),
            'jac_upper_right_max_column_norm':float(np.max(np.linalg.norm(J[:,:10,10],axis=1))),
            'cocycle_upper_right_max_column_norm':float(np.max(np.linalg.norm(P[:,:10,10],axis=1)))}

def rk4_from_j_only(x,J,slope,substeps):
    cubic=CubicSpline(x,J,axis=0)
    matrix=np.eye(11);output=[matrix.copy()]
    for i in range(len(x)-1):
        h=(x[i+1]-x[i])/substeps
        for k in range(substeps):
            t=x[i]+k*h
            L1=slope*cubic(t)@matrix
            L2=slope*cubic(t+h/2)@(matrix+h*L1/2)
            L3=slope*cubic(t+h/2)@(matrix+h*L2/2)
            L4=slope*cubic(t+h)@(matrix+h*L3)
            matrix=matrix+h*(L1+2*L2+2*L3+L4)/6
        output.append(matrix.copy())
    return np.array(output)

def audit():
    data,history=load();rows=[];allrows=[]
    for tol,d in data.items():
        x=d['xi_nodes'];J=d['J_full11D'];P=d['P_full11D'];s=d['slope']
        if x.shape!=(641,) or J.shape!=(641,11,11) or P.shape!=(641,11,11):raise RuntimeError('Unexpected arrays: '+tol)
        if abs(np.max(np.diff(x))-.0125)>1e-8:raise RuntimeError('Unexpected xi grid')
        defects=block_defects(J,P)
        if max(defects[k] for k in ('jac_last_row_max','cocycle_last_row_first10_max','cocycle_last_corner_minus_one_max'))>1e-12:raise RuntimeError('Structural zeros violated '+tol)
        estimates=[]
        for sub in [1,2,4]:
            reconstruction=rk4_from_j_only(x,J,s,sub)
            residual=np.linalg.norm(reconstruction-P,axis=(1,2))/np.maximum(1.,np.linalg.norm(P,axis=(1,2)))
            qrel=np.linalg.norm((reconstruction-P)[:,:10,10],axis=1)/np.maximum(1.,np.linalg.norm(P[:,:10,10],axis=1))
            for xi,err,qerr in zip(x,residual,qrel):
                allrows.append({'tolerance':tol,'substeps':sub,'xi':float(xi),'full_rel_error':float(err),'q_rel_error':float(qerr)})
            estimates.append({'substeps_per_interval':sub,'max_full_cocycle_rel_error':float(np.max(residual)),
                              'terminal_full_cocycle_rel_error':float(residual[-1]),'max_driven_phase_column_rel_error':float(np.max(qrel)),
                              'max_reconstructed_phase_row_defect':float(max(np.max(np.abs(reconstruction[:,10,:10])),np.max(np.abs(reconstruction[:,10,10]-1))))})
        order_ratios=[estimates[k]['max_full_cocycle_rel_error']/estimates[k+1]['max_full_cocycle_rel_error'] for k in range(2)]
        record={'tolerance':tol,'block':defects,'integration':estimates,'refinement_error_ratios':order_ratios,'end_P_condition':float(np.linalg.cond(P[-1]))}
        if estimates[-1]['max_full_cocycle_rel_error']>5e-9 or min(order_ratios)<10:raise RuntimeError('Independent cocycle regression: '+tol+' '+str(record))
        rows.append(record)
    report={'version':'v18.11.1','kind':'INDEPENDENT_FROZEN_J_TO_P_PROPAGATION','original_history':history,'results':rows,
       'frozen_gate':'PASSED','nonlinear_gate':'PENDING_SERVER_OR_SEPARATE_REPORT',
       'physical_constraint_claim':False,'physical_metric_claim':False}
    out=ROOT/'results';out.mkdir(exist_ok=True)
    (out/'frozen_summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    with (out/'frozen_profile.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
    print(json.dumps({'frozen_gate':'PASSED','by_tolerance':[{ 'tol':a['tolerance'],'errors':[b['max_full_cocycle_rel_error'] for b in a['integration']],'refinement_error_ratios':a['refinement_error_ratios']} for a in rows]},indent=2))
if __name__=='__main__':audit()
