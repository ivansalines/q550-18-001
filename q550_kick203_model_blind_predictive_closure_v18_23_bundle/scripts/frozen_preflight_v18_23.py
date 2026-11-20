"""Authentic frozen-tangent PRE-FLIGHT ONLY. Not a nonlinear Q550 server run."""
from pathlib import Path
import sys,zipfile,json
import numpy as np
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
from checkpoint_io_v18_23 import read
from frozen_preflight import observers
from checkpoint_io_v18_21 import PREFIX as PREFIX20
from blind_math import HISTORY,HORIZONS
from blind_experiment import initial_directions,run_experiment,OBSERVERS

def main():
    D,prov,proof=read();d=D['tight'];Hs=observers(d)
    with zipfile.ZipFile(ROOT/'.runtime/v18_20_original_server_sendback.zip') as z:
        records=json.loads(z.read(PREFIX20+'results/matched_hidden_seeds.json'))
    hidden={p['name']:p['seed'] for p in records}
    seeds=initial_directions(hidden)
    xi=np.asarray(d['xi_nodes'],float)
    full=np.asarray(d['P_full11D'][:,:10,:10],float)
    sp=CubicSpline(xi,full,axis=0)
    phi4=sp(4.)
    if np.linalg.cond(phi4)>1e7:raise RuntimeError('Bad propagator conditioning')
    B=sp(np.r_[4+HISTORY,4+HORIZONS])@np.linalg.inv(phi4)
    assert B.shape==(len(HISTORY)+len(HORIZONS),10,10)
    datasets={}
    for key in OBSERVERS:
        H=Hs[key];data={}
        for split,U in seeds[key].items():
            if split in ('matched_seed','provenance'):continue
            reduced=np.einsum('ab,tbc,nc->nta',H,B,U)
            x=reduced[:,:len(HISTORY),:]
            y=reduced[:,len(HISTORY):,:].reshape(len(U),-1)
            data[split]=(x,y)
        datasets[key]=data
    result=run_experiment(datasets,'FROZEN_ORIGINAL_Q550_TANGENT_SENSOR_ONLY_PREFLIGHT',
       {'original_v18_22_sha256':prov['v18_22_server_sha256'],
        'source_nonlinear_status':proof['status'],'construction':'Propagator used ONLY to create samples; estimator sensor-only'},False)
    if not result['predictor_access_to_ODE_Jacobian_cocycle_or_true_state'] is False:raise RuntimeError('Model access violation')
if __name__=='__main__':main()
