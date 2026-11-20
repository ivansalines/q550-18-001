#!/usr/bin/env python3
"""Offline manifest verification, source shape, and original Python syntax gates."""
from __future__ import annotations
import io,json,py_compile,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_sign_crossing_v18_7 import verify_history,read_source

def main():
    verified=verify_history()
    for t in ('standard','tight'):
        d=read_source(t)
        if not np.allclose(d['source_bin_edges'],[0,.25,.5,1,2,4,8]):raise RuntimeError('source age-bin gate')
        if d['source_integrand_11D'].shape!=(641,11):raise RuntimeError('source grid shape gate')
        if d['predicted_quadratic_at_delays_11D'].shape!=(6,11):raise RuntimeError('v18.6 anchor shape gate')
        if not np.all(np.isfinite(d['source_integrand_11D'])):raise RuntimeError('nonfinite sources')
    for file in list((ROOT/'reference').glob('*.py'))+list((ROOT/'scripts').glob('*.py')):
        py_compile.compile(str(file),doraise=True)
    print('Preflight OK:',json.dumps(verified))
if __name__=='__main__':main()
