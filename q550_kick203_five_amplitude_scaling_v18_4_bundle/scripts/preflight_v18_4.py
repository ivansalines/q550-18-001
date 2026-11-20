#!/usr/bin/env python3
"""Local archive/grid/provenance checks; no Q550 original solver needed."""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
from run_scan_v18_4 import verify_source,read_bridge,NEW_AMPS
from scaling_math import AMPS

p=verify_source();z=read_bridge(p)
assert len(z['xi'])==47 and np.all(np.diff(z['xi'])>0)
assert z['centered_stencil_clean'].astype(bool).sum()==36
assert len(AMPS)==5 and all(e not in (1e-8,3e-8) for e in NEW_AMPS)
for tol in ('standard','tight'):
    for eps in (1e-8,3e-8):
        pre=f'{tol}__eps{eps:g}'
        for field in ('centered_delta_11D','centered_interaction_11D','centered_factor_secants_3x11','order2_delta_11D'):
            assert pre+'__'+field in z,(pre,field)
print('PASS: pinned v18.3 upstream ZIP, 47 xi coordinates, 5 eps levels, both tolerance anchors, frozen axis')
print('New ODE integrations: six +/- history pairs (twelve perturbed ODE branches), three new eps x two tolerances.')
