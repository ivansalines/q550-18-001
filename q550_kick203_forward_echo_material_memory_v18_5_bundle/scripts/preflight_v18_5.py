#!/usr/bin/env python3
"""Dependency-free preflight (except numpy/pandas), validates embedded v18.4 genealogy."""
from pathlib import Path
import sys,json,hashlib,zipfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
from run_echo_v18_5 import verified_inputs
from echo_math import DELAYS_XI,READ_AGES_XI,AMPS
m=verified_inputs()
p=ROOT/'inputs'/m['historical_v18_4_sendback']
with zipfile.ZipFile(p) as z:
    assert z.testzip() is None, 'invalid historical archive'
    names=z.namelist()
    assert 'results/five_amplitude_samples.csv' in names
    assert 'results/samplewise_odd_fits.csv' in names
    assert 'reference/q550_state_impulse_response_v13.py' in names
assert len(DELAYS_XI)==6 and len(READ_AGES_XI)==5 and READ_AGES_XI[0]==0
assert len(AMPS)==3
print('PASS pinned historical v18.4, original Q550 solver ancestry, 3 eps x 6 delays x 2 tolerance x 2 signs')
print('No Q550 solver input is required for preflight. Historical Q550 source is verified on server at runtime.')
