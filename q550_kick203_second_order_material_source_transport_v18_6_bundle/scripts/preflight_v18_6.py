#!/usr/bin/env python3
"""Local and server input check, no Q550 orbit recomputation."""
from pathlib import Path
import sys,io,zipfile,numpy as np,json
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_source_transport_v18_6 import verify_history,read_v185,DELAYS,AMPS,TOLS
p=verify_history()
with zipfile.ZipFile(p) as z:
    names=set(z.namelist())
    for tol in TOLS:
        for eps in AMPS:
            for delay in DELAYS:
                name=f'checkpoint/{tol}_eps{eps:g}_delay{delay:g}.npz'
                if name not in names:raise RuntimeError('Missing '+name)
    print('v18.5 archive has all 36 checkpoints')
q=read_v185(p,'standard',AMPS[1],8.)
assert q['source_state_11D'].shape==(11,)
assert q['plus_echo_read_11D'].shape==(5,11)
assert q['minus_echo_read_11D'].shape==(5,11)
assert q['kick_direction_11D'].shape==(11,)
assert np.isclose(np.linalg.norm(q['kick_direction_11D']),1.,rtol=1e-12)
print('Verified v18.5 SHA256:',json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())['input_sha256'])
print('Q550 original modules and orbit are required from their server paths for source generation')
