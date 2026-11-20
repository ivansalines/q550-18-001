#!/usr/bin/env python3
"""Zero-ODE validation of embedded v18.1 genealogy, source fingerprints, and event."""
import io,json,hashlib,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
man=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
for label,name in [('v18_1','q550_kick203_second_order_observer_v18_1_sendback.zip'),
                   ('singular_vectors','singular_direction_vectors.csv')]:
    p=ROOT/'inputs'/name
    if not p.exists():raise SystemExit(f'FAIL missing {p}')
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h != man['sha256'][label]:raise SystemExit(f'FAIL SHA256 {p}: {h}')
    print('PASS embedded SHA256',name,h)
p=ROOT/'inputs/q550_kick203_second_order_observer_v18_1_sendback.zip'
with zipfile.ZipFile(p) as z:
    files=z.namelist()
    needed=['second_order_full11D_arrays.npz','q550_kick203_event_micro_v17_9_sendback.zip',
           'q550_kick203_one_sided_observer_v18_0_sendback.zip',
           'run_kick203_second_order_observer_v18_1.py']
    for name in needed:
        got=[a for a in files if a.endswith('/'+name)]
        if len(got)!=1:raise SystemExit(f'FAIL missing/ambiguous nested v18.1 {name}')
    with np.load(io.BytesIO(z.read(next(a for a in files if a.endswith('/second_order_full11D_arrays.npz')))),allow_pickle=False) as data:
        crossing=float(data['crossing_xi'][0]);xi=data['xi'];print('PASS reference samples',len(xi),'xi_B',crossing)
        if len(xi)!=47 or abs(crossing-210.8611366031777)>5e-11:raise SystemExit('FAIL v18.1 event mismatch')
print('PASS v18.1 lineage integrity and scalar crossing')
print('NOTE full material ODE replay requires original external server dependencies in README')
