#!/usr/bin/env python3
"""Read historical Q550 v18.10->v18.9->v18.7->v18.6 chain with CRC and SHA checks."""
from __future__ import annotations
import hashlib,io,json,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'inputs/v18_10_original_sendback.zip'
sys.path.insert(0,str(ROOT/'reference'))
import analyze_cached_v18_9 as reader

def digest(b): return hashlib.sha256(b).hexdigest()

def load():
    with zipfile.ZipFile(INPUT) as z:
        if z.testzip():raise RuntimeError('v18.10 archive CRC failure')
        checks=json.loads(z.read('SHA256SUMS_OUTPUT.json'))['files']
        for n,h in checks.items():
            if digest(z.read(n))!=h:raise RuntimeError('v18.10 member hash mismatch: '+n)
        prev=z.read('inputs/v18_9_original_sendback.zip')
        past=json.loads(z.read('results/summary.json'))
        ode=json.loads(z.read('results/ode_structure_summary.json')) if 'results/ode_structure_summary.json' in z.namelist() else None
    with zipfile.ZipFile(io.BytesIO(prev)) as z:
        if z.testzip():raise RuntimeError('v18.9 archive CRC failure')
        sums=json.loads(z.read('SHA256SUMS_OUTPUT.json'))['files']
        for n,h in sums.items():
            if digest(z.read(n))!=h:raise RuntimeError('v18.9 member hash mismatch: '+n)
        older=z.read('inputs/v18_7_original_sendback.zip')
    transient=ROOT/'.runtime';transient.mkdir(exist_ok=True)
    p=transient/'v18_7_original_sendback.zip';p.write_bytes(older)
    reader.ARCHIVE=p
    data,checked=reader.get_data()
    return data,{'source_v18_10_sha256':digest(INPUT.read_bytes()),
        'v18_10_verified_members':len(checks),'v18_9_sha256':digest(prev),
        'v18_9_verified_members':len(sums),'v18_10_label':past['version'],
        'v18_10_server_ode_replay':None if ode is None else ode.get('replay_success'),
        'historical_validation':checked}
