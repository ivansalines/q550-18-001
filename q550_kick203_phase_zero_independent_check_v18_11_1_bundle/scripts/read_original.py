#!/usr/bin/env python3
"""Independent nested historical archive reader, not a wrapper of the archived reader."""
from __future__ import annotations
import hashlib, io, json, zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'inputs/v18_11_original_sendback.zip'
def sha(b):return hashlib.sha256(b).hexdigest()
def verified(archive,lab):
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        broken=z.testzip()
        if broken:raise RuntimeError(f'{lab}: corrupt ZIP member {broken}')
        manifests=json.loads(z.read('SHA256SUMS_OUTPUT.json'))['files']
        for name,digest in manifests.items():
            if sha(z.read(name))!=digest:raise RuntimeError(f'{lab}: SHA256 mismatch {name}')
        return {n:z.read(n) for n in z.namelist() if n in ['inputs/v18_10_original_sendback.zip','inputs/v18_9_original_sendback.zip','inputs/v18_7_original_sendback.zip','inputs/v18_6_original_sendback.zip']},len(manifests)
def load():
    b=SOURCE.read_bytes()
    archived=[]
    for label,key in [('v18.11','inputs/v18_10_original_sendback.zip'),('v18.10','inputs/v18_9_original_sendback.zip'),('v18.9','inputs/v18_7_original_sendback.zip'),('v18.7','inputs/v18_6_original_sendback.zip')]:
        parts,count=verified(b,label)
        archived.append({'generation':label,'sha256':sha(b),'verified_members':count})
        before=b
        b=parts[key]
        if label=='v18.7': v18_7_bytes=before
    _,count=verified(b,'v18.6')
    archived.append({'generation':'v18.6','sha256':sha(b),'verified_members':count})
    data={}
    with zipfile.ZipFile(io.BytesIO(v18_7_bytes)) as z7,zipfile.ZipFile(io.BytesIO(b)) as z6:
        for tol in ('standard','tight'):
            with np.load(io.BytesIO(z7.read(f'checkpoint/{tol}_sign_crossing.npz')),allow_pickle=False) as z:
                v={k:z[k].copy() for k in ('xi_nodes','J_full11D','P_full11D')}
            with np.load(io.BytesIO(z6.read(f'checkpoint/{tol}_source_transport.npz')),allow_pickle=False) as z:
                v.update({'slope':float(z['source_scale'][0]),'s0':float(z['source_s'][0]),'source_state_11D':z['source_state_11D'].copy()})
            data[tol]=v
    return data,archived
