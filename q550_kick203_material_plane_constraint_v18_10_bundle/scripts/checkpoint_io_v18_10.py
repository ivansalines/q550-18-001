#!/usr/bin/env python3
"""Validate and read the complete v18.9 -> v18.7 -> v18.6 genealogy.

The v18.7 loader has been byte-preserved from v18.9 and checks all of its
historical SHA256 entries. No external data and no interpolation are required.
"""
from __future__ import annotations
import hashlib,io,json,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'inputs/v18_9_original_sendback.zip'
sys.path.insert(0,str(ROOT/'reference'))
from analyze_cached_v18_9 import get_data as historical_read
import analyze_cached_v18_9 as reader


def sha(data):return hashlib.sha256(data).hexdigest()


def load():
    with zipfile.ZipFile(INPUT) as z:
        corrupt=z.testzip()
        if corrupt:raise ValueError(f'corrupt v18.9 archive member: {corrupt}')
        checks=json.loads(z.read('SHA256SUMS_OUTPUT.json'))['files']
        for n,h in checks.items():
            if sha(z.read(n))!=h:raise ValueError('v18.9 history hash mismatch: '+n)
        inner=z.read('inputs/v18_7_original_sendback.zip')
        artifact=json.loads(z.read('results/summary.json'))
        replay=json.loads(z.read('results/model_replay_summary.json')) if 'results/model_replay_summary.json' in z.namelist() else None
    # The original v18.9 validator expects a path. Keep the nested archive in
    # transient runtime space; never add it to outputs or sendback.
    runtime=ROOT/'.runtime';runtime.mkdir(exist_ok=True)
    nested=runtime/'v18_7_original_sendback.zip';nested.write_bytes(inner)
    reader.ARCHIVE=nested
    data,validated=historical_read()
    return data,{'original_v18_9_sha256':sha(INPUT.read_bytes()),'v18_9_manifest_entries':len(checks),
                 'v18_9_summary_version':artifact.get('version'),
                 'v18_9_ode_replay_success':None if replay is None else replay.get('ode_replay_success'),
                 'historical_validation':validated}
