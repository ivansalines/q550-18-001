#!/usr/bin/env python3
"""Pack independent runnable genealogy plus outputs into a hashed sendback ZIP."""
from __future__ import annotations
import datetime,hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
NAME='q550_kick203_source_birth_sign_crossing_v18_7_sendback.zip'

def sha(file):return hashlib.sha256(file.read_bytes()).hexdigest()

def main():
    expected=[ROOT/'results/summary.json',ROOT/'results/REPORT.md',
              ROOT/'checkpoint/standard_sign_crossing.npz',ROOT/'checkpoint/tight_sign_crossing.npz']
    if not all(f.is_file() for f in expected):
        raise SystemExit('Cannot pack; both actual-Q550 tolerance checkpoints and summary must exist')
    entries=[]
    for d in ['inputs','reference','scripts','tests','checkpoint','results']:
        entries.extend(f for f in (ROOT/d).rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.suffix!='.pyc')
    entries += [ROOT/'README.md',ROOT/'requirements.txt',ROOT/'SOURCE_MANIFEST.json',ROOT/'SHA256SUMS.json']
    entries += [f for f in (ROOT/'v18_7_solver.log',ROOT/'v18_7_analysis.log',ROOT/'v18_7_server.log') if f.exists()]
    entries=sorted(set(entries),key=lambda f:f.relative_to(ROOT).as_posix())
    man={'archive_type':'Q550 v18.7 source birth crossing sendback',
         'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'files':{f.relative_to(ROOT).as_posix():sha(f) for f in entries}}
    (ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(man,indent=2)+'\n')
    out=ROOT/NAME
    with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for f in entries:z.write(f,f.relative_to(ROOT).as_posix())
        z.write(ROOT/'SHA256SUMS_OUTPUT.json','SHA256SUMS_OUTPUT.json')
    with zipfile.ZipFile(out) as z:
        if z.testzip() is not None:raise RuntimeError('Failed sendback ZIP integrity check')
    print(out,'members:',len(entries)+1)
if __name__=='__main__':main()
