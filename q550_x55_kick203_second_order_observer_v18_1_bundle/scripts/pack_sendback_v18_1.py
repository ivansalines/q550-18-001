#!/usr/bin/env python3
"""Package already computed v18.1 results; do not redo the orbit/integration."""
from __future__ import annotations
import argparse, hashlib, json, zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--results',type=Path,default=ROOT/'results/kick203_v18_1')
p.add_argument('--out',type=Path,default=ROOT/'q550_kick203_second_order_observer_v18_1_sendback.zip')
args=p.parse_args()
root=ROOT.resolve()
results=args.results.resolve()
out=args.out.resolve()
if not results.is_dir(): raise SystemExit(f'Missing results folder: {results}')
if not (results/'summary.json').is_file(): raise SystemExit('No summary.json in results folder: did the generator finish?')
files=(list(results.rglob('*'))+list((root/'scripts').glob('*'))+
       list((root/'tests').glob('*.py'))+list((root/'reference').glob('*.py'))+
       list((root/'inputs').glob('*')) +
       [root/n for n in ('README.md','SOURCE_MANIFEST.json','SHA256SUMS.json','requirements.txt')])
files=sorted({f.resolve() for f in files if f.is_file() and f.suffix!='.pyc'})
for f in files:
    if not f.is_relative_to(root): raise RuntimeError(f'Unpackaged path outside bundle: {f}')
checks={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
out.parent.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
    for f in files: archive.write(f,str(f.relative_to(root)))
    archive.writestr('SHA256SUMS_RESULTS.json',json.dumps(checks,indent=2)+'\n')
with zipfile.ZipFile(out) as z:
    problem=z.testzip()
    if problem: raise RuntimeError(f'ZIP corruption: {problem}')
print('SEND_BACK_READY =',out)
print('SHA256 manifest items =',len(checks))
