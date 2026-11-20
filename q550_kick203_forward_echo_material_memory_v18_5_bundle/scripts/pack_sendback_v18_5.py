#!/usr/bin/env python3
"""Self-contained ZIP with runnable code, original lineage, all results and hashes."""
from pathlib import Path
import hashlib,json,zipfile,sys
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'q550_kick203_forward_echo_material_memory_v18_5_sendback.zip'
required=[ROOT/'results'/'summary.json',ROOT/'results'/'forward_echo_samples.csv',ROOT/'checkpoint'/'run_provenance.json']
required += [ROOT/'checkpoint'/f'{tol}_eps{eps:g}_delay{lag:g}.npz'
    for tol in ('standard','tight') for eps in (1e-8,3e-8,1e-7) for lag in (.25,.5,1,2,4,8)]
missing=[str(p) for p in required if not p.exists()]
if missing:sys.exit('Refuse incomplete sendback: '+str(missing[:6]))
files=[]
for folder in ('inputs','reference','scripts','tests','checkpoint','results'):
    files += [p for p in sorted((ROOT/folder).rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc']
for n in ('README.md','requirements.txt','SOURCE_MANIFEST.json','v18_5_solver.log','v18_5_analysis.log','v18_5_server.log'):
    if (ROOT/n).exists():files.append(ROOT/n)
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(manifest,indent=2)+'\n')
files.append(ROOT/'SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(OUT,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,arcname=str(p.relative_to(ROOT)))
with zipfile.ZipFile(OUT) as z: assert z.testzip() is None
print(f'Packed {OUT}, {len(files)} files, {OUT.stat().st_size} bytes')
