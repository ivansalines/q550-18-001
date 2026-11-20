#!/usr/bin/env python3
"""Self-contained sendback with source scripts, reference, checkpoints, and outputs."""
import hashlib,json,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TARGET=ROOT/'q550_kick203_five_amplitude_scaling_v18_4_sendback.zip'
required=['results/summary.json','results/five_amplitude_samples.csv',
          'results/samplewise_odd_fits.csv']
required += [f'checkpoint/{tol}_eps{eps:g}.npz' for tol in ('standard','tight') for eps in (5e-9,1.5e-8,2e-8)]
missing=[p for p in required if not (ROOT/p).exists()]
if missing:sys.exit('Refusing partial sendback; missing '+str(missing))
manifest={}
files=[]
for subtree in ('scripts','reference','inputs','tests','checkpoint','results'):
    for p in sorted((ROOT/subtree).rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and not p.name.endswith('.pyc'):
            files.append(p)
for name in ('README.md','requirements.txt','SOURCE_MANIFEST.json','v18_4_solver.log','v18_4_analysis.log','v18_4_server.log'):
    if (ROOT/name).is_file():files.append(ROOT/name)
for p in files:
    manifest[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(manifest,indent=2)+'\n')
files.append(ROOT/'SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(TARGET,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,arcname=str(p.relative_to(ROOT)))
with zipfile.ZipFile(TARGET) as z:
    assert z.testzip() is None
print('Packed',TARGET,'files=',len(files),'size=',TARGET.stat().st_size)
