#!/usr/bin/env python3
"""Bundle the runnable source, original archive, numerical checkpoints and reports."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,zipfile
ROOT=Path(__file__).resolve().parents[1]
name='q550_kick203_second_order_material_source_transport_v18_6_sendback.zip'
archive=ROOT/name
paths=[]
for folder in ['scripts','tests','reference','inputs','checkpoint','results']:
    paths.extend([x for x in (ROOT/folder).rglob('*') if x.is_file() and '__pycache__' not in x.parts and not x.name.endswith('.pyc')])
for f in ['README.md','requirements.txt','SOURCE_MANIFEST.json','v18_6_server.log','v18_6_solver.log','v18_6_analysis.log']:
    p=ROOT/f
    if p.exists():paths.append(p)
paths=sorted(set(paths)); hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(dict(archive_type='v18.6 Q550 source audit',utc=datetime.now(timezone.utc).isoformat(),files=hashes),indent=2))
paths.append(ROOT/'SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in paths:z.write(p,p.relative_to(ROOT).as_posix())
with zipfile.ZipFile(archive) as z:
    if z.testzip():raise RuntimeError('Bad ZIP CRC')
print('CREATED',archive,'bytes',archive.stat().st_size,'members',len(paths))
