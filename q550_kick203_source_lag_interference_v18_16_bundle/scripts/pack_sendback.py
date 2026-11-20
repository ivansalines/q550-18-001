#!/usr/bin/env python3
"""Package all executables, inputs, products under one exact directory root."""
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
checks=json.loads((root/'SHA256SUMS_OUTPUT.json').read_text())['files']
for rel,dig in checks.items():
    if hashlib.sha256((root/rel).read_bytes()).hexdigest()!=dig:raise RuntimeError('Stale/mismatched SHA256 for '+rel)
out=root.parent/'q550_kick203_source_lag_interference_v18_16_sendback.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
    for p in sorted(root.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(root).as_posix()
        if rel.startswith('.runtime/') or '/__pycache__/' in rel:continue
        z.write(p,arcname=root.name+'/'+rel)
with zipfile.ZipFile(out) as z:
    if z.testzip():raise RuntimeError('Sendback ZIP CRC mismatch')
    roots={p.split('/')[0] for p in z.namelist()};assert roots=={root.name}, roots
print('Sendback verified:',out)
