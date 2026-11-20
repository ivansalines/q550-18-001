#!/usr/bin/env python3
"""Package only files in the certified immutable manifest, under one root."""
from pathlib import Path
import hashlib,json,zipfile
root=Path(__file__).resolve().parents[1]
meta=(root/'SHA256SUMS_OUTPUT.json').read_bytes()
checks=json.loads(meta)['files']
for rel,digest in checks.items():
    if hashlib.sha256((root/rel).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale immutable source: '+rel)
out=root.parent/'q550_kick203_moving_antagonism_ridge_v18_17_sendback.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
    for rel in sorted(checks):z.write(root/rel,arcname=root.name+'/'+rel)
    z.writestr(root.name+'/SHA256SUMS_OUTPUT.json',meta)
with zipfile.ZipFile(out) as z:
    if z.testzip():raise RuntimeError('CRC failure')
    roots={name.split('/')[0] for name in z.namelist()}
    if roots!={root.name}:raise RuntimeError('Unexpected ZIP root')
    for rel,dig in checks.items():
        if hashlib.sha256(z.read(root.name+'/'+rel)).hexdigest()!=dig:raise RuntimeError('ZIP digest mismatch: '+rel)
print('Certified output:',out)
