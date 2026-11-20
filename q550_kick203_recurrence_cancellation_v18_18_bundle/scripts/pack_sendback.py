#!/usr/bin/env python3
"""Check immutable SHA256 before packing, certify root/CRC/checks after."""
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
meta=(root/'SHA256SUMS_OUTPUT.json').read_bytes()
checks=json.loads(meta)['files']
for rel,dig in checks.items():
    if hashlib.sha256((root/rel).read_bytes()).hexdigest()!=dig:
        raise RuntimeError('Stale immutable file: '+rel)
out=root.parent/'q550_kick203_recurrence_cancellation_v18_18_sendback.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
    for rel in sorted(checks):z.write(root/rel,arcname=root.name+'/'+rel)
    z.writestr(root.name+'/SHA256SUMS_OUTPUT.json',meta)
with zipfile.ZipFile(out) as z:
    bad=z.testzip()
    if bad:raise RuntimeError('Sendback CRC failure '+bad)
    roots={s.split('/')[0] for s in z.namelist()}
    if roots!={root.name}:raise RuntimeError('Invalid sendback root')
    for rel,dig in checks.items():
        if hashlib.sha256(z.read(root.name+'/'+rel)).hexdigest()!=dig:
            raise RuntimeError('Sendback hash mismatch '+rel)
print('Certified output:',out)
