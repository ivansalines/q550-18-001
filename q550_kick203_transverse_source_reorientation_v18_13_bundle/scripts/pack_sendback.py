#!/usr/bin/env python3
"""Repack full executable genealogy and validated results, hashed."""
import json,hashlib,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT.parent/'q550_kick203_transverse_source_reorientation_v18_13_sendback.zip'
def contents():
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and not any(s in p.parts for s in ('.runtime','__pycache__','.git')) and p.name!='SHA256SUMS_OUTPUT.json' and p.suffix!='.log')
checks={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in contents()}
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'algorithm':'SHA256','files':checks},indent=2)+'\n')
with zipfile.ZipFile(OUT,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in contents()+[ROOT/'SHA256SUMS_OUTPUT.json']:
        z.write(p,str(p.relative_to(ROOT)))
with zipfile.ZipFile(OUT) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(n)).hexdigest()==h for n,h in checks.items())
print('packaged:',OUT, 'verified files:',len(checks))
