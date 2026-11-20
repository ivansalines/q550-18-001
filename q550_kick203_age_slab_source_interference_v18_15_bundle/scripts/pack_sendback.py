#!/usr/bin/env python3
"""ZIP executable genealogy with reproducible relative-root layout and hashes."""
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT.parent/'q550_kick203_age_slab_source_interference_v18_15_sendback.zip'

def keep(p):
    return p.is_file() and not any(s in ('.runtime','__pycache__','.git') for s in p.parts) and p.name not in ('SHA256SUMS_OUTPUT.json','v18_15_server.log')
files=sorted((p for p in ROOT.rglob('*') if keep(p)),key=lambda p:str(p))
checks={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'algorithm':'SHA256','files':checks},indent=2)+'\n')
with zipfile.ZipFile(OUT,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files+[ROOT/'SHA256SUMS_OUTPUT.json']:z.write(p,arcname=ROOT.name+'/'+str(p.relative_to(ROOT)))
with zipfile.ZipFile(OUT) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(ROOT.name+'/'+name)).hexdigest()==sha for name,sha in checks.items())
print('Verified sendback:',OUT,'members:',len(checks))
