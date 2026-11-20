#!/usr/bin/env python3
"""Package executable full genealogy, analytical outputs and SHA256 manifest."""
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT.parent/'q550_kick203_history_window_echo_v18_14_sendback.zip'
def keep(p):
    return p.is_file() and not any(a in p.parts for a in ('.runtime','__pycache__','.git')) and p.name!='SHA256SUMS_OUTPUT.json' and p.suffix!='.log'
files=sorted(p for p in ROOT.rglob('*') if keep(p))
checks={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'algorithm':'SHA256','files':checks},indent=2)+'\n')
with zipfile.ZipFile(OUT,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files+[ROOT/'SHA256SUMS_OUTPUT.json']:z.write(p,str(p.relative_to(ROOT)))
with zipfile.ZipFile(OUT) as z:
    failed=z.testzip()
    if failed:raise RuntimeError('ZIP CRC failure: '+failed)
    assert all(hashlib.sha256(z.read(n)).hexdigest()==h for n,h in checks.items())
print('Sendback:',OUT,'verified member hashes:',len(checks))
