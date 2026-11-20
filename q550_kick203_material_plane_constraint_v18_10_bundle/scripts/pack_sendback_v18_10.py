#!/usr/bin/env python3
"""Pack scripts, actual inputs, reports and hashes as an executable genealogy."""
from __future__ import annotations
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT.parent/'q550_kick203_material_plane_constraint_v18_10_sendback.zip'
files=[p for p in ROOT.rglob('*') if p.is_file() and not (set(p.relative_to(ROOT).parts)&{'.runtime','__pycache__','.pytest_cache'})
       and p.suffix!='.pyc' and p.name!='SHA256SUMS_OUTPUT.json']
files=sorted(files)
sums={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'version':'v18.10','files':sums},indent=2)+'\n')
files.append(ROOT/'SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(DEST,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,p.relative_to(ROOT))
with zipfile.ZipFile(DEST) as z:
 if z.testzip():raise RuntimeError('sendback zip CRC failure')
print('SEND_BACK_ZIP='+str(DEST))
