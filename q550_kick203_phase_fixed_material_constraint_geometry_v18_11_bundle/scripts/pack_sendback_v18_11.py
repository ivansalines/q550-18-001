#!/usr/bin/env python3
"""Write hash-validated reproducible end-to-end Q550 v18.11 archive."""
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SEND=ROOT.parent/'q550_kick203_phase_fixed_material_constraint_geometry_v18_11_sendback.zip'
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and not ({'.runtime','__pycache__','.pytest_cache'}&set(p.relative_to(ROOT).parts)) and p.name!='SHA256SUMS_OUTPUT.json' and p.suffix!='.pyc')
sums={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'version':'v18.11','files':sums},indent=2)+'\n')
files.append(ROOT/'SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(SEND,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,p.relative_to(ROOT))
with zipfile.ZipFile(SEND) as z:
 if z.testzip():raise RuntimeError('sendback CRC failed')
print('SEND_BACK_ZIP='+str(SEND))
