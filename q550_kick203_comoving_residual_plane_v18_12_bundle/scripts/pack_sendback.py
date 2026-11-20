#!/usr/bin/env python3
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1];out=root.parent/'q550_kick203_comoving_residual_plane_v18_12_sendback.zip'
files=sorted(p for p in root.rglob('*') if p.is_file() and not '.runtime' in p.parts and not '__pycache__' in p.parts)
checks={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.name!='SHA256SUMS_OUTPUT.json'}
(root/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'algorithm':'SHA256','files':checks},indent=2)+'\n')
files=sorted(p for p in root.rglob('*') if p.is_file() and '.runtime' not in p.parts and '__pycache__' not in p.parts)
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files:z.write(f,str(f.relative_to(root)))
print('sendback:',out)
