#!/usr/bin/env python3
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root.parent/'q550_kick203_phase_zero_independent_check_v18_11_1_sendback.zip'
files=sorted(p for p in root.rglob('*') if p.is_file() and not ({'__pycache__','.runtime','.pytest_cache'} & set(p.relative_to(root).parts)) and p.name!='SHA256SUMS_OUTPUT.json' and p.suffix!='.pyc')
manifest={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(root/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'version':'v18.11.1','files':manifest},indent=2)+'\n')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files+[root/'SHA256SUMS_OUTPUT.json']:z.write(p,p.relative_to(root))
with zipfile.ZipFile(out) as z:
    if z.testzip() is not None:raise RuntimeError('corrupt sendback zip')
print('SEND_BACK_ZIP='+str(out))
