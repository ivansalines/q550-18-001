#!/usr/bin/env python3
"""Pack a self-contained v18.9 source+evidence ZIP and new SHA256 manifest."""
from __future__ import annotations
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT.parent/'q550_kick203_ode_intrinsic_metric_constraint_v18_9_sendback.zip'
EXCLUDED={'__pycache__','.pytest_cache'}
FILES=[p for p in sorted(ROOT.rglob('*')) if p.is_file() and not (set(p.parts)&EXCLUDED)
       and p.name!='SHA256SUMS_OUTPUT.json' and not p.name.endswith('.pyc')]
SUMS={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in FILES}
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'version':'v18.9','files':SUMS},indent=2)+'\n')
FILES.append(ROOT/'SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(OUTPUT,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in FILES:z.write(p,p.relative_to(ROOT))
print(OUTPUT)
