#!/usr/bin/env python3
"""Manifest all distributed files without self-reference or transient caches."""
import json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
files={}
for p in sorted(root.rglob('*')):
    if not p.is_file():continue
    rel=p.relative_to(root).as_posix()
    if rel.startswith(('.runtime/','__pycache__/')) or '/__pycache__/' in rel or rel=='SHA256SUMS_OUTPUT.json':continue
    files[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
(root/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'files':files},indent=2)+'\n')
print('SHA256',len(files),'files')
