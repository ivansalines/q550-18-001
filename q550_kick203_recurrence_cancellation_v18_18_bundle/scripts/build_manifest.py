#!/usr/bin/env python3
"""Hash immutable artifacts; active logs, caches, and manifest excluded."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
checks={}
for p in sorted(root.rglob('*')):
    if not p.is_file():continue
    rel=p.relative_to(root).as_posix()
    if (rel=='SHA256SUMS_OUTPUT.json' or rel.startswith('.runtime/') or
        '/__pycache__/' in rel or rel.startswith('__pycache__/') or
        rel.endswith('.log') or rel.endswith('.tmp') or rel.startswith('logs/')):continue
    checks[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
(root/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'files':checks},indent=2)+'\n')
print('SHA256',len(checks),'immutable files (live .log excluded)')
