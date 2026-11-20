"""Immutable-file SHA256 manifest; never include live logs, caches or transient runtime."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];checks={}
for folder in ('inputs','scripts','reference','tests','results'):
 for p in sorted((ROOT/folder).rglob('*')):
  if not p.is_file() or p.suffix in ('.pyc','.log') or '__pycache__' in p.parts:continue
  checks[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
for name in ('README.md','ENVIRONMENT.json','requirements.txt'):
 checks[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'files':checks},indent=2)+'\n')
print('SHA256',len(checks),'immutable files (active logs excluded)')
