"""Hash only immutable files; never hash a tee-managed active log."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
files={}
for folder in ('inputs','scripts','reference','tests','results'):
    for f in sorted((ROOT/folder).rglob('*')):
        if not f.is_file() or f.suffix in ('.pyc','.log') or '__pycache__' in f.parts:continue
        files[f.relative_to(ROOT).as_posix()]=hashlib.sha256(f.read_bytes()).hexdigest()
for name in ('README.md','requirements.txt','ENVIRONMENT.json'):
    f=ROOT/name;files[name]=hashlib.sha256(f.read_bytes()).hexdigest()
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'files':files},indent=2)+'\n')
print('SHA256',len(files),'immutable files, active logs excluded')
