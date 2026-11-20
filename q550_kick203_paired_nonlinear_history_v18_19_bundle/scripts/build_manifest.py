"""Hash immutable inputs/code/finished results, never active *.log or .runtime."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
files={}
for sub in ('inputs','scripts','reference','tests','results'):
 for path in sorted((ROOT/sub).rglob('*')):
  if not path.is_file() or path.suffix in ('.pyc','.log') or '__pycache__' in path.parts:continue
  files[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
for rel in ('README.md','ENVIRONMENT.json','requirements.txt'):
 files[rel]=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'files':files},indent=2)+'\n')
print('SHA256',len(files),'immutable files (active logs excluded)')
