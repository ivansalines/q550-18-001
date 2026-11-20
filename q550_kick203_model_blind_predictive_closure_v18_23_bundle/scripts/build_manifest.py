"""Hash immutable files only; exclude live logs, transient .runtime, and self-reference."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
files={}
for p in ROOT.rglob('*'):
    if not p.is_file():continue
    rel=p.relative_to(ROOT).as_posix()
    if rel=='SHA256SUMS_OUTPUT.json' or rel.startswith('.runtime/') or p.suffix=='.pyc' or '/__pycache__/' in rel or rel.endswith('.log') or rel.endswith('.tmp'):
        continue
    files[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'algorithm':'SHA256','files':dict(sorted(files.items()))},indent=2)+'\n')
print('SHA256',len(files),'immutable files')
