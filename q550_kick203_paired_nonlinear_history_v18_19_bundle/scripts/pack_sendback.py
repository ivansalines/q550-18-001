"""Pack verified immutable server outputs; never hash active console log."""
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
report=ROOT/'results/nonlinear_pair_summary.json'
if not report.is_file():raise RuntimeError('SERVER NONLINEAR AUDIT NOT PRESENT; refusing to package preflight as server result')
checks=json.loads((ROOT/'SHA256SUMS_OUTPUT.json').read_text())['files']
for name,expected in checks.items():
 if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:raise RuntimeError('SHA mismatch '+name)
items=list(checks)+['SHA256SUMS_OUTPUT.json']
parent=ROOT.parent
out=parent/(ROOT.name.replace('_bundle','_sendback')+'.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for rel in items:z.write(ROOT/rel,arcname=ROOT.name+'/'+rel)
with zipfile.ZipFile(out) as z:
 if z.testzip():raise RuntimeError('Sendback CRC failed')
print('SEND_BACK',out)
