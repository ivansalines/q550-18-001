"""Create verified server sendback without hashing running tee logs."""
from pathlib import Path
import hashlib,json,zipfile
ROOT=Path(__file__).resolve().parents[1]
if not (ROOT/'results/nonlinear_jet_summary.json').is_file():
 raise RuntimeError('NO ORIGINAL ODE NONLINEAR RESULT: cannot package frozen preflight as server result')
if json.loads((ROOT/'results/nonlinear_jet_summary.json').read_text()).get('status')!='ORIGINAL_Q550_NONLINEAR_JET_VALIDATION_COMPLETED':
 raise RuntimeError('Not a completed nonlinear experiment')
checks=json.loads((ROOT/'SHA256SUMS_OUTPUT.json').read_text())['files']
for rel,sha in checks.items():
 if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha:raise RuntimeError('SHA256 mismatch '+rel)
out=ROOT.parent/(ROOT.name.replace('_bundle','_sendback')+'.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for rel in sorted(checks):z.write(ROOT/rel,arcname=ROOT.name+'/'+rel)
 z.write(ROOT/'SHA256SUMS_OUTPUT.json',arcname=ROOT.name+'/SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(out) as z:
 if z.testzip():raise RuntimeError('sendback ZIP CRC failure')
print('SEND_BACK',out)
