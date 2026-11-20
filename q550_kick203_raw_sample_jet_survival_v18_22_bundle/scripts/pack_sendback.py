"""Refuse server-certified archive unless fresh original-Q550 nonlinear proof exists."""
from pathlib import Path
import json,hashlib,zipfile
ROOT=Path(__file__).resolve().parents[1]
proof=ROOT/'results/original_ode_raw_summary.json'
if not proof.exists() or json.loads(proof.read_text()).get('status')!='ORIGINAL_Q550_RAW_SAMPLE_NONLINEAR_VALIDATION_COMPLETED':
    raise RuntimeError('Cannot call a frozen-only result original-server-certified')
manifest=json.loads((ROOT/'SHA256SUMS_OUTPUT.json').read_text())['files']
for rel,digest in manifest.items():
    if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=digest:
        raise RuntimeError('SHA mismatch '+rel)
out=ROOT.parent/(ROOT.name.replace('_bundle','_sendback')+'.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for rel in sorted(manifest):
        f=ROOT/rel
        z.write(f,arcname=ROOT.name+'/'+rel,compress_type=zipfile.ZIP_STORED if f.suffix=='.zip' else zipfile.ZIP_DEFLATED)
    z.write(ROOT/'SHA256SUMS_OUTPUT.json',arcname=ROOT.name+'/SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(out) as z:
    if z.testzip():raise RuntimeError('CRC failed')
    if {n.split('/')[0] for n in z.namelist()}!={ROOT.name}:raise RuntimeError('Unexpected root')
print('SEND_BACK',out)
