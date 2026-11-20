"""Require actual nonlinear server proof before packaging certified results."""
from pathlib import Path
import json,zipfile,hashlib
ROOT=Path(__file__).resolve().parents[1]
proof=ROOT/'results/original_ode_robustness_proof.json'
if not proof.exists() or json.loads(proof.read_text()).get('status')!='ORIGINAL_Q550_ENDPOINT_PRESERVED_BLIND_NONLINEAR_VALIDATION_COMPLETED':
    raise RuntimeError('Cannot package frozen preflight as server-certified nonlinear v18.24')
manifest=json.loads((ROOT/'SHA256SUMS_OUTPUT.json').read_text())['files']
for rel,sha in manifest.items():
    if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha:raise RuntimeError('SHA mismatch '+rel)
output=ROOT.parent/(ROOT.name.replace('_bundle','_sendback')+'.zip')
with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as z:
    for rel in sorted(manifest):
        f=ROOT/rel
        z.write(f,arcname=ROOT.name+'/'+rel,compress_type=zipfile.ZIP_STORED if f.suffix=='.zip' else zipfile.ZIP_DEFLATED)
    z.write(ROOT/'SHA256SUMS_OUTPUT.json',arcname=ROOT.name+'/SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(output) as z:
    if z.testzip():raise RuntimeError('ZIP CRC failure')
    if {x.split('/')[0] for x in z.namelist()}!={ROOT.name}:raise RuntimeError('Unexpected ZIP root')
print('SEND_BACK',output)
