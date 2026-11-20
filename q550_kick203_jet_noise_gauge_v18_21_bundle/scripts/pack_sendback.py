"""Pack only after confirmed original-ODE validation; no inclusion of live logs."""
from pathlib import Path
import json,hashlib,zipfile
ROOT=Path(__file__).resolve().parents[1]
proof=ROOT/'results/original_ode_check_summary.json'
if not proof.exists() or json.loads(proof.read_text()).get('status')!='ORIGINAL_Q550_INDEPENDENT_NONLINEAR_SPOT_CHECK_COMPLETED':
    raise RuntimeError('Original Q550 nonlinear verification not completed; frozen-only archive cannot be called server sendback')
checks=json.loads((ROOT/'SHA256SUMS_OUTPUT.json').read_text())['files']
for rel,digest in checks.items():
    if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=digest:
        raise RuntimeError('SHA256 mismatch '+rel)
out=ROOT.parent/(ROOT.name.replace('_bundle','_sendback')+'.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for rel in sorted(checks):
        f=ROOT/rel
        # Pre-compressed original original-server genealogy can be safely stored, not recompressed.
        if f.suffix=='.zip':z.write(f,arcname=ROOT.name+'/'+rel,compress_type=zipfile.ZIP_STORED)
        else:z.write(f,arcname=ROOT.name+'/'+rel)
    z.write(ROOT/'SHA256SUMS_OUTPUT.json',arcname=ROOT.name+'/SHA256SUMS_OUTPUT.json')
with zipfile.ZipFile(out) as z:
    if z.testzip():raise RuntimeError('CRC failure in sendback')
    roots={n.split('/')[0] for n in z.namelist()}
    if roots!={ROOT.name}:raise RuntimeError('Wrong sendback ZIP hierarchy')
print('SEND_BACK',out)
