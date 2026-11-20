#!/usr/bin/env python3
"""Pack the full executable genealogy and measured server outputs."""
import argparse,hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/'q550_kick203_on_shell_full11d_material_transport_v18_2_sendback.zip')
    args=ap.parse_args(); target=args.out.resolve();target.parent.mkdir(parents=True,exist_ok=True)
    files=[f for f in ROOT.rglob('*') if f.is_file() and f.resolve()!=target and
           '__pycache__' not in f.parts and f.suffix not in ('.pyc',) and
           not f.name.endswith('.zip.part')]
    result=ROOT/'results/kick203_v18_2/summary.json'
    if not result.exists():raise RuntimeError('no completed full11D material audit; refuse to package empty sendback')
    catalog={f.relative_to(ROOT).as_posix():{'sha256':sha(f),'size':f.stat().st_size} for f in files}
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for f in files:z.write(f,f.relative_to(ROOT).as_posix())
        z.writestr('SHA256SUMS_SENDBACK.json',json.dumps(catalog,indent=2)+'\n')
    print(f'packed {len(files)} files -> {target} ({target.stat().st_size:,} bytes)',flush=True)
