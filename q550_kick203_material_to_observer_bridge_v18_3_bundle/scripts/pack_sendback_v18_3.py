#!/usr/bin/env python3
"""Create a provenance-complete v18.3 return archive without re-running numerical data."""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED,ZIP_STORED
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
def main():
    outfile=ROOT/'q550_kick203_material_to_observer_bridge_v18_3_sendback.zip'
    files=[p for p in ROOT.rglob('*') if p.is_file() and p!=outfile and p.name !='SHA256SUMS.json' and '__pycache__' not in str(p)]
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (ROOT/'SHA256SUMS.json').write_text(json.dumps(hashes,indent=2)+'\n')
    with ZipFile(outfile,'w',ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(files+[ROOT/'SHA256SUMS.json']):
            name=str(p.relative_to(ROOT))
            z.write(p,name,compress_type=ZIP_STORED if name.endswith('.zip') or name.endswith('.npz') else ZIP_DEFLATED)
    with ZipFile(outfile) as z:
        bad=z.testzip()
        if bad:raise RuntimeError(f'invalid archive entry {bad}')
    print('Sendback archive:',outfile,'files:',len(files)+1)
if __name__=='__main__':main()
