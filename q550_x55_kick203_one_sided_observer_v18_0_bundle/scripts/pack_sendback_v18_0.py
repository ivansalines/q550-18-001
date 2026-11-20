#!/usr/bin/env python3
"""Pack all new measurements plus the actual source scripts and dependency manifest."""
import argparse,zipfile,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
root=Path(__file__).resolve().parents[1]
files=list(args.results.glob('*'))+list((root/'scripts').glob('*'))+list((root/'reference').glob('*.py'))+list((root/'inputs').glob('*'))+[root/'requirements.txt',root/'README.md',root/'SOURCE_MANIFEST.json',root/'SHA256SUMS.json']
files=[f for f in files if f.is_file() and f.suffix!='.pyc']
hashes={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
with zipfile.ZipFile(args.out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files:z.write(f,str(f.relative_to(root)))
 z.writestr('SHA256SUMS_RESULTS.json',json.dumps(hashes,indent=2)+'\n')
print('sendback:',args.out)
