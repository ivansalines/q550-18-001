"""Complete reproducibility sendback with SHA256 manifest and original Q550 ancestry."""
from pathlib import Path
import zipfile,hashlib,json
root=Path(__file__).resolve().parents[1]
name='q550_kick203_material_metric_constraint_sensitivity_v18_8_sendback.zip'
path=root/name
ignore={name,'SHA256SUMS.json','v18_8_server.log'}
files=sorted(p for p in root.rglob('*') if p.is_file() and p.name not in ignore and '__pycache__' not in p.parts and not p.name.endswith('.pyc'))
man={'version':'2026-10-17-q550-v18.8','files':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(root/'SHA256SUMS.json').write_text(json.dumps(man,indent=2)+'\n')
with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in files+[root/'SHA256SUMS.json']:z.write(p,p.relative_to(root))
with zipfile.ZipFile(path) as z:
    damaged=z.testzip()
    if damaged:raise RuntimeError('bad zip: '+damaged)
    for k,h in man['files'].items():
        if hashlib.sha256(z.read(k)).hexdigest()!=h:raise RuntimeError('bad sha '+k)
print(f'PACKED {path} ({len(files)+1} entries, verified)')
