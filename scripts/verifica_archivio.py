"""Verifica con la libreria standard l'archivio e i suoi file estratti."""
from pathlib import Path
import hashlib,json,zipfile

ROOT=Path(__file__).resolve().parent.parent
archive=ROOT/'outputs/grafene-modelli-e-sorgenti.zip'
manifest=json.loads((ROOT/'MANIFEST-ARCHIVIO.json').read_text())

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while block:=f.read(1024*1024):h.update(block)
    return h.hexdigest()

with zipfile.ZipFile(archive) as z:
    stored=json.loads(z.read('MANIFEST.json'))
    assert stored==manifest
    assert z.testzip() is None
    for row in manifest['files']:
        label=Path(row['path']);assert not label.is_absolute() and '..' not in label.parts
        if label.parts[0]=='paper':path=ROOT/'outputs'/Path(*label.parts[1:])
        elif label.parts[0]=='molecular':path=ROOT/'work/molecular'/Path(*label.parts[1:])
        elif label.parts[0]=='theory':path=ROOT/'work'/Path(*label.parts[1:])
        else:raise RuntimeError(label)
        assert digest(path)==row['sha256'],path
        assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256'],row['path']
print(f"Verificati {len(manifest['files'])} file; revisione {manifest['revision']}.")
