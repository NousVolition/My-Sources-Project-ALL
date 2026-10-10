"""Verify and restore recorded arrays without overwriting different files."""
from pathlib import Path
import hashlib,io,json,zipfile
ROOT=Path(__file__).resolve().parent
def digest(b):return hashlib.sha256(b).hexdigest()
def main():
    index=json.loads((ROOT/'recordings.json').read_text(encoding='utf-8'));pieces=[]
    for item in index['parts']:
        b=(ROOT/item['path']).read_bytes()
        if len(b)!=item['size'] or digest(b)!=item['sha256']:raise ValueError('Corrupt part '+item['path'])
        pieces.append(b)
    b=b''.join(pieces)
    if len(b)!=index['archive_size'] or digest(b)!=index['archive_sha256']:raise ValueError('Corrupt archive')
    expected={e['path']:e for e in index['files']}
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        if set(z.namelist())!=set(expected):raise ValueError('Archive file list differs')
        for name in z.namelist():
            path=(ROOT/name).resolve()
            if not path.is_relative_to(ROOT.resolve()):raise ValueError('Unsafe archive path')
            raw=z.read(name)
            if len(raw)!=expected[name]['size'] or digest(raw)!=expected[name]['sha256']:raise ValueError('Corrupt entry '+name)
            if path.exists() and digest(path.read_bytes())!=expected[name]['sha256']:raise ValueError('Different existing file '+name)
        for name in z.namelist():
            path=ROOT/name;path.parent.mkdir(parents=True,exist_ok=True)
            if not path.exists():path.write_bytes(z.read(name))
    print(f"Restored and verified {len(expected)} recorded files")
if __name__=='__main__':main()
