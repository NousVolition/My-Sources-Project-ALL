"""Restore the original completed pilot from locally included archive parts.

Uses only Python's standard library. Existing files must match the archive;
different files are preserved and restoration stops before writing anything.
"""
from pathlib import Path,PurePosixPath
import hashlib,json,tempfile,zipfile

BASE=Path(__file__).resolve().parent

def restore(target=BASE):
    target=Path(target).resolve()
    index=json.loads((BASE/'recorded-package.json').read_text())
    total=hashlib.sha256();size=0
    with tempfile.TemporaryFile() as tmp:
        for part in index['parts']:
            data=(BASE/part['path']).read_bytes()
            if len(data)!=part['size'] or hashlib.sha256(data).hexdigest()!=part['sha256']:
                raise ValueError('Archive part does not match: '+part['path'])
            total.update(data);size+=len(data);tmp.write(data)
        if size!=index['size'] or total.hexdigest()!=index['sha256']:
            raise ValueError('Complete archive does not match its checksum')
        tmp.seek(0)
        with zipfile.ZipFile(tmp) as archive:
            members=[]
            for info in archive.infolist():
                rel=PurePosixPath(info.filename)
                if rel.is_absolute() or '..' in rel.parts or not rel.parts or rel.parts[0]!='fluid_pilot':
                    raise ValueError('Unexpected archive path: '+info.filename)
                if info.is_dir():continue
                path=target.joinpath(*rel.parts[1:]).resolve()
                if not path.is_relative_to(target):raise ValueError('Path escapes destination')
                data=archive.read(info)
                if path.exists() and path.read_bytes()!=data:
                    raise FileExistsError('Preserved different existing file: '+str(path))
                members.append((path,data))
            manifest=json.loads(archive.read('fluid_pilot/manifest.json'))['files']
            for path,data in members:
                rel=path.relative_to(target).as_posix()
                if rel in manifest and hashlib.sha256(data).hexdigest()!=manifest[rel]:
                    raise ValueError('Original file checksum mismatch: '+rel)
            for path,data in members:
                if not path.exists():
                    path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    print(f'Verified and restored {len(members)} original files to {target}')

if __name__=='__main__':restore()
