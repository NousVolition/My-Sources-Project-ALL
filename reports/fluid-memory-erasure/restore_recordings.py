"""Verify and restore the released raw recordings; no network by default."""
import argparse,hashlib,json,tarfile
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for value in iter(lambda:f.read(8*1024**2),b''):h.update(value)
    return h.hexdigest()


def main():
    root=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser();ap.add_argument('--parts',type=Path,default=root);ap.add_argument('--out',type=Path,default=root);ap.add_argument('--download',action='store_true');args=ap.parse_args()
    manifest=json.loads((root/'raw-recordings.json').read_text())
    args.parts.mkdir(parents=True,exist_ok=True);args.out.mkdir(parents=True,exist_ok=True)
    if args.download:
        import urllib.request
        base=f'https://github.com/{manifest["repository"]}/releases/download/{manifest["release_tag"]}/'
        for part in manifest['parts']:
            path=args.parts/part['file']
            if not path.exists():urllib.request.urlretrieve(base+part['file'],path)
    combined=args.parts/manifest.get('archive_name','fluid-memory-recordings.tar.xz')
    if not combined.exists() or sha(combined)!=manifest['archive_sha256']:
        with combined.open('wb') as out:
            for part in manifest['parts']:
                path=args.parts/part['file']
                if sha(path)!=part['sha256']:raise RuntimeError('Part checksum failed: '+str(path))
                with path.open('rb') as src:
                    for value in iter(lambda:src.read(8*1024**2),b''):out.write(value)
    if sha(combined)!=manifest['archive_sha256']:raise RuntimeError('Joined checksum failed')
    target=args.out.resolve()
    with tarfile.open(combined,'r:*') as archive:
        for member in archive:
            path=(target/member.name).resolve()
            if not path.is_relative_to(target) or not member.isfile():raise RuntimeError('Unsafe archive member')
            if path.exists():
                existing=path.read_bytes();proposed=archive.extractfile(member).read()
                if existing!=proposed:raise RuntimeError('Refusing to replace differing existing file: '+str(path))
                continue
            path.parent.mkdir(parents=True,exist_ok=True)
            with archive.extractfile(member) as source,path.open('wb') as out:
                for value in iter(lambda:source.read(8*1024**2),b''):out.write(value)
    for record in manifest['records']:
        if sha(target/record['file'])!=record['sha256']:raise RuntimeError('Recording checksum failed')
    print(f'Restored and verified {len(manifest["records"])} exact-array recordings. No simulations were rerun.')

if __name__=='__main__':main()
