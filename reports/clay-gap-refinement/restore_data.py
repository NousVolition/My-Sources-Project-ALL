"""Restore exact release arrays, refusing traversal and conflicting files."""
from pathlib import Path,PurePosixPath
import argparse,json,hashlib,urllib.request,zipfile,tempfile
R=Path(__file__).resolve().parent
def digest(p):
 with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path);ap.add_argument('--destination',type=Path,default=R);ap.add_argument('--verify-only',action='store_true');a=ap.parse_args();catalog=json.loads((R/'data-release.json').read_text());dest=a.destination.resolve();dest.mkdir(parents=True,exist_ok=True);rows=catalog['files'];names=[r['member'] for r in rows];assert len(set(names))==len(rows)
 for r in rows:
  rel=PurePosixPath(r['path']);assert not rel.is_absolute() and '..' not in rel.parts and '\\' not in r['path'];target=(dest/r['path']).resolve();assert target.is_relative_to(dest)
  if target.exists() and (target.stat().st_size!=r['bytes'] or digest(target)!=r['sha256']):raise ValueError('Refusing to overwrite a conflicting file: '+str(target))
 missing=[r for r in rows if not (dest/r['path']).exists()]
 if missing and a.verify_only:raise FileNotFoundError(f'{len(missing)} arrays missing; run restore_data.py first')
 if missing:
  with tempfile.TemporaryDirectory(dir=dest,prefix='hug-data-download-') as folder:
   assert Path(folder).resolve().is_relative_to(dest)
   archive=a.archive.resolve() if a.archive else Path(folder)/catalog['archive']['name']
   if not a.archive:urllib.request.urlretrieve(catalog['archive']['url'],archive)
   assert archive.stat().st_size==catalog['archive']['bytes'] and digest(archive)==catalog['archive']['sha256'],'Archive hash/size mismatch'
   with zipfile.ZipFile(archive) as z:
    assert z.namelist()==names,'Unexpected or duplicate ZIP members'
    # Validate every member before creating any result file.
    for r in rows:
     with z.open(r['member']) as stream:assert hashlib.file_digest(stream,'sha256').hexdigest()==r['sha256']
     assert z.getinfo(r['member']).file_size==r['bytes']
    for r in missing:
     target=dest/r['path'];target.parent.mkdir(parents=True,exist_ok=True)
     with target.open('xb') as out,z.open(r['member']) as source:
      while chunk:=source.read(1024*1024):out.write(chunk)
 for r in rows:assert (dest/r['path']).stat().st_size==r['bytes'] and digest(dest/r['path'])==r['sha256']
 print(json.dumps(dict(passed=True,files=len(rows),restored=len(missing),destination=str(dest),archive_sha256=catalog['archive']['sha256'])))
if __name__=='__main__':main()
