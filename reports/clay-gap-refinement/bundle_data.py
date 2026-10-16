"""Build and independently read back the release archive of exact raw arrays."""
from pathlib import Path
import argparse,json,hashlib,zipfile,shutil
R=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('destination',type=Path);a=ap.parse_args();dest=a.destination.resolve();dest.mkdir(exist_ok=False)
analysis=json.loads((R/'analysis.json').read_text());assert analysis['attempts']==32 and analysis['completed']==32,'Wait for all declared experiments'
files=[];prefix='reports/clay-gap-refinement/';archive=dest/'hug-gap-refinement-data.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_STORED) as out:
 for p in sorted(R.rglob('*.npz')):
  assert p.with_suffix('.json').exists();b=p.read_bytes();rel=p.relative_to(R).as_posix();info=zipfile.ZipInfo(prefix+rel,(1980,1,1,0,0,0));info.external_attr=0o100644<<16;out.writestr(info,b)
  files.append(dict(path=rel,member=prefix+rel,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),git_blob=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()))
with zipfile.ZipFile(archive) as z:
 assert z.namelist()==[f['member'] for f in files]
 for f in files:
  with z.open(f['member']) as stream:h=hashlib.file_digest(stream,'sha256').hexdigest()
  assert h==f['sha256'] and z.getinfo(f['member']).file_size==f['bytes']
b=archive.read_bytes();catalog=dict(schema=1,repository='NousVolition/My-Sources-Project-ALL',study='reports/clay-gap-refinement',tag='hug-gap-refinement-2026-10-15',archive=dict(name=archive.name,url='https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/hug-gap-refinement-2026-10-15/'+archive.name,bytes=len(b),sha256=hashlib.sha256(b).hexdigest()),files=files,archive_readback_passed=True)
(R/'data-release.json').write_text(json.dumps(catalog,indent=2)+'\n');shutil.copy2(R/'data-release.json',dest/'data-release.json')
(dest/'SHA256SUMS.txt').write_text(f"{catalog['archive']['sha256']}  {archive.name}\n{hashlib.sha256((dest/'data-release.json').read_bytes()).hexdigest()}  data-release.json\n")
print(json.dumps(dict(archive=str(archive),bytes=len(b),files=len(files),sha256=catalog['archive']['sha256'])))
