"""Preserve all original NPZ bytes in a split, verified raw archive."""
import argparse,json,tarfile,zipfile,hashlib
from pathlib import Path
from final_experiment import ROOT,dump
from verify_extension import FOLDERS
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as src:
        for block in iter(lambda:src.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()

def compact(out):
    code=out/'fluid-memory-extension-code-and-results.zip'
    with zipfile.ZipFile(code,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(ROOT.iterdir()):
            if p.is_file() and (p.suffix in ['.py','.json','.txt','.md','.html'] or p.name in ['.gitattributes','.gitignore']):z.write(p,p.name)
        for folder in ['results','figures']:
            for p in sorted((ROOT/folder).glob('*')):
                if p.is_file():z.write(p,p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(code) as z:
        if z.testzip():raise RuntimeError('Code archive CRC failed')
    dump(out/'compact-checksum.json',dict(file=code.name,bytes=code.stat().st_size,sha256=sha(code)))
    return code
def package(out):
    verified=json.loads((ROOT/'results'/'verification.json').read_text())
    if not verified['completion_and_integrity_passed']:raise RuntimeError('Integrity verification required')
    out.mkdir(parents=True,exist_ok=True);archive=out/'fluid-memory-extension-recordings.tar.gz'
    with tarfile.open(archive,'w:gz',compresslevel=0) as tar:
        for record in verified['records']:
            p=ROOT/record['file']
            if sha(p)!=record['sha256']:raise RuntimeError('Recording changed after verification')
            tar.add(p,arcname=record['file']);tar.add(p.with_suffix('.json'),arcname=Path(record['file']).with_suffix('.json').as_posix())
    parts=[]
    with archive.open('rb') as src:
        i=1
        while True:
            value=src.read(512*1024**2)
            if not value:break
            p=out/f'fluid-memory-extension-recordings.part{i:03d}.tar.gz';p.write_bytes(value)
            parts.append(dict(file=p.name,bytes=len(value),sha256=sha(p)));i+=1
    manifest=dict(format='Byte-split tar.gz; join listed parts in order',archive_name=archive.name,archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,
        parts=parts,records=verified['records'],recording_count=len(verified['records']),repository='NousVolition/My-Sources-Project-ALL',release_tag='fluid-memory-extension-2026-10-11',
        precision='Original float64/complex128 NPY and NPZ container bytes preserved',compression='Gzip level 0 around original compressed NPZ containers',request_id='2300ecd1-d208-4695-9b39-62b93d705892')
    dump(ROOT/'raw-recordings.json',manifest);dump(out/'raw-recordings.json',manifest)
    compact(out)
    print('Packaged '+str(len(manifest['records']))+' exact recordings into '+str(len(parts))+' release parts',flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();package(a.out)
