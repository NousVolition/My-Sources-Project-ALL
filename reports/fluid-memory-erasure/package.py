"""Publish the original compressed NPZ bytes unchanged in a split tar.gz.

Gzip uses stored blocks: the NPZ recordings are already compressed. Avoiding
redundant recompression speeds publication and preserves container hashes too.
"""
import argparse,hashlib,json,lzma,tarfile,tempfile,time,zipfile
from pathlib import Path
from run import ROOT,digest,dump


def reencode_npz(source,destination):
    members={}
    with zipfile.ZipFile(source,'r') as zin,zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_STORED,allowZip64=True) as zout:
        for name in zin.namelist():
            value=zin.read(name);members[name]=hashlib.sha256(value).hexdigest()
            info=zipfile.ZipInfo(name,date_time=(2026,10,10,0,0,0));info.compress_type=zipfile.ZIP_STORED
            zout.writestr(info,value)
    with zipfile.ZipFile(destination) as z:
        if any(hashlib.sha256(z.read(k)).hexdigest()!=v for k,v in members.items()):raise RuntimeError('Lossless reencoding failed')
    return members


def raw_archive(destination):
    records=[];arch=destination/'fluid-memory-recordings.tar.gz'
    folders=['data','pilot-data','matched-sham-data','delayed-probe-data']
    with tempfile.TemporaryDirectory(prefix='raw-reencode-',dir=ROOT) as td:
        tmp=Path(td)
        if not tmp.resolve().is_relative_to(ROOT.resolve()):raise RuntimeError('Temporary directory outside study')
        with tarfile.open(arch,'w:gz',compresslevel=0) as tar:
            for folder in folders:
                if folder in ['matched-sham-data','delayed-probe-data']:
                    deadline=time.monotonic()+1200;execution=ROOT/folder/'execution.json'
                    while not (execution.exists() and (lambda d:d['completed']==d['planned'])(json.loads(execution.read_text()))):
                        if time.monotonic()>deadline:raise RuntimeError('Incomplete simulation folder: '+folder)
                        time.sleep(2)
                for source in sorted((ROOT/folder).glob('s*.npz')):
                    old=json.loads(source.with_suffix('.json').read_text())
                    if digest(source)!=old['raw_sha256']:raise RuntimeError('Original recording failed SHA check')
                    with zipfile.ZipFile(source) as z:members={k:hashlib.sha256(z.read(k)).hexdigest() for k in z.namelist()}
                    packaged=old['raw_sha256']
                    relative=f'{folder}/{source.name}'
                    tar.add(source,arcname=relative);tar.add(source.with_suffix('.json'),arcname=str(Path(relative).with_suffix('.json')).replace('\\','/'))
                    records.append(dict(file=relative,sha256=packaged,original_npz_sha256=digest(source),array_members_sha256=members,
                                        original_bytes=source.stat().st_size,published_npz_bytes=source.stat().st_size))
                    print('Archived '+relative,flush=True)
                for source in sorted((ROOT/folder).glob('*.json')):
                    if not source.with_suffix('.npz').exists():tar.add(source,arcname=f'{folder}/{source.name}')
    parts=[];part_size=128*1024**2
    with arch.open('rb') as src:
        number=1
        while True:
            value=src.read(part_size)
            if not value:break
            name=f'fluid-memory-recordings.part{number:02d}.tar.gz';path=destination/name;path.write_bytes(value)
            parts.append(dict(file=name,bytes=len(value),sha256=digest(path)));number+=1
    manifest=dict(format='Byte-split tar.gz; join parts in listed order',archive_name=arch.name,archive_sha256=digest(arch),archive_bytes=arch.stat().st_size,
                  parts=parts,records=records,recording_count=len(records),precision='Original float64/complex128 NPY bytes, unchanged',
                  original_npz_bytes=sum(r['original_bytes'] for r in records),array_content_verification='Original NPZ container bytes unchanged; each original NPY member SHA256 also recorded',
                  release_tag='fluid-memory-erasure-2026-10-10',repository='NousVolition/My-Sources-Project-ALL')
    manifest['outer_compression']='Gzip stored blocks, level 0; original NPZ ZIP_DEFLATE compression unchanged'
    dump(ROOT/'raw-recordings.json',manifest);dump(destination/'raw-recordings.json',manifest)
    return manifest


def compact(destination):
    path=destination/'fluid-memory-code-and-results.zip'
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for file in sorted(ROOT.iterdir()):
            if file.is_file() and file.suffix in ['.py','.json','.txt','.md','.html'] and file.name!='github-publication.json':z.write(file,file.name)
        for folder in ['results','pilot-results','figures','freeze-source']:
            for file in sorted((ROOT/folder).rglob('*')):
                if file.is_file():z.write(file,str(file.relative_to(ROOT)).replace('\\','/'))
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:raise RuntimeError('Compact archive CRC failed')
    return dict(file=path.name,bytes=path.stat().st_size,sha256=digest(path))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--raw',action='store_true');ap.add_argument('--raw-only',action='store_true');args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if args.raw:raw_archive(args.out)
    if not args.raw_only:
        info=compact(args.out);dump(args.out/'compact-checksum.json',info);print(json.dumps(info,indent=2))
