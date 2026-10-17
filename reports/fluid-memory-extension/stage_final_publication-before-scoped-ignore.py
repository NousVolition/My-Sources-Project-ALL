"""Freeze reviewable source/results and release assets after the raw roundtrip."""
import json,shutil,hashlib
from pathlib import Path
from final_experiment import ROOT,digest,dump
from verify_extension import FOLDERS

def main():
    outputs=ROOT.parents[1]/'outputs'/'fluid-memory-extension';bundle=outputs/'final-release'
    verified=json.loads((ROOT/'results'/'verification.json').read_text())
    ready=json.loads((ROOT/'results'/'release-ready.json').read_text())
    proof=json.loads((ROOT/'results'/'raw-roundtrip.json').read_text())
    if not (verified['completion_and_integrity_passed'] and ready['raw_roundtrip_passed'] and proof['passed']):raise RuntimeError('Verified finished study required')
    stage=outputs/'final-publication'
    if stage.exists():raise RuntimeError('Final publication staging directory already exists; do not overwrite a frozen snapshot')
    stage.mkdir(parents=True)
    files=[p for p in ROOT.iterdir() if p.is_file() and (p.suffix in ['.py','.json','.txt','.md','.html'] or p.name=='.gitattributes')]
    for folder in ['results','figures']+FOLDERS:
        files.extend(p for p in (ROOT/folder).glob('*') if p.is_file() and p.suffix in ['.json','.png'])
    manifest=[]
    for p in sorted(files):
        rel=p.relative_to(ROOT).as_posix();dest=stage/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        value=dest.read_bytes()
        if len(value)>=100*1024**2:raise RuntimeError('File belongs in release assets instead: '+rel)
        gitsha=hashlib.sha1(b'blob '+str(len(value)).encode()+b'\0'+value).hexdigest()
        manifest.append(dict(file=rel,bytes=len(value),sha256=digest(dest),git_blob_sha=gitsha))
    raw=json.loads((bundle/'raw-recordings.json').read_text());assets=list(raw['parts'])
    for name in ['fluid-memory-extension-code-and-results.zip','raw-recordings.json','raw-roundtrip.json','compact-checksum.json']:
        p=bundle/name;assets.append(dict(file=name,bytes=p.stat().st_size,sha256=digest(p)))
    for asset in assets:
        p=bundle/asset['file']
        if p.stat().st_size!=asset['bytes'] or digest(p)!=asset['sha256']:raise RuntimeError('Asset integrity changed: '+asset['file'])
    dump(outputs/'final-publication-manifest.json',dict(stage=str(stage),files=manifest,recording_count=verified['recording_count'],request_id=ready['request_id']))
    dump(outputs/'release-assets.json',dict(bundle=str(bundle),assets=assets,tag=raw['release_tag'],total_bytes=sum(a['bytes'] for a in assets),request_id=ready['request_id']))
    print(f'Frozen {len(manifest)} GitHub files and {len(assets)} release assets for {verified["recording_count"]} verified recordings.',flush=True)

if __name__=='__main__':main()
