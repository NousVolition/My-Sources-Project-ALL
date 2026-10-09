"""Package completed scientific artifacts without bundling runtime dependencies."""
from pathlib import Path
import argparse, hashlib, json, shutil, time, zipfile
ROOT=Path(__file__).resolve().parent
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2),encoding='utf-8')
def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    verify=json.loads((ROOT/'verification.json').read_text());assert verify['passed']
    assert verify['completed_runs']+len(verify['failed_runs'])==verify['planned_runs'],'Every planned case must have an explicit completed or failed record'
    completed=[];interrupted=[]
    for folder in sorted((ROOT/'data').glob('*')):
        if not folder.is_dir():continue
        if (folder/'result.json').exists():completed.append(folder.name)
        elif (folder/'progress.json').exists():interrupted.append(dict(run=folder.name,last_saved=json.loads((folder/'progress.json').read_text()),status='intentionally interrupted when switching to the validated GPU backend'))
    dump(ROOT/'cpu-execution.json',dict(completed_controls=completed,interrupted_runs=interrupted,primary_matrix=verify['primary_data'],note='For the published GPU experiment, CPU runs are extra controls, not additional independent starts; partial CPU checkpoints are retained but not counted as completed trajectories.'))
    raw=[]
    for f in sorted(ROOT.rglob('*.npz')):
        with zipfile.ZipFile(f) as z:
            bad=z.testzip()
            if bad:raise ValueError(f'Corrupt array archive {f}: {bad}')
        raw.append(dict(path=f.relative_to(ROOT).as_posix(),bytes=f.stat().st_size,sha256=sha(f)))
    dump(ROOT/'raw-manifest.json',dict(files=raw,count=len(raw),bytes=sum(f['bytes'] for f in raw),all_npz_zip_crc_passed=True))
    dump(ROOT/'provenance.json',dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),repository='NousVolition/My-Sources-Project-ALL',study='reports/navier-stokes-stress',
        primary_equations='unforced positive-viscosity 3D incompressible Navier-Stokes',source_hashes={f.name:sha(f) for f in ROOT.glob('*.py')},protocol_sha256=sha(ROOT/'protocol.json'),gpu_amendment_sha256=sha(ROOT/'gpu-amendment.json'),
        verification=verify,molecular_runs_in_this_experiment=0,untrusted_late_results='Retained; resolution failures cannot support continuum singularity claims.'))
    files=[f for f in ROOT.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.name not in ['manifest.json']]
    compact=[f for f in files if f.suffix not in ['.npz','.log'] and f.name!='progress.json']
    archive=a.out/'navier-stokes-code-and-results.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for f in compact:z.write(f,f.relative_to(ROOT).as_posix())
    # Bound each upload to about 500 MiB; each part contains complete files.
    parts=[];groups=[];current=[];size=0
    for f in files:
        if size+f.stat().st_size>500*1024**2 and current:groups.append(current);current=[];size=0
        current.append(f);size+=f.stat().st_size
    if current:groups.append(current)
    for i,group in enumerate(groups):
        path=a.out/(f'navier-stokes-full-part{i+1:02d}.zip' if len(groups)>1 else 'navier-stokes-full.zip')
        with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=3) as z:
            for f in group:z.write(f,f.relative_to(ROOT).as_posix(),compress_type=zipfile.ZIP_STORED if f.suffix=='.npz' else zipfile.ZIP_DEFLATED)
        parts.append(path)
    shutil.copy2(ROOT/'report.html',a.out/'report.html')
    for name in ['stress.png','convergence.png','perturbations.png','taylor-green.png','budgets.png','vorticity-slices.png','README.md','summary.json','verification.json']:shutil.copy2(ROOT/name,a.out/name)
    delivery=dict(primary_runs=verify['completed_runs'],raw_arrays=len(raw),assets=[dict(file=f.name,bytes=f.stat().st_size,sha256=sha(f)) for f in [archive,*parts,a.out/'report.html']],
        extract='Extract all full archive parts into the same folder; parts contain separate files, not a split byte stream. Compact archive excludes NPZ arrays.',
        study_url='https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/navier-stokes-stress',release_url='https://github.com/NousVolition/My-Sources-Project-ALL/releases/tag/navier-stokes-stress-2026-10-09')
    dump(a.out/'delivery.json',delivery);print(json.dumps(delivery,indent=2))
if __name__=='__main__':main()
