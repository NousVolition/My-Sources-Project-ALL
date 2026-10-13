"""Hosted acquisition; immutable ZIPs connect pilot, freeze and confirmation."""
import argparse,json,os,platform,sys,zipfile,subprocess
from pathlib import Path
from common import ROOT,save,now,protocol,digest

def archive(name,files):
    dest=ROOT/(name+'.zip')
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(set(files)):z.write(p,p.relative_to(ROOT).as_posix())
    print(json.dumps(dict(archive=dest.name,size=dest.stat().st_size,sha256=digest(dest))),flush=True)
    return dest

def unpack(folder):
    for path in sorted(Path(folder).rglob('*.zip')):
        with zipfile.ZipFile(path) as z:
            for name in z.namelist():
                dest=(ROOT/name).resolve();assert dest.is_relative_to(ROOT.resolve()),name
                data=z.read(name);dest.parent.mkdir(parents=True,exist_ok=True)
                if dest.exists():assert dest.read_bytes()==data,name
                else:dest.write_bytes(data)

def run(command):
    tag='_'.join(map(str,command));log=ROOT/'data/logs'/(tag+'.log');log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w',encoding='utf-8') as fh:
        result=subprocess.run([sys.executable,*map(str,command)],cwd=ROOT,stdout=fh,stderr=subprocess.STDOUT)
    print(log.read_text(encoding='utf-8')[-2500:],flush=True)
    if result.returncode:raise RuntimeError('Failed '+tag)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['pilot','pilot-review','calibration','freeze','test','control','evaluate']);ap.add_argument('index',type=int,nargs='?',default=0);ap.add_argument('--input',default='downloaded');a=ap.parse_args();p=protocol()
    environment=dict(started_utc=now(),stage=a.stage,index=a.index,github_sha=os.environ.get('GITHUB_SHA'),github_run_id=os.environ.get('GITHUB_RUN_ID'),platform=platform.platform(),python=sys.version,cpus=os.cpu_count())
    if a.stage!='pilot':unpack(a.input)
    if a.stage=='pilot':
        run(['pilot.py',a.index]);files=list((ROOT/'data/pilot').glob(f'pilot_{a.index}*'))+list((ROOT/'data/logs').glob('*'))
        if a.index==p['response_pilot']['case']:files+=list((ROOT/'data/pilot').glob('response_*'))
    elif a.stage=='pilot-review':
        from checks import screening
        result=screening();save('results/pilot_response_checks.json',result);assert result['passed'],result
        files=list((ROOT/'data/pilot').glob('*.json'))+list((ROOT/'data/pilot').glob('response_*.npz'))+[ROOT/'results/pilot_response_checks.json']
    elif a.stage=='calibration':
        seeds=p['splits']['train']+p['splits']['validation'];seed=seeds[a.index];run(['simulate.py',seed]);files=list((ROOT/'data/runs').glob(f'base_{seed}*'))+list((ROOT/'data/logs').glob('*'))
    elif a.stage=='freeze':
        from forecast import freeze
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=1):freeze()
        files=list((ROOT/'data/runs').glob('*'))+list((ROOT/'data/pilot').glob('*.json'))+list((ROOT/'data/pilot').glob('response_*.npz'))+list((ROOT/'results').glob('*.json'))+list((ROOT/'data/cloud').glob('*'))+list((ROOT/'data/logs').glob('*'))
    elif a.stage=='test':
        from forecast import verify_freeze
        verify_freeze();seed=p['splits']['test'][a.index];environment['frozen_models_sha256']=digest(ROOT/'results/frozen_models.json');run(['simulate.py',seed]);files=list((ROOT/'data/runs').glob(f'base_{seed}*'))+list((ROOT/'data/logs').glob('*'))
    elif a.stage=='control':
        from forecast import verify_freeze
        verify_freeze();seed=p['refinement_seeds'][a.index//2];label=['half','fine'][a.index%2];run(['simulate.py',seed,'--label',label]);files=list((ROOT/'data/runs').glob(f'{label}_{seed}*'))+list((ROOT/'data/logs').glob('*'))
    else:
        from forecast import evaluate
        from checks import verify_all
        from threadpoolctl import threadpool_limits
        checks=verify_all()
        with threadpool_limits(limits=1):evaluate()
        files=[q for folder in ['data/runs','data/pilot','data/cloud','data/logs','results'] for q in (ROOT/folder).rglob('*') if q.is_file()]+[ROOT/'data/predictions.npz']
    environment['completed_utc']=now();receipt=ROOT/f'data/cloud/{a.stage}_{a.index}_environment.json';save(str(receipt.relative_to(ROOT)),environment);files.append(receipt)
    archive(f'higher_re_{a.stage}_{a.index}',files)

if __name__=='__main__':main()
