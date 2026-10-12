"""Bounded two-run worker for the explicitly listed calibration seed groups."""
import argparse,concurrent.futures,json,os,platform,subprocess,sys
from pathlib import Path
from common import ROOT,save,now,protocol
from batch import child

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('group',type=int);a=ap.parse_args()
 config=json.loads((ROOT/'cloud_calibration.json').read_text(encoding='utf-8'));seeds=config['groups'][str(a.group)]
 assert len(seeds)==2 and set(seeds)<=set(protocol()['splits']['train']+protocol()['splits']['validation'])
 save(f'data/cloud/group_{a.group}_environment.json',dict(started_utc=now(),platform=platform.platform(),python=sys.version,logical_cpus=os.cpu_count(),group=a.group,seeds=seeds,github_sha=os.environ.get('GITHUB_SHA'),github_run_id=os.environ.get('GITHUB_RUN_ID'),workers=2))
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  for command in pool.map(child,[[s] for s in seeds]):print('Completed',command,flush=True)
 import zipfile,hashlib
 target=ROOT/f'calibration_group_{a.group}.zip'
 with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for seed in seeds:
   for p in sorted((ROOT/'data/runs').glob(f'base_{seed}*')):z.write(p,p.relative_to(ROOT).as_posix())
   log=ROOT/'data/logs'/f'{seed}.log';z.write(log,log.relative_to(ROOT).as_posix())
  p=ROOT/f'data/cloud/group_{a.group}_environment.json';z.write(p,p.relative_to(ROOT).as_posix())
 print(json.dumps(dict(archive=target.name,size=target.stat().st_size,sha256=hashlib.sha256(target.read_bytes()).hexdigest())))

