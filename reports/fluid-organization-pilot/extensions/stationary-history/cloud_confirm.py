"""Acquire the declared test ensemble and dependent controls without model fitting."""
import argparse,concurrent.futures,hashlib,json,os,platform,sys,zipfile
from pathlib import Path
from common import ROOT,save,now,protocol,digest
from forecast import verify_freeze
from batch import child

def archive(name,paths):
 out=ROOT/(name+'.zip')
 with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in sorted(set(paths)):z.write(p,p.relative_to(ROOT).as_posix())
 print(json.dumps(dict(archive=out.name,size=out.stat().st_size,sha256=digest(out))),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['base','control']);ap.add_argument('group',type=int);ap.add_argument('--download-directory');a=ap.parse_args();p=protocol();verify_freeze()
 receipt=ROOT/f'data/cloud/{a.action}_{a.group}_environment.json'
 environment=dict(started_utc=now(),action=a.action,group=a.group,platform=platform.platform(),python=sys.version,logical_cpus=os.cpu_count(),github_sha=os.environ.get('GITHUB_SHA'),github_run_id=os.environ.get('GITHUB_RUN_ID'),frozen_models_sha256=digest(ROOT/'results/frozen_models.json'))
 files=[receipt]
 if a.action=='base':
  seeds=p['splits']['test'][2*a.group:2*a.group+2];assert len(seeds)==2;environment.update(seeds=seeds,workers=2);save(str(receipt.relative_to(ROOT)),environment)
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   for command in pool.map(child,[[s] for s in seeds]):print('Completed',command,flush=True)
  for seed in seeds:files.extend((ROOT/'data/runs').glob(f'base_{seed}*'));files.append(ROOT/f'data/logs/{seed}.log')
 else:
  seed=p['refinement_seeds'][a.group//2];label=['half','fine'][a.group%2];environment.update(seeds=[seed],label=label,workers=1)
  folder=Path(a.download_directory);archives=list(folder.rglob('confirmation_base_*.zip'));assert len(archives)==1,archives
  with zipfile.ZipFile(archives[0]) as z:
   required=[f'data/runs/base_{seed}{suffix}' for suffix in ['.json','.npz','_checkpoints.npz']]
   for name in required:
    dest=ROOT/name;dest.parent.mkdir(parents=True,exist_ok=True);b=z.read(name)
    if dest.exists():assert dest.read_bytes()==b,name
    else:dest.write_bytes(b)
  meta=json.loads((ROOT/required[0]).read_text());assert digest(ROOT/required[1])==meta['sha256'];assert digest(ROOT/required[2])==meta['checkpoints']['sha256']
  environment['base_archive_sha256']=digest(archives[0]);save(str(receipt.relative_to(ROOT)),environment)
  child([seed,'--label',label]);files += [ROOT/f'data/runs/{label}_{seed}.json',ROOT/f'data/runs/{label}_{seed}.npz',ROOT/f'data/logs/{seed}_--label_{label}.log']
 environment['completed_utc']=now();save(str(receipt.relative_to(ROOT)),environment)
 archive(f'confirmation_{a.action}_{a.group}',files)
