import argparse,concurrent.futures,subprocess,sys,os
from common import ROOT,protocol

def child(command):
 env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join(sys.path);env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1'
 log=ROOT/'data/logs'/('_'.join(map(str,command))+'.log');log.parent.mkdir(parents=True,exist_ok=True)
 with log.open('w',encoding='utf-8') as fh:r=subprocess.run([sys.executable,str(ROOT/'simulate.py'),*map(str,command)],env=env,stdout=fh,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError('Failed '+str(command)+'; '+log.read_text(encoding='utf-8')[-1800:])
 return command

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['calibration','test','controls']);ap.add_argument('--workers',type=int,default=8);a=ap.parse_args();p=protocol()
 if a.stage in ['test','controls']:
  from forecast import verify_freeze
  verify_freeze()
 if a.stage=='calibration':jobs=[[s] for kind in ['train','validation'] for s in p['splits'][kind]]
 elif a.stage=='test':jobs=[[s] for s in p['splits']['test']]
 else:jobs=[[s,'--label',label] for s in p['refinement_seeds'] for label in ['half','fine']]
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  futures=[pool.submit(child,job) for job in jobs]
  for i,future in enumerate(concurrent.futures.as_completed(futures)):
   print(f'{a.stage}: {i+1}/{len(jobs)} complete {future.result()}',flush=True)
