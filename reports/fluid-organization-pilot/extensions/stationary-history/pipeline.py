"""Run from this directory after installing requirements. Never refits an existing freeze."""
import argparse,subprocess,sys,os
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def run(*args):
 env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join(sys.path)
 subprocess.run([sys.executable,*map(str,args)],cwd=ROOT,env=env,check=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('phase',choices=['calibration','confirmation']);p.add_argument('--workers',type=int,default=8);a=p.parse_args()
 if a.phase=='calibration':
  run('batch.py','calibration','--workers',a.workers);run('forecast.py','freeze')
 else:
  run('batch.py','test','--workers',a.workers);run('batch.py','controls','--workers',a.workers)
  run('checks.py');run('forecast.py','evaluate');run('analyze.py');run('report.py');run('manifest.py');run('audit.py')
