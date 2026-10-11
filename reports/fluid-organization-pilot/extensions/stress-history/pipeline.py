"""Fresh run sequence; test fields are generated only after model freeze."""
import os,subprocess,sys
from common import ROOT
env=os.environ.copy();env.update(PYTHONPATH=os.pathsep.join(sys.path),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONIOENCODING='utf-8')
for script,args in [('batch.py',['calibration']),('forecast.py',['freeze']),('batch.py',['test']),('batch.py',['controls']),('verify.py',[]),('forecast.py',['evaluate']),('analyze.py',[])]:
    print('Starting '+script+' '+' '.join(args),flush=True);subprocess.run([sys.executable,str(ROOT/script)]+args,check=True,env=env)
