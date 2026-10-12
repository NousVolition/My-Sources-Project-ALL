"""Execute reciprocal long runs and independent accuracy controls sequentially."""
import sys,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def execute(*args):subprocess.run([sys.executable,*map(str,args)],cwd=ROOT,check=True)
for parent,name,end in [('hard48-two','long48-two',.5),('hard72-two','long72-two',.5),('stress128-two','long24-two',1.5)]:
    execute('continue_run.py',parent,name,'--end',end,'--seconds',2400)
execute('solver/experiment.py','fresh24-one','--n',24,'--layers',3,'--seconds',600)
execute('solver/experiment.py','fresh24-two','--n',24,'--layers',3,'--two','--seconds',600)
for mode in ['one','two']:
    execute('continue_run.py','stress128-'+mode,'radau24-'+mode,'--end',.6,'--seconds',1200,'--method','Radau','--rtol',2e-10)
