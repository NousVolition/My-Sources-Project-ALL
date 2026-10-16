import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,platform,json,hashlib,ctypes,re
from pathlib import Path
import numpy,scipy,scipy.linalg,matplotlib
R=Path(__file__).resolve().parent
libraries=[]
for owner in [numpy,scipy]:
 for p in (Path(owner.__file__).parent.parent/(owner.__name__+'.libs')).glob('*openblas*'):
  b=p.read_bytes();names=sorted(set(re.findall(rb'[A-Za-z0-9_]*openblas_get_num_threads[A-Za-z0-9_]*',b)))
  if not names:continue
  lib=ctypes.CDLL(str(p));getter=getattr(lib,names[0].decode());getter.restype=ctypes.c_int;getter.argtypes=[]
  libraries.append(dict(owner=owner.__name__,name=p.name,threads=getter(),sha256=hashlib.sha256(b).hexdigest()))
assert libraries and all(v['threads']==1 for v in libraries)
out=dict(python=platform.python_version(),numpy=numpy.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,platform=platform.platform(),machine=platform.machine(),libraries=libraries)
(R/'data/environment.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
