"""Recreate the short reference and control trajectories in a fresh directory."""
from pathlib import Path
import argparse,json,shutil,subprocess,sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'source'))
from numerics import Flow

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    target=root/'study/join-isolation';target.mkdir(parents=True)
    for name in ('run.py','analyze.py','protocol.json','source/numerics.py'):
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/name,p)
    for half in (False,True):
        folder=root/'study/runs'/('baseline-n64-'+('half' if half else 'base'));folder.mkdir(parents=True)
        f=Flow(64,.001,workers=1);h=f.initial();cached=f.rhs(h);acc=np.zeros(5);rows=[]
        first=f.observe(h,0);substeps=103*(2 if half else 1);dt=.01/substeps
        for oi in range(5):
            if oi:
                for step in range(substeps):
                    previous=cached[1][0];h,v,_=f.step(h,dt,cached);cached=f.rhs(h)
                    v[0]=.5*dt*(previous+cached[1][0]);acc+=v
                    if not np.isfinite(h).all():raise FloatingPointError('Nonfinite reference field')
            rows.append(f.observe(h,oi*.01,acc,first));np.save(folder/f'field-{oi:03d}.npy',h)
        (folder/'result.json').write_text(json.dumps(dict(status='complete',dt=dt,series=rows),indent=2)+'\n')
    subprocess.run([sys.executable,'-u',str(target/'run.py')],check=True)
    subprocess.run([sys.executable,str(target/'analyze.py')],check=True)
    print('Reproduction saved in',target)

if __name__=='__main__':main()
