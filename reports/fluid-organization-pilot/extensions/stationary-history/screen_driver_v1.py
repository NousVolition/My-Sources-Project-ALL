import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
from solver import Flow,initial,observe
ROOT=Path(__file__).resolve().parent
def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=48);p.add_argument('--nu',type=float,default=.025);p.add_argument('--power',type=float,default=.2);p.add_argument('--dt',type=float,default=.0125);p.add_argument('--end',type=float,default=60);p.add_argument('--seed',type=int,default=300000);a=p.parse_args()
    f=Flow(a.n,a.nu,a.power);h=initial(f,a.seed)[None];E0=float(f.inner(h,h)[0]/2);budget=np.zeros((1,2));rows=[];spectra=[];start=time.time();cfl=0
    folder=ROOT/'data/screen';folder.mkdir(parents=True,exist_ok=True);tag=f'n{a.n}_nu{a.nu}_dt{a.dt}_seed{a.seed}';path=folder/tag
    for step in range(round(a.end/a.dt)+1):
        if step%round(.5/a.dt)==0:
            _,s,sp=observe(f,h[0]);s.update(t=step*a.dt,cfl_max=cfl,balance=float(s['energy']-E0-budget[0,0]+budget[0,1]));rows.append(s);spectra.append(sp)
            if step%round(5/a.dt)==0:
                Path(str(path)+'.json').write_text(json.dumps({'config':vars(a),'elapsed':time.time()-start,'rows':rows},indent=2),encoding='utf-8')
                print(json.dumps({'seed':a.seed,'n':a.n,**{k:round(s[k],6) for k in ['t','energy','dissipation','re_lambda','keta','flux_k4','cfl_max','balance']},'elapsed':round(time.time()-start,1)}),flush=True)
        if step==round(a.end/a.dt):break
        h,b,_,c=f.step(h,a.dt);budget+=b;cfl=max(cfl,c)
        if not np.isfinite(h).all() or cfl>.6:raise RuntimeError('Unstable screen')
    np.savez_compressed(str(path)+'.npz',final=h[0],spectra=np.array(spectra),time=np.array([r['t'] for r in rows]))
    Path(str(path)+'.json').write_text(json.dumps({'config':vars(a),'elapsed':time.time()-start,'source_sha256':hashlib.sha256((ROOT/'solver.py').read_bytes()).hexdigest(),'rows':rows},indent=2),encoding='utf-8')
if __name__=='__main__':main()
