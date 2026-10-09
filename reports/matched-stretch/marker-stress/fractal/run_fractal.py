"""Finite random-carpet labels advected by the unchanged saved fluid flow."""
from pathlib import Path
import sys,os,json,hashlib
import numpy as np
from scipy import fft,ndimage
ROOT=Path(__file__).resolve().parent;PARENT=ROOT.parent;STUDY=PARENT.parent
sys.path.insert(0,str(PARENT))
from run import interp,rk4,wrap,sha,save
def carpet(depth=4):
    rng=np.random.default_rng(227);cells=np.zeros((1,2),dtype=int);stages={}
    for level in range(1,depth+1):
        out=[]
        for cell in cells:
            omit=int(rng.integers(9))
            out.extend(cell*3+[i,j] for i in range(3) for j in range(3) if i*3+j!=omit)
        cells=np.array(out,dtype=int);stages[level]=cells.copy()
    return stages
def main():
    assert not (ROOT/'state.json').exists(),'Inspect existing state before any rerun'
    stages=carpet();initial={d:np.c_[.375*((c+.5)/3**d-.5),np.zeros(len(c))] for d,c in stages.items() if d in (3,4)};seed=np.vstack([initial[3],initial[4]])
    np.savez_compressed(ROOT/'construction.npz',**{'level'+str(d):c for d,c in stages.items()})
    out={'protocol':json.loads((ROOT/'protocol.json').read_text()),'source_hashes':{n:sha(STUDY/n) for n in ('numerics.py','run_suite.py')},'run_script_sha256':sha(Path(__file__)),'construction_sha256':sha(ROOT/'construction.npz'),'counts':{'depth3':512,'depth4':4096},'fields':{},'runs':{}}
    configs=[('n64','baseline-n64-base',1,64),('n64_half','baseline-n64-base',1,128),('n128','baseline-n128-base',1,128),('coarse_saved_times','baseline-n64-base',2,256)]
    for name,rid,stride,steps in configs:
        save(ROOT/'state.json',{'status':'running','pid':os.getpid(),'active':name,'completed':list(out['runs'])});record=json.loads((STUDY/'runs'/rid/'result.json').read_text());n=record['job']['n'];assert record['series'][-1]['t']>=.1;cache={}
        def field(i):
            if i not in cache:
                p=STUDY/'runs'/rid/f'field-{i:03d}.npy';h=np.load(p,allow_pickle=False);assert np.isfinite(h).all();key=rid+'/'+p.name;digest=sha(p)
                if key in out['fields']:assert out['fields'][key]==digest
                out['fields'][key]=digest;cache[i]=fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1),workers=1)
            return cache[i]
        x=seed.copy();pos=[x.copy()];times=[0.];u=field(0);check=np.stack([ndimage.map_coordinates(a,(((seed+3)%6)*n/6).T,order=1,mode='grid-wrap') for a in u],axis=1);err=float(np.max(abs(check-interp(u,seed))));assert err<1e-11
        for i in range(0,10,stride):
            u0=field(i);u1=field(i+stride);duration=.01*stride;dt=duration/steps
            def v(p,t):
                a=np.clip(t/duration,0,1);return (1-a)*interp(u0,p)+a*interp(u1,p)
            for k in range(steps):x=rk4(x,dt,v,k*dt)
            assert np.isfinite(x).all();pos.append(x.copy());times.append(round(.01*(i+stride),8));cache.pop(i,None)
        np.savez_compressed(ROOT/(name+'.npz'),positions=pos,time=times)
        out['runs'][name]={'source':rid,'n':n,'dt':dt,'snapshot_dt':duration,'trajectory_sha256':sha(ROOT/(name+'.npz')),'initial_interpolation_check':err};save(ROOT/'results.json',out);print(name+' fractal labels complete through0.10',flush=True)
    assert all(sha(STUDY/n)==h for n,h in out['source_hashes'].items());save(ROOT/'state.json',{'status':'complete','configurations':4,'labels_per_configuration':4608,'fields':len(out['fields'])});print('Fractal tracking complete.',flush=True)
if __name__=='__main__':main()
