"""Menger-only reproduction entry point; requires the original saved fluid fields.

The archived tracks were extracted losslessly from a larger batch of independent
passive labels. This script reproduces only its valid Menger subset. Its hash is
not the hash of the historical runner recorded in menger-results.json.
"""
from pathlib import Path
import sys,os,json
import numpy as np
from scipy import fft,ndimage
ROOT=Path(__file__).resolve().parent;STUDY=ROOT.parents[1]
sys.path.insert(0,str(ROOT.parent))
from run import interp,rk4,sha,save
def construction():
    cells=np.zeros((1,3),int);stages={}
    offsets=np.array([v for v in np.ndindex(3,3,3) if sum(a==1 for a in v)<=1])
    assert len(offsets)==20
    for depth in range(1,4):
        cells=(3*cells[:,None]+offsets[None]).reshape(-1,3);stages[depth]=cells.copy()
    return stages
def initial_positions():
    stages=construction()
    return np.vstack([.375*((stages[d]+.5)/3**d-.5) for d in (2,3)])
def main():
    assert not (ROOT/'menger-results.json').exists(),'Existing Menger results: inspect before any reproduction'
    stages=construction();np.savez_compressed(ROOT/'menger-construction.npz',**{'menger'+str(d):p for d,p in stages.items()})
    seed=initial_positions();groups=[{'name':'menger2','family':'menger','depth':2,'base':3,'start':0,'stop':400,'count':400},{'name':'menger3','family':'menger','depth':3,'base':3,'start':400,'stop':8400,'count':8000}]
    out={'scope':'Menger-only passive-label reproduction','groups':groups,'geometry_sha256':sha(ROOT/'menger-construction.npz'),'run_script_sha256':sha(Path(__file__)),'source_hashes':{n:sha(STUDY/n) for n in ('numerics.py','run_suite.py')},'runs':{},'fields':{}}
    for name,rid,stride,steps in [('n64','baseline-n64-base',1,64),('n64_half','baseline-n64-base',1,128),('n128','baseline-n128-base',1,128),('coarse_saved_times','baseline-n64-base',2,256)]:
        save(ROOT/'menger-state.json',{'status':'running','pid':os.getpid(),'active':name,'completed':list(out['runs'])})
        n=json.loads((STUDY/'runs'/rid/'result.json').read_text())['job']['n'];cache={}
        def field(i):
            if i not in cache:
                p=STUDY/'runs'/rid/f'field-{i:03d}.npy';h=np.load(p,allow_pickle=False);assert np.isfinite(h).all();key=rid+'/'+p.name;digest=sha(p)
                if key in out['fields']:assert out['fields'][key]==digest
                out['fields'][key]=digest;cache[i]=fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1),workers=1)
            return cache[i]
        x=seed.copy();pos=[x.copy()];times=[0.];u=field(0)
        independent=np.stack([ndimage.map_coordinates(a,(((seed+3)%6)*n/6).T,order=1,mode='grid-wrap') for a in u],axis=1)
        err=float(np.max(abs(independent-interp(u,seed))));assert err<1e-11
        for i in range(0,10,stride):
            u0=field(i);u1=field(i+stride);duration=.01*stride;dt=duration/steps
            def v(p,t):
                a=np.clip(t/duration,0,1);return (1-a)*interp(u0,p)+a*interp(u1,p)
            for k in range(steps):x=rk4(x,dt,v,k*dt)
            assert np.isfinite(x).all();pos.append(x.copy());times.append(round(.01*(i+stride),8));cache.pop(i,None)
        np.savez_compressed(ROOT/('menger-'+name+'.npz'),positions=pos,time=times)
        out['runs'][name]={'source':rid,'n':n,'dt':dt,'snapshot_dt':duration,'trajectory_sha256':sha(ROOT/('menger-'+name+'.npz')),'independent_interpolation_error':err}
        save(ROOT/'menger-results.json',out);print(name+' Menger tracks complete through 0.10',flush=True)
    assert all(sha(STUDY/n)==s for n,s in out['source_hashes'].items())
    save(ROOT/'menger-state.json',{'status':'complete','labels_per_configuration':len(seed),'configurations':4})
if __name__=='__main__':main()
