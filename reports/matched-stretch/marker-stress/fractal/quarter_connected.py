"""Correct screenshot recipe: two connected survivors per interval, one coin per stage."""
from pathlib import Path
import sys,os,json
import numpy as np
from scipy import fft,ndimage
ROOT=Path(__file__).resolve().parent;STUDY=ROOT.parents[1]
sys.path.insert(0,str(ROOT.parent))
from run import interp,rk4,sha,save
def build():
    rng=np.random.default_rng(551);intervals=np.array([[0.,1.]]);stages={};flips=[]
    for depth in range(1,11):
        quarter=int(rng.integers(2,4));flips.append(quarter);out=[]
        for a,b in intervals:
            length=b-a;out.extend([[a,a+(quarter-1)*length/4],[a+quarter*length/4,b]])
        intervals=np.array(out);stages[depth]=intervals.copy()
        assert len(intervals)==2**depth and abs(np.sum(intervals[:,1]-intervals[:,0])-(.75)**depth)<1e-13
    return stages,flips
def main():
    assert not (ROOT/'quarter-state.json').exists(),'Inspect existing state first'
    stages,flips=build();np.savez_compressed(ROOT/'quarter-construction.npz',**{'level'+str(d):x for d,x in stages.items()})
    groups=[{'name':'quarter8','start':0,'stop':256,'depth':8,'base':2},{'name':'quarter10','start':256,'stop':1280,'depth':10,'base':2}];seed=np.vstack([np.c_[.375*(stages[d].mean(axis=1)-.5),np.zeros((2**d,2))] for d in (8,10)])
    out={'scope':'Corrected connected-interval quarter-removal set; no fluid change','groups':groups,'seed':551,'deleted_quarter_by_stage':flips,'rule':'One coin per stage, delete second or third quarter from each connected interval, leaving two intervals of lengths1/4 and1/2','limiting_dimension':float(np.log((1+np.sqrt(5))/2)/np.log(2)),'discarded_interpretation':'other_patterns.py treated adjacent quarters as separate children. Its quarter subset is retained locally as an excluded construction error; its independent Menger tracks are unaffected.','geometry_sha256':sha(ROOT/'quarter-construction.npz'),'source_hashes':{n:sha(STUDY/n) for n in ('numerics.py','run_suite.py')},'run_script_sha256':sha(Path(__file__)),'runs':{},'fields':{}}
    for name,rid,stride,steps in [('n64','baseline-n64-base',1,64),('n64_half','baseline-n64-base',1,128),('n128','baseline-n128-base',1,128),('coarse_saved_times','baseline-n64-base',2,256)]:
        save(ROOT/'quarter-state.json',{'status':'running','pid':os.getpid(),'active':name,'completed':list(out['runs'])});n=json.loads((STUDY/'runs'/rid/'result.json').read_text())['job']['n'];cache={}
        def field(i):
            if i not in cache:
                p=STUDY/'runs'/rid/f'field-{i:03d}.npy';h=np.load(p,allow_pickle=False);assert np.isfinite(h).all();key=rid+'/'+p.name;digest=sha(p)
                if key in out['fields']:assert out['fields'][key]==digest
                out['fields'][key]=digest;cache[i]=fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1),workers=1)
            return cache[i]
        x=seed.copy();times=[0.];pos=[x.copy()];u=field(0);check=np.stack([ndimage.map_coordinates(a,(((seed+3)%6)*n/6).T,order=1,mode='grid-wrap') for a in u],axis=1);err=float(np.max(abs(check-interp(u,seed))));assert err<1e-11
        for i in range(0,10,stride):
            u0=field(i);u1=field(i+stride);duration=.01*stride;dt=duration/steps
            def v(p,t):
                a=np.clip(t/duration,0,1);return (1-a)*interp(u0,p)+a*interp(u1,p)
            for k in range(steps):x=rk4(x,dt,v,k*dt)
            assert np.isfinite(x).all();pos.append(x.copy());times.append(round(.01*(i+stride),8));cache.pop(i,None)
        np.savez_compressed(ROOT/('quarter-'+name+'.npz'),positions=pos,time=times);out['runs'][name]={'source':rid,'n':n,'dt':dt,'snapshot_dt':duration,'trajectory_sha256':sha(ROOT/('quarter-'+name+'.npz')),'initial_interpolation_check':err};save(ROOT/'quarter-results.json',out);print(name+' corrected quarter tracks complete through0.10',flush=True)
    assert all(sha(STUDY/n)==s for n,s in out['source_hashes'].items());save(ROOT/'quarter-state.json',{'status':'complete','labels_per_configuration':1280,'configurations':4})
if __name__=='__main__':main()
