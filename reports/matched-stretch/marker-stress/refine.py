"""One additional tracking-step control after the 128-grid path screen failed."""
from run import *
def main():
    out=json.loads((ROOT/'results.json').read_text());assert json.loads((ROOT/'state.json').read_text())['status']=='complete'
    name='grid128_refined';assert name not in out['runs'] and not (ROOT/(name+'.npz')).exists()
    reason='The 128-grid 0.0003125-to0.00015625 tracking comparison exceeded the declared path-discrepancy screen. Add exactly one half step to distinguish tracking error from the much larger grid sensitivity.'
    protocol=json.loads((ROOT/'protocol.json').read_text());protocol['additional_tracking_control']={'dt':.000078125,'reason':reason};save(ROOT/'protocol.json',protocol);out['protocol']=protocol
    save(ROOT/'state.json',{'status':'running','pid':os.getpid(),'active':name,'completed':list(out['runs'])})
    seed,edges=seeds(json.loads((ROOT/'source-records/baseline-n64-base.json').read_text())['series'][0]['peak_location'])
    neighbors=[]
    for g in (0,1):
        p=seed[g*125:(g+1)*125];dd=np.linalg.norm(wrap(p[:,None]-p[None,:]),axis=2);np.fill_diagonal(dd,np.inf);neighbors.append(np.argsort(dd,axis=1,kind='stable')[:,:6])
    folder=STUDY/'runs/baseline-n128-base';cache={}
    def field(i):
        if i not in cache:
            path=folder/f'field-{i:03d}.npy';assert sha(path)==out['fields']['baseline-n128-base/'+path.name];h=np.load(path,allow_pickle=False);cache[i]=fft.irfftn(h,s=(128,)*3,axes=(-3,-2,-1),workers=1)
        return cache[i]
    x=seed.copy();times=[0.];pos=[x.copy()];vel=[interp(field(0),x)];rows=[metrics(x,vel[-1],edges,neighbors)];dt=.000078125
    for i in range(30):
        u0=field(i);u1=field(i+1)
        def v(p,t):
            a=np.clip(t/.01,0,1);return (1-a)*interp(u0,p)+a*interp(u1,p)
        for k in range(128):x=rk4(x,dt,v,k*dt)
        assert np.isfinite(x).all() and np.array_equal(x[500:],x[:5]);times.append(round(.01*(i+1),8));pos.append(x.copy());vel.append(interp(u1,x));rows.append(metrics(x,vel[-1],edges,neighbors));cache.pop(i,None)
    np.savez_compressed(ROOT/(name+'.npz'),time=times,positions=pos,velocities=vel)
    record=out['runs']['grid128'];out['runs'][name]={**record,'tracking_dt':dt,'times':times,'metrics':rows,'trajectory_sha256':sha(ROOT/(name+'.npz')),'refinement_reason':reason}
    out['refine_script_sha256']=sha(Path(__file__));assert all(sha(STUDY/n)==s for n,s in out['source_hashes'].items());save(ROOT/'results.json',out);save(ROOT/'state.json',{'status':'complete','pid':os.getpid(),'completed':list(out['runs']),'saved_fields_verified':len(out['fields'])});print(name+' complete through0.30',flush=True)
if __name__=='__main__':main()
