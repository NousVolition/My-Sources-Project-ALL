"""Short unforced counterfactuals; original trajectories are reused, never overwritten."""
from pathlib import Path
import hashlib, json, os, sys, time
import numpy as np

ROOT = Path(__file__).resolve().parent
STUDY = ROOT.parent
sys.path.insert(0, str(ROOT/'source'))
from numerics import Flow, L

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    for attempt in range(100):
        try: tmp.replace(path); return
        except PermissionError:
            if attempt == 99: raise
            time.sleep(.1)

def periodic_background(f):
    k = 2*np.pi/L
    a = np.arange(f.n)*f.dx-L/2
    xyz = np.meshgrid(a,a,a,indexing='ij',sparse=True)
    r2 = sum(2*(1-np.cos(k*x))/k**2 for x in xyz)
    env = np.exp(-.04*r2)
    raw = np.stack([c*np.sin(k*x)/k*env for c,x in zip((-40.,-40.,80.),xyz)])
    return f.project(f.hat(raw))

def initial(f, original, variant):
    bg = f.initial(spin_factor=0.)
    bg[:,0,0,0] = original[:,0,0,0]
    if variant == 'without-tube': return bg
    # Exact same sampled tube and mean as the original at this grid.
    h = periodic_background(f)+(original-bg)
    h[:,0,0,0] = original[:,0,0,0]
    return h

def extra(f,h):
    a = np.arange(f.n)*f.dx-L/2
    wh = f.curl(h); w = f.real(wh); mag = np.sqrt(np.sum(w*w,axis=0))
    ix = np.unravel_index(mag.argmax(),mag.shape)
    x,y,z = np.meshgrid(a,a,a,indexing='ij',sparse=True)
    central = (x*x+y*y<=.4**2)&(abs(z)<=2.5)
    band = (abs(x)>=2.625)|(abs(y)>=2.625)|(abs(z)>=2.625)
    c = f.n//2
    grad = [[float(f.real(1j*f.k[j]*h[i])[c,c,c]) for j in range(3)] for i in range(3)]
    return dict(peak_distance_to_nearest_join=float(min(3-abs(a[j]) for j in ix)),
        interior_central_cylinder_W=float(mag[central].max()),
        fixed_join_band_W=float(mag[band].max()), center_gradient=grad)

def run(job, protocol):
    folder=ROOT/'runs'/job['id'];folder.mkdir(parents=True,exist_ok=True)
    path=folder/'result.json'
    hashes={name:sha(ROOT/name) for name in ('run.py','source/numerics.py','protocol.json')}
    if path.exists():
        old=json.loads(path.read_text())
        if old['source_hashes']!=hashes: raise RuntimeError('Source changed; refusing resume')
        if old['status'] in ('complete','failed_check'): return old
    f=Flow(64,.001,workers=1)
    orig=np.load(STUDY/'runs/baseline-n64-base/field-000.npy')
    assert np.sqrt(f.inner(orig-f.initial(),orig-f.initial())/f.inner(orig,orig))<1e-13
    h=initial(f,orig,job['variant'])
    first=f.observe(h,0); first.update(extra(f,h))
    acc=np.zeros(5);series=[]; start=0;max_cfl=0.
    substeps=103*(2 if job['half'] else 1);dt=.01/substeps
    if path.exists():
        acc=np.array(old['accum']);series=old['series'];start=old['last_output']+1
        h=np.load(folder/old['field']);max_cfl=old['max_cfl']
    cached=f.rhs(h)
    for oi in range(start,5):
        if oi:
            for sub in range(substeps):
                W0=cached[1][0]
                h,integ,speed=f.step(h,dt,cached)
                if not np.isfinite(h).all(): raise FloatingPointError('Nonfinite field')
                cached=f.rhs(h);integ[0]=.5*dt*(W0+cached[1][0]);acc+=integ
                max_cfl=max(max_cfl,speed*dt/f.dx)
        row=f.observe(h,oi*.01,acc,first);row.update(extra(f,h))
        row['resolved_screen']=(row['width_at_global_peak']['minimum_chord_cells'] or 0)>=6 and row['high_band_energy_fraction']<.001 and row['high_band_enstrophy_fraction']<.01
        fail=max_cfl>.5 or row['divergence_max']*f.dx/max(row['peak_speed'],1e-15)>1e-10 or abs(row['energy_budget_relative_residual'])>.005
        field=f'field-{oi:03d}.npy';np.save(folder/field,h)
        row['field']=field;row['field_sha256']=sha(folder/field);series.append(row)
        result=dict(status='failed_check' if fail else ('complete' if oi==4 else 'running'),job=job,
            n=64,nu=.001,external_force=0,dt=dt,end_time=.04,last_output=oi,field=field,accum=acc.tolist(),
            source_hashes=hashes,series=series,max_cfl=max_cfl)
        save(path,result)
        save(ROOT/'state.json',dict(pid=os.getpid(),active=job['id'],t=oi*.01,end=.04,status=result['status']))
        print(json.dumps(dict(job=job['id'],t=oi*.01,W=row['Wmax'],central=row['interior_central_cylinder_W'],status=result['status'])),flush=True)
        if fail:return result
    return result

def main():
    lock=ROOT/'controller.lock'
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    os.write(fd,str(os.getpid()).encode());os.close(fd)
    try:
        protocol=json.loads((ROOT/'protocol.json').read_text())
        assert sha(ROOT/'source/numerics.py')==protocol['numerics_sha256']
        results=[]
        for variant in ('without-tube','periodic-join'):
            for half in (False,True):
                results.append(run(dict(id=variant+('-half' if half else '-base'),variant=variant,half=half),protocol))
        save(ROOT/'state.json',dict(pid=os.getpid(),status='complete',completed=[r['job']['id'] for r in results if r['status']=='complete'],failed=[r['job']['id'] for r in results if r['status']!='complete']))
    finally:lock.unlink(missing_ok=True)

if __name__=='__main__':main()
