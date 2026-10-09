"""Separate user-authorized methods controls. Never edits the supplied solver."""
import argparse,hashlib,importlib.util,json,os,time
from pathlib import Path
import numpy as np
from diagnostics import Meter,initial,mirror,transform,norm,dot

ROOT=Path(__file__).resolve().parent
def save(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_bytes((json.dumps(v,indent=2,allow_nan=False)+'\n').encode())
def load_solver():
    p=ROOT/'source/clay_hug.py'
    assert hashlib.sha256(p.read_bytes()).hexdigest()=='422b7160f9e96badfc68fb76958c3cb512982905b73b37f8db92e5fd37379c01'
    spec=importlib.util.spec_from_file_location('original_clay_hug',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m
def jobs():
    out=[]
    for n in (33,49,65):
        for sign in (0,1,-1):out.append(dict(n=n,sign=sign,half=False,shape=0,transform='none'))
    for sign in (0,1,-1):out.append(dict(n=65,sign=sign,half=True,shape=0,transform='none'))
    for shape in (1,2):
        for sign in (1,-1):out.append(dict(n=33,sign=sign,half=False,shape=shape,transform='none'))
    for tr in ('shift1','shifthalf','rotate'):
        for sign in (1,-1):out.append(dict(n=33,sign=sign,half=False,shape=0,transform=tr))
    for j in out:j['id']=f"n{j['n']}-s{j['sign']}-shape{j['shape']}-{j['transform']}-{'half' if j['half'] else 'base'}"
    return out
def protocol():
    return dict(L=6,nu=.01,T=.4,jobs=jobs(),forcing='none',method='Unchanged supplied Heun with pressure projection',
        timestep='ceil(.4 / (.04*dx/maxspeed(u0_reference))); dt=.4/steps; controls share reference plus dt. Zero state retains own dt. Half means twice as many steps.',
        output='Every step W, energy, enstrophy and I; full diagnostics and full field at roughly .02 and both endpoints.',
        same_energy='All shape variants have odd unit-L2 phi; use same 0.01*norm(base) amplitude. Base and odd template are orthogonal.',
        translations='1 and 0.5 grid cells along x; Fourier shift all vector components. Diagnostics use conjugated mirror and shifted template, equivalently pull back the whole field.',
        rotation='Proper 90 degree z rotation, rotating positions AND components. Pull back before E,D.',
        one_step='At every grid, both dt and half dt; even and plus starts; mirror, cell shift, half-cell shift and 90 degree rotation.',
        symmetry_tolerance=1e-10,mirror_pair_tolerance=1e-10,
        symmetry_scale='Relative L2 field residual. Conservative roundoff audit ceiling, not a continuum error bound.',
        resolution_flags={'tail_energy_fraction':.001,'tail_enstrophy_fraction':.01},
        curve_tolerance=.01,curve_metric='relative L2 in time against finer curve on common output interval',
        persistent_state='Not tested by this finite .4 batch. A future plateau candidate requires an unforced interval [T/2,T] of length >=1, nonzero |D| >100 times numerical symmetry floor, <=1% variation, matching signed outcome under 3 grids and half dt. Even a pass is a finite-window plateau, not an attractor proof.',
        stop='No additional jobs or longer intervals are created automatically.')
def controls(m):
    out=[]
    for n in (33,49,65):
        b,p,ops,c=initial(m,n);u=b+.01*norm(b)*p
        steps=int(np.ceil(.4/(.04*ops[-1]/np.sqrt(np.sum(u*u,axis=0)).max())))
        for half in (False,True):
            dt=.4/steps/(2 if half else 1)
            for name,v in [('even',b),('plus',u)]:
                step=np.asarray(m.advance(v,dt,ops)); transforms={'mirror':mirror,'shift1':lambda a:transform(a,'shift1'),'shifthalf':lambda a:transform(a,'shifthalf'),'rotate':lambda a:transform(a,'rotate')}
                for kind,fn in transforms.items():
                    other=np.asarray(m.advance(fn(v),dt,ops))
                    residual=norm(other-fn(step))/max(norm(step),1e-300)
                    out.append(dict(n=n,half=half,dt=dt,start=name,transform=kind,relative_step_residual=residual))
    save(ROOT/'one-step-controls.json',out)
def run(m,j):
    dest=ROOT/'runs'/j['id'];rp=dest/'result.json'
    if rp.exists() and json.loads(rp.read_text())['status']=='complete':return
    b,p,ops,c=initial(m,j['n'],j['shape']);u0=b+.01*norm(b)*p
    if j['sign']<0:u0=mirror(u0)
    if j['sign']==0:u0=b.copy()
    tr=j['transform'];u=transform(u0,tr)
    # dt uses standard-shape plus reference for all nonzero starts.
    br,pr,_,_=initial(m,j['n']);ref=br+.01*norm(br)*pr if j['sign'] else br
    steps=int(np.ceil(.4/(.04*ops[-1]/np.sqrt(np.sum(ref*ref,axis=0)).max())))
    if j['half']:steps*=2
    dt=.4/steps;meter=Meter(m,ops,p,tr);dest.mkdir(parents=True,exist_ok=True)
    capture=set(np.rint(np.linspace(0,steps,21)).astype(int));rows=[];dense=[];manifest=[];old=None;I=0;zint=0;prevloc=None
    result=dict(job=j,status='running',dt=dt,steps=steps,L=m.L,nu=m.NU,T=.4,rows=rows,dense=dense,fields=manifest)
    print('START',j['id'],steps,flush=True)
    for step in range(steps+1):
        if step:u=np.asarray(m.advance(u,dt,ops))
        assert np.isfinite(u).all(),j
        basic=meter.basic(u);basic.update(t=step*dt,step=step)
        if old:I+=.5*dt*(old['W']+basic['W']);zint+=.5*dt*(old['enstrophy']+basic['enstrophy'])
        basic['I']=I;basic['viscous_energy_balance_residual']=(basic['energy']-dense[0]['energy']+2*m.NU*zint)/dense[0]['energy'] if dense else 0.
        dense.append(basic);old=basic
        if step in capture:
            r=meter.full(u,prevloc);prevloc=r['peak_index'];r.update(basic)
            r['cfl']=float(np.sqrt(np.sum(u*u,axis=0)).max()*dt/ops[-1]);rows.append(r)
            assert r['divergence_max']<1e-9 and r['cfl']<.5,(j,r['cfl'],r['divergence_max'])
            f=dest/f'field-{step:04d}.npz';np.savez_compressed(f,u=u)
            manifest.append(dict(file=f.name,step=step,t=step*dt,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
            save(rp,result)
    result['status']='complete';save(rp,result)
    print('DONE',j['id'],json.dumps({k:rows[-1][k] for k in ('W','E','D','high_band_enstrophy_fraction')}),flush=True)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--prepare',action='store_true');args=ap.parse_args()
    p=protocol();path=ROOT/'protocol.json'
    if path.exists():assert json.loads(path.read_text())==p,'Protocol changed; do not alter a running study.'
    else:save(path,p)
    if args.prepare:return
    lock=ROOT/'controller.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    state={'pid':os.getpid(),'started_epoch':time.time(),'status':'running','complete':[],'planned':len(jobs())}
    try:
        m=load_solver()
        if not (ROOT/'one-step-controls.json').exists():controls(m)
        for j in jobs():
            state['active']=j['id'];save(ROOT/'state.json',state)
            run(m,j);state['complete'].append(j['id']);save(ROOT/'state.json',state)
        state.update(status='complete',active=None);save(ROOT/'state.json',state)
    except Exception as exc:
        state.update(status='failed',error=repr(exc));save(ROOT/'state.json',state);raise
    finally:lock.unlink()
if __name__=='__main__':main()
