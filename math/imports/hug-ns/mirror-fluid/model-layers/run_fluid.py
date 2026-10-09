"""Separate driven fluid controls for the reduced-model comparison."""
import hashlib,importlib.util,json,os,sys,time
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source'))
from diagnostics import Meter,initial,mirror,norm,dot
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.tmp');tmp.write_bytes((json.dumps(x,indent=2,allow_nan=False)+'\n').encode());os.replace(tmp,p)
def solver():
    p=ROOT/'source/clay_hug.py';assert hashlib.sha256(p.read_bytes()).hexdigest()=='422b7160f9e96badfc68fb76958c3cb512982905b73b37f8db92e5fd37379c01'
    spec=importlib.util.spec_from_file_location('clay_hug_layers',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def jobs():
    js=[]
    for parity in ('even-gap','odd-lean'):
      for A,P in ((.15,.2),(.5,.2),(.5,.1)):
       for sign in (1,-1):js.append(dict(n=49,parity=parity,A=A,P=P,sign=sign,half=False))
      js += [dict(n=65,parity=parity,A=.5,P=.2,sign=1,half=False),dict(n=49,parity=parity,A=.5,P=.2,sign=1,half=True)]
    js.append(dict(n=49,parity='even-gap',A=.5,P=.2,sign=0,half=False))
    for j in js:j['id']=f"n{j['n']}-{j['parity']}-A{j['A']}-P{j['P']}-s{j['sign']}-{'half' if j['half'] else 'base'}"
    return js
def driven_step(m,u,dt,t,ops,F,P):
    force1=np.sin(2*np.pi*t/P)*F;force2=np.sin(2*np.pi*(t+dt)/P)*F
    r1=np.asarray(m.rhs(u,*ops[:5]))+force1;y=u+dt*r1
    r2=np.asarray(m.rhs(y,*ops[:5]))+force2
    out=np.asarray(m.project(u+.5*dt*(r1+r2),*ops[:5]))
    return out,dt*.5*m.L**3*(dot(u,force1)+dot(y,force2))
def line_gap(w,dx):
    n=w.shape[1];mag=np.sqrt(np.sum(w*w,axis=0));x=np.linspace(-3,3,4*n,endpoint=False)
    line=map_coordinates(mag,np.array([(x+3)/dx,np.full(len(x),n/2),np.full(len(x),n/2)]),order=1,mode='grid-wrap')
    left=np.where(x<0)[0];right=np.where(x>=0)[0];li=left[np.argmax(line[left])];ri=right[np.argmax(line[right])]
    lc=np.where(line[li:ri+1]<line[li]*.5)[0];rc=np.where(line[li:ri+1]<line[ri]*.5)[0]
    if not len(lc) or not len(rc):return 0.
    return float(max(0,x[li+rc[-1]]-x[li+lc[0]]))
def check(m):
    b,p,ops,c=initial(m,17);u=b+.01*norm(b)*p;dt=.0001
    v,_=driven_step(m,u,dt,0.,ops,np.zeros_like(u),.2);ref=np.asarray(m.advance(u,dt,ops))
    error=norm(v-ref)/norm(ref);assert error<1e-12
    # Reflection must transform the whole force, not only the starting state.
    F=p*.1
    v,_=driven_step(m,u,dt,.013,ops,F,.2);vm,_=driven_step(m,mirror(u),dt,.013,ops,mirror(F),.2)
    equiv=norm(vm-mirror(v))/norm(v);assert equiv<1e-12
    save(ROOT/'implementation-check.json',dict(no_force_step_relative=error,forced_mirror_step_relative=equiv,solver_sha256=hashlib.sha256((ROOT/'source/clay_hug.py').read_bytes()).hexdigest()))
def run(m,j):
    folder=ROOT/'fluid-runs'/j['id'];rp=folder/'result.json'
    if rp.exists() and json.loads(rp.read_text())['status']=='complete':return
    b,p,ops,c=initial(m,j['n']);U0=norm(b);u=b+j['sign']*.01*U0*p;meter=Meter(m,ops,p)
    x,y,z,*_=m.grids(j['n'])
    breath=np.exp(-8*(x*x+y*y));raw=np.array([breath*np.sign(x),np.zeros_like(breath),np.zeros_like(breath)])
    gap=np.asarray(m.project(raw,*ops[:5]));gap=(gap+mirror(gap))*.5;gap-=gap.mean(axis=(1,2,3),keepdims=True);gap/=norm(gap)
    force_shape=gap if j['parity']=='even-gap' else p*(1 if j['sign']>=0 else -1)
    force_shape-=force_shape.mean(axis=(1,2,3),keepdims=True);force_shape=np.asarray(m.project(force_shape,*ops[:5]));force_shape/=norm(force_shape)
    F=j['A']*U0/.4*force_shape
    ref=b+.01*U0*p;maxspeed=float(np.sqrt(np.sum(ref*ref,axis=0)).max())
    steps=int(np.ceil(.4/min(.04*ops[-1]/maxspeed,j['P']/80)))
    if j['half']:steps*=2
    dt=.4/steps;capture=set(np.rint(np.linspace(0,steps,21)).astype(int));rows=[];dense=[];fields=[];start=0;I=0.;viscint=0.;workint=0.
    result=dict(job=j,status='running',L=6,nu=.01,T=.4,dt=dt,steps=steps,initial_norm=U0,forcing='F(t)=A*U0/0.4*sin(2*pi*t/P)*unit_template; included at both Heun stages',force_parity_residual=norm(mirror(F)-(F if j['parity']=='even-gap' else -F))/norm(F),rows=rows,dense=dense,fields=fields)
    folder.mkdir(parents=True,exist_ok=True)
    if rp.exists():
      result=json.loads(rp.read_text());assert result['job']==j
      rows=result['rows'];dense=result['dense'];fields=result['fields'];u=np.load(folder/fields[-1]['file'])['u'];start=dense[-1]['step']+1;I=dense[-1]['I'];viscint=dense[-1]['integrated_viscous_loss'];workint=dense[-1]['integrated_force_work']
    print('START',j['id'],steps,flush=True)
    for step in range(start,steps+1):
      if step:
        prev=dense[-1];u,work=driven_step(m,u,dt,(step-1)*dt,ops,F,j['P']);workint+=work
      assert np.isfinite(u).all()
      obs=meter.basic(u);a=u-mirror(u);obs.update(t=step*dt,step=step,D=dot(a,p),E=norm(a)/norm(u))
      if step:I+=dt*.5*(prev['W']+obs['W']);viscint+=dt*m.NU*(prev['enstrophy']+obs['enstrophy'])
      obs.update(I=I,integrated_viscous_loss=viscint,integrated_force_work=workint,energy_balance=(obs['energy']-(dense[0]['energy'] if dense else obs['energy'])+viscint-workint)/(dense[0]['energy'] if dense else obs['energy']))
      dense.append(obs)
      if step in capture:
        full=meter.full(u);full.update(obs);full['cfl']=float(np.max(np.sqrt(np.sum(u*u,axis=0)))*dt/ops[-1]);full['gap_center_line_half_peak']=line_gap(meter.real(meter.curlh(meter.fft(u))),ops[-1])
        assert full['divergence_max']<1e-9 and full['cfl']<.5 and abs(full['energy_balance'])<.005
        rows.append(full);path=folder/f'field-{step:04d}.npz';np.savez_compressed(path,u=u);fields.append(dict(file=path.name,step=step,t=step*dt,sha256=hashlib.sha256(path.read_bytes()).hexdigest()));save(rp,result)
    result['status']='complete';save(rp,result);print('DONE',j['id'],json.dumps({k:rows[-1][k] for k in ('E','D','W','energy_balance')}),flush=True)
def main():
    protocol=dict(jobs=jobs(),T=.4,nu=.01,L=6,base='Same methods-control mirrored start, plus/minus one percent unit antisymmetric template',drive='Separate driven experiment. Use actual dt in both Heun stages; no per-step unscaled kick.',force_parities='Even gap preserves reflection. Odd lean drive reverses in the paired mirrored test.',reference='Unforced controls reused from methods-controls, not rerun.',scope='17 runs only; no automatic extensions',comparison='49 versus65 and base versus half dt for strong P=.2 positive starts; whole-field mirrored pairs at49; even zero-bias control',interpretation='A scalar signed additive drive and a mirror-even gap drive are distinct mechanisms.')
    path=ROOT/'fluid-protocol.json'
    if path.exists():assert json.loads(path.read_text())==protocol
    else:save(path,protocol)
    lock=ROOT/'fluid-controller.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    state=dict(pid=os.getpid(),status='running',complete=[],planned=len(jobs()),active=None);save(ROOT/'fluid-state.json',state)
    try:
      m=solver();check(m)
      for j in jobs():
        state['active']=j['id'];save(ROOT/'fluid-state.json',state);run(m,j);state['complete'].append(j['id']);save(ROOT/'fluid-state.json',state)
      state.update(status='complete',active=None);save(ROOT/'fluid-state.json',state)
    except Exception as exc:
      state.update(status='failed',error=repr(exc));save(ROOT/'fluid-state.json',state);raise
    finally:lock.unlink()
if __name__=='__main__':main()
