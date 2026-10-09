"""Predefined reduced models, evaluated separately from Navier-Stokes."""
import hashlib,json,os
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parent
def save(name,x):(ROOT/name).write_bytes((json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item())+'\n').encode())
def integrate(fun,y0,t,rtol=1e-10):
    sol=solve_ivp(fun,(float(t[0]),float(t[-1])),y0,t_eval=t,method='DOP853',rtol=rtol,atol=rtol*.01)
    assert sol.success and np.isfinite(sol.y).all()
    return sol.y.T
def crossings(x):return int(np.sum(x[1:]*x[:-1]<0))
def main():
    protocol=dict(status='fixed before runs',ball=dict(r=[-1,1],damping=.4,amplitudes=[0,.15,.5],period=4,end=80,initial=[.1,0],mirror='x,v,force -> negatives',method='DOP853 rtol1e-10 atol1e-12'),
        four_node=dict(r=.5,k_AB=.6,k_CD=.4,k_AC_and_BD=.4,k_AD_and_BC=.1,deltas=[0,.05,.2],end=40,initial=[.1,.1,.2,.2],mirror='Swap A/B and C/D, conjugate full coupling matrix and source presence. Reciprocal edges used.'),
        listening=dict(r=-1,k=.3,gamma=.2,beta=.8,g=1,h=[0,.1,.4],initial_bias=[0,.0001,-.0001],end=80,mirror='A/B swap, q -> -q',feedback='h*q enters A, -h*q enters B; this is a reduced ODE, not a new fluid force.'),
        spatial=dict(L=float(2*np.pi),r=.5,kappa=.05,period=4,end=32,n=[128,256],dt=.005,half_dt=.0025,amplitudes=[0,.15,.5],mirror='M[a](x)=-a(-x), forcing sin(x) invariant under M',method='Strang exact spectral diffusion + RK4 reaction/forcing, strict |mode|<N/4 for cubic dealiasing',initial='.05+0.2*sin(x); optional .1*exp(4*(cos(x-.8)-1))'),
        interpretation='Finite-window diagnostics and model tests. No inferred fluid pitchfork, attractor, chaos, or coefficient equivalence.')
    save('protocol.json',protocol)
    out={'ball':{},'four_node':{},'listening':{},'spatial':{},'checks':{}}
    t=np.linspace(0,80,1601)
    for r in (-1,1):
      for A in (0,.15,.5):
        for sign in (1,-1):
          fn=lambda tt,y: [y[1],r*y[0]-y[0]**3-.4*y[1]+sign*A*np.sin(2*np.pi*tt/4)]
          y=integrate(fn,[sign*.1,0],t);key=f'r{r}-A{A}-s{sign}'
          st=y[::80];late=y[t>=60,0]
          out['ball'][key]=dict(r=r,A=A,sign=sign,t=t.tolist(),x=y[:,0].tolist(),v=y[:,1].tolist(),strobe=st.tolist(),
            crossings=crossings(y[:,0]),late_crossings=crossings(late),last_x=float(y[-1,0]),last_five_strobe_range=float(np.ptp(st[-5:,0])))
    out['checks']['ball_mirror_max']=max(float(np.max(abs(np.array(out['ball'][f'r{r}-A{a}-s1']['x'])+np.array(out['ball'][f'r{r}-A{a}-s-1']['x'])))) for r in (-1,1) for a in (0,.15,.5))
    check=integrate(lambda tt,y:[y[1],y[0]-y[0]**3-.4*y[1]+.5*np.sin(2*np.pi*tt/4)],[.1,0],t,1e-12)
    out['checks']['ball_tolerance_refinement_max']=float(np.max(abs(check[:,0]-out['ball']['r1-A0.5-s1']['x'])))
    print('BALL complete',flush=True)
    t=np.linspace(0,40,801);perm=np.array([1,0,3,2])
    for delta in (0,.05,.2):
      for absent in (False,True) if delta==0 else (False,):
        K=np.array([[0,.6,.4+delta,.1],[.6,0,.1,.4-delta],[.4+delta,.1,0,.4],[.1,.4-delta,.4,0.]])
        present=np.ones(4);y0=np.array([.1,.1,.2,.2])
        if absent:K[3,:]=0;K[:,3]=0;present[3]=0;y0[3]=0
        for sign in (1,-1):
          kk=K if sign==1 else K[np.ix_(perm,perm)];pr=present if sign==1 else present[perm];yy=y0 if sign==1 else y0[perm]
          y=integrate(lambda tt,z:pr*(.5*z-z**3+kk@z-kk.sum(axis=1)*z),yy,t)
          qab=(y[:,0]-y[:,1])/2;qcd=(y[:,2]-y[:,3])/2
          out['four_node'][f'd{delta}-absent{absent}-s{sign}']=dict(delta=delta,missing_D=absent,sign=sign,K=kk.tolist(),t=t.tolist(),states=y.tolist(),q_AB=qab.tolist(),q_CD=qcd.tolist(),E_model=np.hypot(qab,qcd).tolist())
    out['checks']['four_node_mirror_max']=max(float(np.max(abs(np.array(v['states'])[:,perm]-np.array(out['four_node'][key.replace('-s1','-s-1')]['states'])))) for key,v in out['four_node'].items() if key.endswith('-s1'))
    print('FOUR NODE complete',flush=True)
    t=np.linspace(0,80,1601)
    for h in (0,.1,.4):
      for bias in (0,.0001,-.0001):
        def fn(tt,z):
          A,B,q=z
          return [-A-A**3+.3*(B-A)+h*q,-B-B**3+.3*(A-B)-h*q,-.2*q+.8*(A-B)-q**3]
        y=integrate(fn,[bias,-bias,0],t)
        out['listening'][f'h{h}-b{bias}']=dict(h=h,bias=bias,t=t.tolist(),states=y.tolist(),q_AB=((y[:,0]-y[:,1])/2).tolist())
    out['checks']['listening_zero_control_max']=max(float(np.max(abs(np.array(out['listening'][f'h{h}-b0']['states'])))) for h in (0,.1,.4))
    out['checks']['listening_mirror_max']=max(float(np.max(abs(np.array(out['listening'][f'h{h}-b0.0001']['states'])+np.array(out['listening'][f'h{h}-b-0.0001']['states'])))) for h in (0,.1,.4))
    print('LISTENING complete',flush=True)
    for n,half in ((128,False),(256,False),(256,True)):
      x=2*np.pi*np.arange(n)/n;mode=np.fft.fftfreq(n)*n;keep=abs(mode)<n/4;ref=(-np.arange(n))%n
      for A,local in (((.5,False),) if half else ((0,False),(.15,False),(.5,False),(.15,True))):
       for sign in ((1,) if half else (1,-1)):
        dt=.0025 if half else .005;steps=int(round(32/dt))
        base=.05+.2*np.sin(x)
        if local:base+=.1*np.exp(4*(np.cos(x-.8)-1))
        a=base if sign==1 else -base[ref]
        a=np.fft.ifft(np.fft.fft(a)*keep).real
        diffusion=np.exp(-.05*mode**2*dt/2);rows=[];images=[]
        def reaction(z,tt):
          rhs=.5*z-z**3+A*np.sin(2*np.pi*tt/4)*np.sin(x)
          return np.fft.ifft(np.fft.fft(rhs)*keep).real
        for step in range(steps+1):
          tt=step*dt
          if step%40==0:
            antis=a+a[ref];rows.append(dict(t=tt,q=float(a.mean()),unsigned=float(np.sqrt(np.mean(antis**2))),rms=float(np.sqrt(np.mean(a*a)))))
          if step%int(round(4/dt))==0:images.append(a.tolist())
          if step==steps:break
          a=np.fft.ifft(np.fft.fft(a)*diffusion).real
          k1=reaction(a,tt);k2=reaction(a+dt*k1/2,tt+dt/2);k3=reaction(a+dt*k2/2,tt+dt/2);k4=reaction(a+dt*k3,tt+dt)
          a=a+dt*(k1+2*k2+2*k3+k4)/6
          a=np.fft.ifft(np.fft.fft(a)*diffusion*keep).real
        assert np.isfinite(a).all()
        key=f'n{n}-A{A}-loc{local}-s{sign}'+('-half' if half else '')
        out['spatial'][key]=dict(n=n,A=A,local=local,sign=sign,dt=dt,rows=rows,snapshot_times=list(range(0,33,4)),snapshots=images,final=a.tolist())
        print('SPATIAL',n,A,local,sign,flush=True)
    out['checks']['spatial_mirror_max']=max(float(np.max(abs(np.array(v['final'])+np.array(out['spatial'][key.replace('-s1','-s-1')]['final'])[(-np.arange(v['n']))%v['n']]))) for key,v in out['spatial'].items() if key.endswith('-s1'))
    out['checks']['spatial_grid_relative_errors']={}
    for A,local in ((0,False),(.15,False),(.5,False),(.15,True)):
      a=np.array(out['spatial'][f'n128-A{A}-loc{local}-s1']['final']);b=np.array(out['spatial'][f'n256-A{A}-loc{local}-s1']['final'])[::2]
      out['checks']['spatial_grid_relative_errors'][f'A{A}-loc{local}']=float(np.linalg.norm(a-b)/np.linalg.norm(b))
    a=np.array(out['spatial']['n256-A0.5-locFalse-s1']['final']);b=np.array(out['spatial']['n256-A0.5-locFalse-s1-half']['final'])
    out['checks']['spatial_dt_relative_error']=float(np.linalg.norm(a-b)/np.linalg.norm(b))
    out['status']='complete';out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    save('reduced-results.json',out);print('REDUCED COMPLETE',json.dumps(out['checks']),flush=True)
if __name__=='__main__':main()
