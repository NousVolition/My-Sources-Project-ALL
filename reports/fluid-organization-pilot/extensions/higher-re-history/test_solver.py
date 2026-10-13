import importlib.util
from pathlib import Path
import numpy as np
from solver import Flow,initial,impulse,resize,L,observe

def test_exact_energy_injection():
 f=Flow(24);h=initial(f,2);assert abs(f.inner(h,f.force(h))-f.power)<1e-14
def test_advection_energy_and_divergence():
 f=Flow(24,0,0);h=initial(f,2)[None];r,*_=f.parts(h)
 assert abs(f.inner(h,r)[0])<1e-13
 assert abs(f.real(1j*sum(k*a for k,a in zip(f.k,h[0])))).max()<1e-13
 assert not np.any(h[0,:,~f.keep])
def test_difference_budget_includes_forcing():
 f=Flow(24);h=initial(f,4);d=impulse(f,.75);hs=np.stack((h,h+d));r,b,q,_=f.parts(hs)
 grad=np.stack([f.real(1j*k*h) for k in f.k],axis=1);v=f.real(d)
 p=-np.mean(np.einsum('iabc,jabc,ijabc->abc',v,v,grad))
 assert abs(q[0,0]-p)<1e-13
 assert abs(f.inner(d,r[1]-r[0])-(q[0,0]-q[0,1]+q[0,2]))<1e-13
 assert abs(q[0,2])>1e-4
def test_identical_clones():
 f=Flow(24);h=initial(f,9);out,_,q,_=f.step(np.stack((h,h.copy())),.01)
 assert np.array_equal(out[0],out[1]);assert np.array_equal(q,np.zeros((1,3)))
def test_exact_beltrami_solution():
 f=Flow(24,.2,0);a=np.arange(f.n)*L/f.n;z=np.broadcast_to(a[None,None,:],(f.n,)*3)
 u=np.stack((np.sin(2*np.pi*z/L),np.cos(2*np.pi*z/L),np.zeros_like(z)));h=f.hat(u)[None];start=h.copy()
 for _ in range(10):h,*_=f.step(h,.02)
 exact=start*np.exp(-.2*(2*np.pi/L)**2*.2)
 assert np.sqrt(f.inner(h-exact,h-exact)/f.inner(exact,exact))[0]<1e-9
def test_integrated_budgets():
 f=Flow(24);h=initial(f,4);hs=np.stack((h,h+impulse(f,.75)));e0=.5*f.inner(hs,hs);d0=.5*f.inner(hs[1]-hs[0],hs[1]-hs[0]);b=0;q=0
 for _ in range(10):hs,db,dq,_=f.step(hs,.005);b+=db;q+=dq
 assert np.max(abs(.5*f.inner(hs,hs)-e0-b[:,0]+b[:,1]))<1e-8
 assert abs(.5*f.inner(hs[1]-hs[0],hs[1]-hs[0])-d0-q[0,0]+q[0,1]-q[0,2])<1e-8
def test_same_initial_polynomial_and_impulse():
 a=Flow(24);b=Flow(48)
 for x,y in ((initial(a,4),initial(b,4)),(impulse(a,.75),impulse(b,.75))):
  assert np.allclose(a.real(x),b.real(y)[:,::2,::2,::2],atol=1e-13)
  assert np.allclose(resize(x,48),y,atol=1e-10)
def test_local_gradient():
 f=Flow(24);h=initial(f,2);v,_,_=observe(f,h)
 g=np.stack([f.real(1j*k*h)[:,:,12,12] for k in f.k],axis=1)[:,:,12]
 assert np.allclose(v[11:20].reshape(3,3),g,atol=1e-12)
def test_optimized_matches_original():
 spec=importlib.util.spec_from_file_location('old',Path(__file__).with_name('screen_original_solver.py'));old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
 f=Flow(24);g=old.Flow(24);h=initial(f,4);hs=np.stack((h,h+impulse(f,.75)))
 for a,b in zip(f.step(hs,.01),g.step(hs,.01)):assert np.allclose(a,b,rtol=1e-12,atol=1e-12)

def test_forced_beltrami_equilibrium_has_no_cascade():
 f=Flow(24);a=np.arange(f.n)*L/f.n;z=np.broadcast_to(a[None,None,:],(f.n,)*3);k=2*np.pi/L;amplitude=np.sqrt(f.power/(f.nu*k*k))
 u=amplitude*np.stack((np.sin(k*z),np.cos(k*z),np.zeros_like(z)));h=f.project(f.hat(u));rhs,*_=f.parts(h[None]);_,s,_=observe(f,h)
 assert np.sqrt(f.inner(rhs,rhs))[0]<1e-11;assert abs(s['flux_k4'])<1e-12
 assert abs(s['dissipation']-f.power)<1e-12
