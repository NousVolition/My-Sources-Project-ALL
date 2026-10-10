import numpy as np
from coupled_model import CoupledFlow, run_pair


def test_equal_histogram_and_impulse():
    m=CoupledFlow(48)
    a=m.start(100,'x');b=m.start(100,'y')
    np.testing.assert_allclose(m.physical(a)[:2],m.physical(b)[:2],atol=1e-15)
    np.testing.assert_allclose(np.sort(m.physical(a)[2].ravel()),np.sort(m.physical(b)[2].ravel()),atol=1e-15)
    np.testing.assert_allclose(m.start(100,'x',True)-a,m.start(100,'y',True)-b,atol=1e-14)


def test_viscous_stress_gradient_term():
    m=CoupledFlow(32);k=2*np.pi
    f=np.stack([np.sin(k*m.y),np.zeros_like(m.x),.02+.015*np.cos(k*m.x)])
    h=m.project(np.fft.fft2(f));nu=.01*(1+2.5*f[2]);nux=-.01*2.5*.015*k*np.sin(k*m.x)
    expected=np.stack([-nu*k*k*np.sin(k*m.y),nux*k*np.cos(k*m.y),np.zeros_like(m.x)])
    np.testing.assert_allclose(m.rhs(h)[:2],m.project(np.fft.fft2(expected))[:2],atol=1e-10)


def test_semidiscrete_energy_and_mass():
    m=CoupledFlow(48);h=m.start(102,'y',True);f=m.physical(h);r=m.physical(m.rhs(h))
    assert abs(np.mean(np.sum(f[:2]*r[:2],axis=0))+m.diagnostics(h)['dissipation'])<1e-14
    assert abs(r[2].mean())<1e-14


def test_constant_viscosity_shear_decay():
    m=CoupledFlow(32,feedback=False);k=2*np.pi
    h=m.project(np.fft.fft2(np.stack([np.sin(k*m.y),np.zeros_like(m.x),np.full_like(m.x,.02)])))
    for _ in range(50):h,_=m.step(h,.002)
    expected=np.sin(k*m.y)*np.exp(-.0105*k*k*.1)
    np.testing.assert_allclose(m.physical(h)[0],expected,atol=1e-12)


def test_feedback_off_and_uniform_nulls():
    _,x=run_pair(n=32,duration=.1,organization='x',feedback=False)
    _,y=run_pair(n=32,duration=.1,organization='y',feedback=False)
    np.testing.assert_array_equal(x[:,:,:2],y[:,:,:2])
    _,a=run_pair(n=32,duration=.1,organization='uniform',feedback=True)
    _,b=run_pair(n=32,duration=.1,organization='uniform',feedback=False)
    np.testing.assert_allclose(a,b,atol=1e-13)


def test_initial_grid_polynomial():
    coarse=CoupledFlow(32);fine=CoupledFlow(64)
    np.testing.assert_allclose(coarse.physical(coarse.start(100,'x',True)),fine.physical(fine.start(100,'x',True))[:,::2,::2],atol=1e-15)
