"""Meaningful analytic, conservation, probability and refinement tests."""
import json
from pathlib import Path
import numpy as np
import pytest
from core import (initial_scalar,diffusion,shear,scalar_run,particle_run,column_run,
                  Boussinesq,hazard_model)


def test_analytic_diffusion_variance_and_mass():
    n=64; D=.001; t=8.
    c=initial_scalar(n); got=diffusion(c,D,t)
    expected=.5+(c-.5)*np.exp(-D*(2*np.pi)**2*t)
    assert np.max(abs(got-expected))<1e-13
    assert abs(got.mean()-.5)<1e-14
    assert abs(got.var()/c.var()-np.exp(-2*D*(2*np.pi)**2*t))<1e-13


def test_analytic_shear_advection():
    n=64; x,z=np.meshgrid(np.arange(n)/n,np.arange(n)/n,indexing='ij')
    dt=.123; phase=.45
    expected=.5+.45*np.cos(2*np.pi*(x-dt*np.sin(2*np.pi*z+phase)))
    got=shear(initial_scalar(n),dt,0,phase)
    assert np.max(abs(got-expected))<2e-14
    assert abs(got.var()-initial_scalar(n).var())<1e-14


def test_analytic_translation_diffusion_composition():
    n=64; D=.002; t=.3; c=initial_scalar(n)
    k=2*np.pi*np.fft.fftfreq(n,d=1/n)
    got=diffusion(np.fft.ifft(np.fft.fft(c,axis=0)*np.exp(-1j*k[:,None]*.7*t),axis=0).real,D,t)
    expected=.5+.45*np.cos(2*np.pi*(np.arange(n)/n-.7*t))[:,None]*np.exp(-D*(2*np.pi)**2*t)
    assert np.max(abs(got-expected))<1e-13


def test_scalar_bounds_and_mass():
    r,_=scalar_run(n=64,duration=2.)
    assert np.max(abs(r[:,1]-.5))<1e-12
    assert r[:,5].min()>=0 and r[:,6].max()<=1
    assert np.all(np.diff(r[:,2])<1e-12)


def test_particle_count_and_exact_shear_timestep():
    a,_,xa,ea=particle_run(count=1024,dt=.02,duration=1.)
    b,_,xb,eb=particle_run(count=1024,dt=.01,duration=1.)
    assert np.max(abs(xa[-1]-xb[-1]))<1e-12
    assert np.max(abs(ea-eb))<1e-12
    assert np.all(a[:,-1]==1024) and np.all(b[:,-1]==1024)


def test_brownian_msd():
    D=1e-5; t=.1
    _,_,x,_=particle_run(count=32768,duration=t,flow='still',brownian=D)
    displacement=(x[-1]-x[0]+.5)%1-.5
    msd=np.mean(np.sum(displacement**2,axis=1))
    assert abs(msd/(4*D*t)-1)<.03


def test_stokes_column_budget_and_expected_deposition():
    v,b,z0=column_run(42,count=100000)
    assert b.sum()==100000
    assert abs(b[1]/100000-v*3600/.001)<.004
    _,b2,_=column_run(42,count=10000,adsorption_rate=1/3600)
    assert b2.sum()==10000 and np.all(b2>=0)
    _,b0,_=column_run(42,count=1000,excess_density=0)
    assert b0.tolist()==[1000,0,0]


def test_boussinesq_divergence_mass_and_stable_energy():
    r,_=Boussinesq(32,N2=4.).run(duration=1.)
    assert r[:,5].max()<1e-12
    assert np.max(abs(r[:,1]-.5))<1e-13
    assert np.all(np.diff(r[:,3]+r[:,4])<1e-12)


def test_buoyancy_exchange_cancels_in_energy():
    m=Boussinesq(32,nu=0,kappa=0,diffusivity=0,N2=4)
    h=m.start(42); h[2]=.3*h[1]
    f=np.fft.ifft2(h).real; rhs=np.fft.ifft2(m.rhs(h)).real
    rate=np.mean(f[0]*rhs[0]+f[1]*rhs[1]+f[2]*rhs[2]/4)
    assert abs(rate)<1e-14


def test_hazard_quadrature_and_survival():
    x=np.linspace(0,4,8001); cooling=2.
    rate,H,S=hazard_model(x,bio=1,cooling=cooling)
    numerical=np.concatenate(([0],np.cumsum((rate[:-1]+rate[1:])/2*np.diff(x)/cooling)))
    assert np.max(abs(H-numerical))/H[-1]<2e-6
    assert S[0]==1 and np.all(np.diff(S)<=0) and np.all((S>=0)&(S<=1))


def test_hazard_controls_and_cooling():
    x=np.linspace(0,10,2001)
    _,H,S=hazard_model(x)
    assert np.array_equal(S,hazard_model(x,bio=0,mineral=0)[2])
    assert np.all(hazard_model(x,bio=1)[2]<=S)
    assert np.all(hazard_model(x,cooling=2)[2]>=S)
    assert np.all(hazard_model(x,volume=2)[2]<=S)
    with pytest.raises(ValueError): hazard_model(x,cooling=0)


def test_exponential_survival_monte_carlo():
    rng=np.random.default_rng(81); thresholds=rng.exponential(size=100000)
    for H in [.2,1.,3.]:
        assert abs(np.mean(thresholds>H)-np.exp(-H))<.006


def test_completed_refinement_gates():
    p=Path(__file__).parent/'summary.json'
    if not p.exists(): pytest.skip('run_study.py required for completed-data checks')
    q=json.loads(p.read_text())['numerics']
    assert q['scalar_grid64_to128_L2']<1e-5
    assert q['scalar_grid64_to128_L2']<q['scalar_grid32_to64_L2']
    assert 1.7<q['temporal_observed_order']<2.3
    assert q['scalar_dt01_to005_L2']<1e-4
    assert q['max_solute_mass_error']<1e-12
