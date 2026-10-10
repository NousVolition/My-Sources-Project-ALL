"""Independent PDE, conservation, boundary, reference and integration checks."""
from pathlib import Path
import sys
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

sys.path.insert(0,str(Path(__file__).resolve().parent))
from convection import Layer, RAC, KC, Q, modal_matrix, run


def test_critical_point_from_neutral_curve_minimization():
    optimum=minimize_scalar(lambda k:(k*k+np.pi**2)**3/(k*k),bounds=(.2,8),method='bounded',options={'xatol':1e-12})
    assert abs(optimum.x-KC)<1e-7
    assert abs(optimum.fun-RAC)<1e-9


@pytest.mark.parametrize('pr',[.7,7.])
@pytest.mark.parametrize('ratio',[.8,1.,1.2])
def test_full_pde_linear_mode_matches_separate_matrix(pr,ratio):
    layer=Layer(24,ratio*RAC,pr=pr,nonlinear=False)
    h=layer.initial(amplitude=1e-4,mixed=False)
    r,_=layer.rhs(h,diffusion=True)
    expected=modal_matrix(ratio*RAC,pr)@h[[1,2],1,1]
    # At the neutral mode large diffusion and buoyancy terms cancel. Scale
    # roundoff by those terms, not by the nearly zero net derivative.
    term_scale=np.max(abs(modal_matrix(ratio*RAC,pr))@abs(h[[1,2],1,1]))
    np.testing.assert_allclose(r[[1,2],1,1],expected,rtol=2e-13,atol=32*np.finfo(float).eps*term_scale)


def test_nonlinear_energy_and_temperature_budget_identity():
    layer=Layer(24,1100.,pr=3.)
    h=layer.initial(amplitude=.02)
    f=np.fft.ifft2(h).real
    r,_=layer.rhs(h,diffusion=True)
    r=np.fft.ifft2(r).real
    dx=np.fft.ifft2(1j*layer.kx*h).real
    dz=np.fft.ifft2(1j*layer.kz*h).real
    expected=3*1100*np.mean(f[1]*f[2])-3*np.mean(np.sum(dx[:2]**2+dz[:2]**2,axis=0))
    observed=np.mean(f[0]*r[0]+f[1]*r[1])
    assert abs(observed-expected)<1e-12
    expected_heat=np.mean(f[1]*f[2])-np.mean(dx[2]**2+dz[2]**2)
    assert abs(np.mean(f[2]*r[2])-expected_heat)<1e-14


def test_wall_conditions_and_incompressibility_survive_nonlinear_steps():
    layer=Layer(24,1.2*RAC)
    h=layer.initial(amplitude=.02)
    for _ in range(40):
        h,_,_=layer.step(h,.0025)
    _,speed,div,wall,_=layer.diagnostics(h)
    assert div/speed<1e-12
    assert wall/speed<1e-12


def test_forbidden_wall_modes_are_removed_not_just_divergence_projected():
    layer=Layer(24,1.2*RAC)
    h=layer.initial()
    h[1,1,0]=1e-20  # Divergence-free, but violates impermeable plates.
    h[2,1,0]=1e-20  # Violates fixed plate temperature.
    h=layer.project(h)
    assert h[1,1,0]==0 and h[2,1,0]==0
    for _ in range(125):
        h,_,_=layer.step(h,.01)
    _,speed,div,wall,_=layer.diagnostics(h)
    assert wall/speed<1e-12
    assert div/speed<1e-12


@pytest.mark.parametrize('ratio',[-4.,.95,1.05])
def test_matrix_exponential_against_independent_adaptive_integration(ratio):
    matrix=modal_matrix(ratio*RAC)
    y0=np.array([Q*ratio,1.])
    t=np.linspace(0,.8,101)
    reference=solve_ivp(lambda t,y:matrix@y,(0,.8),y0,t_eval=t,method='DOP853',rtol=2e-12,atol=1e-13)
    exact=np.array([expm(matrix*ti)@y0 for ti in t])
    np.testing.assert_allclose(reference.y.T,exact,rtol=2e-9,atol=1e-10)


def test_linear_pde_trajectory_against_exact_mode_solution():
    layer=Layer(16,1.05*RAC,nonlinear=False)
    h=layer.initial(amplitude=1e-7,mixed=False)
    initial=h[[1,2],1,1].copy()
    for _ in range(100):
        h,_,_=layer.step(h,.001)
    expected=expm(modal_matrix(1.05*RAC)*.1)@initial
    np.testing.assert_allclose(h[[1,2],1,1],expected,rtol=1e-8,atol=1e-17)


def test_zero_disturbance_and_buoyancy_ablation():
    result,rows,_,_=run(n=16,duration=.2,amplitude=0,ratio=1.2)
    assert not np.any(rows[:,1:])
    layer=Layer(16,1.2*RAC,buoyancy=False,nonlinear=False)
    h=layer.initial(amplitude=1e-7,mixed=False)
    initial=h[1,1,1]
    for _ in range(100):
        h,_,_=layer.step(h,.001)
    np.testing.assert_allclose(h[1,1,1],initial*np.exp(-Q*.1),rtol=1e-12)
