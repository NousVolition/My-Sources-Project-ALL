"""Independent identities for the new periodic evolution experiment."""
import numpy as np
import pytest

from hidden_flow_experiment import PeriodicFlow, relative_field_error, run


@pytest.mark.parametrize("name", ("A", "B", "C", "C3D"))
def test_initial_conditions_share_energy_and_coarse_view(name):
    sim = PeriodicFlow(18)
    h = sim.initial(name)
    assert sim.norm2(h)/2 == pytest.approx(1.5, abs=1e-13)
    assert sim.norm2(h*sim.low) < 1e-27
    assert sim.norm2(np.sum(sim.k*h, axis=0)) < 1e-27


def test_c_initial_acceleration_has_the_predicted_low_mode():
    sim = PeriodicFlow(18)
    rhs = sim.real(sim.rhs(sim.initial("C"))*sim.low)
    _, y, _ = sim.xyz
    assert np.max(np.abs(rhs[:2])) < 1e-13
    assert np.max(np.abs(rhs[2]+np.sin(y))) < 1e-13


def test_shear_evolution_matches_exact_heat_solution():
    report, got = run("A", 12, .02, end=.4, viscosity=.1)
    sim = PeriodicFlow(12, .1)
    exact = np.exp(-4*.1*.4)*sim.initial("A")
    assert relative_field_error(got, exact) < 1e-10
    assert report["summary"]["max_abs_energy_balance_error"] < 1e-10


def test_c_horizontal_flow_matches_the_exact_decaying_shear():
    _, got = run("C", 18, .01, end=.4)
    sim = PeriodicFlow(18)
    exact = np.exp(-5*sim.nu*.4)*sim.initial("C")
    assert np.max(np.abs(got[:2]-exact[:2])) < 1e-12
    assert sim.norm2(1j*sim.k[2]*got) < 1e-27


def test_three_dimensional_rhs_preserves_divergence_and_energy_law():
    sim = PeriodicFlow(18)
    h = sim.initial("C3D")
    rhs = sim.rhs(h)
    assert sim.norm2(1j*sim.k[2]*h) > .1
    assert sim.norm2(np.sum(sim.k*rhs, axis=0)) < 1e-25
    energy_rate = float(np.real(np.sum(np.conj(h)*rhs*sim.weights)))
    assert energy_rate == pytest.approx(-sim.dissipation(h), abs=1e-12)


def test_projection_removes_a_pure_pressure_gradient():
    sim = PeriodicFlow(18)
    x, y, z = sim.xyz
    scalar = np.cos(2*x+y)+np.sin(2*z+2*y)
    gradient = 1j*sim.k*sim.fft(scalar)
    assert sim.norm2(sim.project(gradient)) < 1e-25


def test_comparison_counts_modes_omitted_by_the_coarser_grid():
    coarse, fine = PeriodicFlow(18), PeriodicFlow(30)
    c, f = coarse.initial("C"), fine.initial("C")
    assert relative_field_error(c, f) < 1e-14
    f[2, 7, 0, 0] = .1
    f[2, -7, 0, 0] = .1
    assert relative_field_error(c, f) == pytest.approx(np.sqrt(.02/3.02))
