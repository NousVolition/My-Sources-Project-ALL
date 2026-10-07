"""Independent Fourier modes and real-space differences check the diagnostic."""
import itertools

import numpy as np
import pytest

from hug_resolution import grid_spectrum


def test_constant_field_has_only_dc_energy():
    u = np.ones((3,8,8,8))
    got, spectrum = grid_spectrum(u,box=4.)
    assert got['energy'] == pytest.approx(96.)
    for key in ('checkerboard_energy_fraction','nyquist_plane_energy_fraction',
                'short_wave_energy_fraction','centered_gradient_rms','neighbor_gradient_rms'):
        assert got[key] == 0.
    assert spectrum['energy_fraction'][0] == 1.
    assert sum(spectrum['neighbor_gradient_fraction']) == 0.


@pytest.mark.parametrize('axes',[a for a in itertools.product((0,1),repeat=3) if any(a)])
def test_all_seven_joint_checkerboards_are_detected(axes):
    n, box = 8, 4.
    x = np.indices((n,n,n))
    u = np.zeros((3,n,n,n))
    u[0] = (-1.)**sum(a*b for a,b in zip(axes,x))
    got, _ = grid_spectrum(u,box)
    assert got['checkerboard_energy_fraction'] == pytest.approx(1.)
    assert got['nyquist_plane_energy_fraction'] == pytest.approx(1.)
    assert got['short_wave_energy_fraction'] == pytest.approx(1.)
    assert got['centered_gradient_rms'] == 0.
    assert got['neighbor_gradient_rms'] == pytest.approx(2*np.sqrt(sum(axes))/(box/n))
    assert got['gradient_rms_gap_percent'] == pytest.approx(100.)


@pytest.mark.parametrize('axis',(0,1,2))
def test_sinusoid_has_known_energy_and_derivative_symbols(axis):
    n, mode, box = 16, 3, 6.
    u = np.zeros((3,n,n,n))
    u[axis] = np.sin(2*np.pi*mode*np.indices((n,n,n))[axis]/n)
    got, spectrum = grid_spectrum(u,box)
    assert got['energy'] == pytest.approx(box**3/4)
    assert got['centered_gradient_rms'] == pytest.approx(np.sin(2*np.pi*mode/n)/(box/n)/np.sqrt(2))
    assert got['neighbor_gradient_rms'] == pytest.approx(2*np.sin(np.pi*mode/n)/(box/n)/np.sqrt(2))
    assert spectrum['energy_fraction'][mode] == pytest.approx(1.)
    assert got['short_wave_energy_fraction'] < 1e-28


def test_near_nyquist_wave_is_not_mislabeled_as_joint_blind_mode():
    n, mode = 16, 7
    u = np.zeros((3,n,n,n))
    u[0] = np.sin(2*np.pi*mode*np.indices((n,n,n))[2]/n)
    got, _ = grid_spectrum(u)
    assert got['short_wave_energy_fraction'] == pytest.approx(1.)
    assert got['checkerboard_energy_fraction'] < 1e-28
    assert got['nyquist_plane_energy_fraction'] < 1e-28
    assert got['gradient_rms_gap_percent'] == pytest.approx(100*(1-np.cos(np.pi*mode/n)))


def test_nyquist_plane_can_have_variation_in_another_direction():
    n = 16
    x,y,z = np.indices((n,n,n))
    u = np.zeros((3,n,n,n))
    u[0] = (-1.)**x*np.sin(2*np.pi*y/n)
    got, _ = grid_spectrum(u)
    assert got['nyquist_plane_energy_fraction'] == pytest.approx(1.)
    assert got['checkerboard_energy_fraction'] < 1e-28
    assert got['centered_gradient_rms'] > 0.


def test_random_field_matches_independent_real_space_differences_and_translation():
    u = np.random.default_rng(24).normal(size=(3,8,8,8))
    before = u.copy()
    dx = 6/8
    center = sum(np.mean(((np.roll(c,-1,d)-np.roll(c,1,d))/(2*dx))**2)
                 for c in u for d in range(3))
    neighbor = sum(np.mean(((np.roll(c,-1,d)-c)/dx)**2) for c in u for d in range(3))
    got, spectrum = grid_spectrum(u)
    shifted, shifted_spectrum = grid_spectrum(np.roll(u,(1,2,3),axis=(1,2,3)))
    assert got['energy'] == pytest.approx(.5*dx**3*np.sum(u*u))
    assert got['centered_gradient_rms'] == pytest.approx(np.sqrt(center))
    assert got['neighbor_gradient_rms'] == pytest.approx(np.sqrt(neighbor))
    assert got == pytest.approx(shifted)
    np.testing.assert_allclose(spectrum['energy_fraction'],shifted_spectrum['energy_fraction'],rtol=1e-13)
    np.testing.assert_array_equal(u,before)
    assert sum(spectrum['energy_fraction']) == pytest.approx(1.)
    assert sum(spectrum['neighbor_gradient_fraction']) == pytest.approx(1.)


def test_invalid_fields_are_rejected():
    for u,box in [(np.zeros((3,8,8,8)),6.), (np.ones((3,7,7,7)),6.),
                  (np.ones((3,8,8,8)),0.), (np.ones((2,8,8,8)),6.),
                  (np.full((3,8,8,8),np.nan),6.), (np.ones((3,8,8,8),complex),6.)]:
        with pytest.raises(ValueError):
            grid_spectrum(u,box)
