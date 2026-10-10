"""Regression checks that reject stable-but-inaccurate oscillation settings."""
import cmath
import numpy as np
import pytest
from oscillation_accuracy import factor, direct_step, measure


@pytest.mark.parametrize("method", ["rk4","implicit_midpoint","gauss_legendre4","exact_rotation"])
def test_direct_step_matches_independent_amplification_formula(method):
    for omega in (2.,376.):
        for h in (.00025,.001,.004):
            y=.3-.7j
            actual=direct_step(method,y,h,omega)
            expected=factor(method,-1j*omega*h)*y
            assert abs(actual-expected)<2e-15


def test_stable_but_damped_rk4_is_rejected():
    result=measure("rk4",.004,10.)
    assert result["numerically_stable"]
    assert result["relative_amplitude_error"]>.99
    assert not result["accuracy_accepted"]


def test_conservative_but_out_of_phase_midpoint_is_rejected():
    result=measure("implicit_midpoint",.00025,10.)
    assert result["relative_amplitude_error"]==0
    assert result["accumulated_phase_error_degrees"]>150
    assert not result["accuracy_accepted"]


@pytest.mark.parametrize("method,h", [("rk4",.000125),("gauss_legendre4",.00025),("implicit_midpoint",.00000390625)])
def test_refined_candidates_pass_both_budgets_over_100_units(method,h):
    result=measure(method,h,100.)
    assert result["accuracy_accepted"]
    assert result["relative_amplitude_error"]<=.001
    assert result["accumulated_phase_error_degrees"]<=1


def test_rk4_tiny_attenuation_is_not_lost_to_roundoff():
    larger=measure("rk4",.0000078125,100.)["relative_amplitude_error"]
    smaller=measure("rk4",.00000390625,100.)["relative_amplitude_error"]
    assert smaller>0
    assert 31.9<larger/smaller<32.1 # global amplitude loss scales as h^5


def test_phase_error_accumulates_without_wrapping_away_lost_cycles():
    one=measure("implicit_midpoint",.0005,1.)["accumulated_phase_error_degrees"]
    hundred=measure("implicit_midpoint",.0005,100.)["accumulated_phase_error_degrees"]
    assert hundred>360
    assert np.isclose(hundred/one,100.,rtol=1e-9)


def test_exact_oracle_is_not_accepted_with_inadequate_output_sampling():
    result=measure("exact_rotation",.008,100.)
    assert result["relative_amplitude_error"]==0
    assert result["accumulated_phase_error_degrees"]==0
    assert not result["phase_sampling_adequate"]
    assert not result["accuracy_accepted"]
