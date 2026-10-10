"""Sequence interpretation, equal work and independent amplification audit."""
from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sequence_transport as st


def test_user_sequence_exactly_decoded():
    expected = [1, 2, 1, 2, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 3, 2, 2, 3]
    np.testing.assert_array_equal(st.multipliers("user"), expected)
    assert sum(expected) == 35
    np.testing.assert_array_equal(st.multipliers("reverse"), expected[::-1])
    np.testing.assert_array_equal(st.multipliers("uniform"), np.full(20, 1.75))


@pytest.mark.parametrize("bad", ["", "12 7", "12 20", "123"])
def test_rejects_undefined_pairs(bad):
    with pytest.raises(ValueError):
        st.decode_pairs(bad)


@pytest.mark.parametrize("schedule", ["user", "reverse", "uniform"])
def test_matched_duration_work_and_step_multiples(schedule):
    for cycles in (1, 80, 1280):
        times = st.schedule_times(schedule, cycles)
        assert len(times) == 20*cycles+1
        assert times[0] == 0 and times[-1] == 8
        assert np.all(np.diff(times) > 0)
        expected = np.tile(st.multipliers(schedule), cycles)*8/(35*cycles)
        np.testing.assert_allclose(np.diff(times), expected, atol=3e-15, rtol=0)
        np.testing.assert_allclose(times[::20], np.linspace(0, 8, cycles+1), atol=2e-15, rtol=0)


def test_constant_rotation_matches_independent_amplification_product(monkeypatch):
    theta = 2*np.pi*np.arange(64)/64
    k = np.fft.fftfreq(64, d=1/64)
    omega = .2
    def independent_field(grid, profile):
        assert grid == 64
        return theta, lambda t, q: -omega*np.fft.ifft(1j*k*np.fft.fft(q)).real
    monkeypatch.setattr(st, "field", independent_field)
    finals = []
    for schedule in ("user", "reverse"):
        _, times, states = st.integrate_case("constant-audit", schedule, 40, 64)
        z = -17j*omega*np.diff(times)
        gain = np.prod(1+z+z*z/2+z**3/6+z**4/24)
        expected = (gain*np.exp(17j*theta)).real
        np.testing.assert_allclose(states[-1], expected, rtol=0, atol=1e-10)
        finals.append(states[-1])
    np.testing.assert_allclose(finals[0], finals[1], rtol=0, atol=1e-10)


@pytest.mark.parametrize("profile", ["slow-fast-slow", "fast-slow-fast"])
def test_sequence_recovers_fourth_order_time_convergence(profile):
    errors = []
    for cycles in (160, 320):
        theta, times, states = st.integrate_case(profile, "user", cycles, 64)
        result, _ = st.diagnostics(theta, times, states, profile)
        assert result["conservation_passed"] and result["boundedness_guard_passed"]
        errors.append(result["endpoint_profile_error"])
    assert 14 < errors[0]/errors[1] < 18
