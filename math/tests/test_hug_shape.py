"""Independent small distributions check the spatial-moment diagnostic."""
import numpy as np
import pytest

from hug_shape import energy_shape


def test_two_point_distribution_has_known_center_and_widths():
    u = np.zeros((3, 4, 4, 4))
    u[0, 0, 1, 2] = 1.
    u[2, 3, 2, 0] = 2.
    result = energy_shape(u, box=4.)
    # Positions (-1.5,-.5,.5) and (1.5,.5,-1.5), weights 1 and 4.
    assert result['energy'] == 2.5
    np.testing.assert_allclose([result['center_x'], result['center_y'], result['center_z']], [.9, .3, -1.1])
    assert result['radial_width'] == pytest.approx(np.sqrt(1.6))
    assert result['axial_width'] == pytest.approx(.8)
    assert result['rms_radius'] == pytest.approx(np.sqrt(2.24))


def test_speed_rescaling_preserves_width_and_counts_edge_union_once():
    u = np.zeros((3, 20, 20, 20))
    u[0, 0, 0, 0] = 1.
    u[1, 10, 10, 10] = 2.
    before, after = energy_shape(u), energy_shape(-3*u)
    assert before['edge_energy_fraction'] == pytest.approx(.2)
    assert after['energy'] == pytest.approx(9*before['energy'])
    for key in before.keys() - {'energy'}:
        assert after[key] == pytest.approx(before[key], abs=1e-14)


def test_undefined_or_invalid_inputs_are_rejected():
    for u, box in [(np.zeros((3,4,4,4)),6), (np.ones((3,4,4,4)),0),
                   (np.ones((3,4,4,3)),6), (np.full((3,4,4,4),np.nan),6)]:
        with pytest.raises(ValueError):
            energy_shape(u, box)
