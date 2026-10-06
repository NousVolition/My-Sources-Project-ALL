"""Independent checks of the requested flip and its interface effects."""
import math
import pytest

from reflection_probe import (FIELD, divergence, face_jump, joined_velocity,
                              neighboring_copy, normal_slopes)


def test_original_box_is_preserved():
    for point in ((-2.5, .4, -.2), (0., 0., .4), (1., 1., 1.), (3., .3, .2)):
        assert joined_velocity(*point) == FIELD.velocity(*point)


def test_matched_flip_joins_full_vector_values_across_face():
    for y, z in ((0., 0.), (.2, -.3), (1., .5), (-1.1, .8)):
        assert neighboring_copy(3., y, z) == pytest.approx(FIELD.velocity(3., y, z), abs=1e-14)
    assert face_jump('matched_flip') < 1e-14


def test_plain_sign_flip_only_hides_the_gap_on_the_center_line():
    assert neighboring_copy(3., 0., 0., 'sign_flip') == pytest.approx(FIELD.velocity(3., 0., 0.), abs=1e-14)
    assert face_jump('sign_flip') > 1e-3


def test_matched_flip_has_opposing_limiting_slopes():
    # Independent analytic derivative of -x exp(-(x^2-2.25)^2/9) at x=3.
    expected = 26*math.exp(-5.0625)
    errors = []
    for h in (1e-3, 1e-4, 1e-5):
        left, right = normal_slopes('matched_flip', h)
        assert left > 0 and right < 0
        errors.append(max(abs(left-expected), abs(right+expected)))
    assert errors[1] < errors[0]/5
    assert errors[2] < errors[1]/5
    assert errors[-1] < 1e-5


def test_matched_flip_creates_divergence_away_from_the_join():
    # div(C(y,z)-u(x-L,y,z)) = d_y C_y+d_z C_z.
    # On y=z=0 the independent closed form is -52 exp(-5.0625).
    expected = -52*math.exp(-5.0625)
    for x in (3.5, 4.2, 6.5, 8.5):
        got = divergence(neighboring_copy, (x, 0., 0.), h=1e-5)
        assert got == pytest.approx(expected, abs=2e-8)


def test_repeat_and_sign_flip_preserve_interior_divergence():
    # Does not cover the jump distribution on an interface.
    for mode in ('repeat', 'sign_flip'):
        for point in ((4.2, .3, -.4), (6., 0., .4)):
            velocity = lambda x, y, z: neighboring_copy(x, y, z, mode)
            assert abs(divergence(velocity, point)) < 1e-7
