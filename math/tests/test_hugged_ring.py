"""Check that the existing hug actually prepares the original ring's boundary."""
import math

import numpy as np
import pytest

from box_experiment import SMOOTH, boundary_jump, hug_for_box, run_case
from hugged_ring import HuggedRingField, flat_gate
from hug_envelope import state_at
from navier import PurePythonNavierStokes3D


def div(velocity, point, h):
    result = 0.
    for axis in range(3):
        plus, minus = list(point), list(point)
        plus[axis] += h
        minus[axis] -= h
        result += (velocity(*plus)[axis]-velocity(*minus)[axis])/(2*h)
    return result


def test_enclosure_uses_the_existing_closed_hug_geometry():
    field = hug_for_box(6.)
    pose = state_at(4.)
    assert field.radial_extent/field.axial_extent == pytest.approx(pose["eigenvalues"][0]/pose["eigenvalues"][1])
    for name, points in field.arms().items():
        for (x,z), (px,pz) in zip(points, pose[name]):
            assert (x,z) == pytest.approx((field.arm_scale*px,field.arm_scale*pz))
            assert (x/field.radial_extent)**2+(z/field.axial_extent)**2 == pytest.approx(1., abs=2e-15)


def test_original_ring_is_unchanged_inside_the_hug_core():
    field = hug_for_box(6.)
    for point in ((0,0,.3), (1.5,0,.7), (1.2,.4,-.4), (0,1.5,-.7)):
        assert field.gate(*point) == (1.,0.)
        assert field.velocity(*point) == SMOOTH.velocity(*point)


def test_faces_and_a_neighborhood_are_exactly_zero():
    field = hug_for_box(6.)
    for axis in range(3):
        tangents = [i for i in range(3) if i != axis]
        for side in (-1,1):
            for a in (-3,-1,0,.7,3):
                for b in (-3,-.5,0,1,3):
                    for inward in (0,.01,.05,.1):
                        p = [0.,0.,0.]
                        p[axis] = side*(3-inward)
                        p[tangents[0]],p[tangents[1]] = a,b
                        assert field.velocity(*p) == (0.,0.,0.)
    assert boundary_jump("hug",6.,17) == 0.
    assert boundary_jump("smooth",6.,17) > 0.


def test_gate_has_flat_joins_and_stays_monotone():
    core = .7**2
    assert flat_gate(core) == (1.,0.)
    assert flat_gate(1.) == (0.,0.)
    for distance in (1e-2,1e-3,1e-4):
        inside, di = flat_gate(core+distance*(1-core))
        outside, do = flat_gate(1-distance*(1-core))
        assert abs(1-inside) < 1e-35 and abs(di) < 1e-35
        assert abs(outside) < 1e-35 and abs(do) < 1e-35
    values = [flat_gate(s)[0] for s in np.linspace(core,1,101)]
    assert all(a >= b for a,b in zip(values,values[1:]))


def test_transition_velocity_matches_stream_function_derivatives():
    field = hug_for_box(6.)
    for r,z in ((2.1,.5),(1.5,1.4),(.8,1.9)):
        h = 1e-5
        numeric_r = -(field.stream_function(r,z+h)-field.stream_function(r,z-h))/(2*h*r)
        numeric_z = (field.stream_function(r+h,z)-field.stream_function(r-h,z))/(2*h*r)
        ur,uth,uz = field.velocity_cylindrical(r,z)
        assert ur == pytest.approx(numeric_r, abs=2e-8)
        assert uz == pytest.approx(numeric_z, abs=2e-8)
        assert uth == pytest.approx(field.gate(r,0,z)[0]*SMOOTH.velocity(r,0,z)[1], abs=1e-14)


def test_divergence_converges_to_zero_in_the_transition_layer():
    field = hug_for_box(6.)
    points = ((2.1,.2,.5),(1.2,1.1,1.3),(.7,.6,1.9),(0,0,.8))
    errors = [max(abs(div(field.velocity,p,h)) for p in points) for h in (1e-3,5e-4,1e-5)]
    assert errors[1] < errors[0]*.3
    assert errors[-1] < 2e-7


def test_naive_velocity_mask_fails_the_same_divergence_check():
    field = hug_for_box(6.)
    def naive(x,y,z):
        return tuple(field.gate(x,y,z)[0]*v for v in SMOOTH.velocity(x,y,z))
    assert abs(div(naive,(2.1,.2,.5),1e-5)) > .01
    assert abs(div(field.velocity,(2.1,.2,.5),1e-5)) < 2e-7


def test_open_hug_is_not_silently_used_as_a_closed_enclosure():
    for seconds in (0,2,2.999999,5.000001,6,8):
        with pytest.raises(ValueError, match="closed hug"):
            HuggedRingField(hug_seconds=seconds)


def test_existing_runner_uses_hug_without_any_external_force(monkeypatch):
    original = PurePythonNavierStokes3D.step
    calls = []
    def checked_step(self,dt,P_U):
        assert self.sigma == 0. and P_U == 0.
        calls.append(True)
        return original(self,dt,P_U)
    monkeypatch.setattr(PurePythonNavierStokes3D,"step",checked_step)
    rows = run_case("hug","projected",6.,8,dt=.002,steps=1)
    assert len(calls) == 2
    assert all(r["face_jump"] == 0 for r in rows)
    assert rows[-1]["max_div"] < 1e-10
    assert math.isfinite(rows[-1]["energy"])
