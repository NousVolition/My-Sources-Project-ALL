"""Test an upside-down neighboring copy of the maintained smooth field.

This probes initial data and an interface; it is not a time evolution.
The endpoint-matched version is the scalar sketch's vertical reflection,
applied component by component to the velocity vector. It is not a spatial
coordinate reflection or a uniquely implied interpretation of a 3-D mirror.
"""
import argparse
import json
import math
from pathlib import Path

from initial_field import SmoothRingField

FIELD = SmoothRingField(radial_scale=math.sqrt(3.0))
BOX = 6.0


def neighboring_copy(x, y, z, mode='matched_flip', field=FIELD, box=BOX):
    """Evaluate a copy to the right of x=box/2; the original box is untouched.

    repeat:       u(x-L,y,z)
    sign_flip:   -u(x-L,y,z)
    matched_flip: u(-L/2,y,z) + u(L/2,y,z) - u(x-L,y,z)

    matched_flip reflects around the midpoint of the endpoint values for
    each component and each (y,z). This matches values at the joining face.
    The copy occupies [L/2,3L/2]; y/z boundary treatment is unchanged.
    """
    if mode not in ('repeat', 'sign_flip', 'matched_flip'):
        raise ValueError('Unknown copy mode')
    value = field.velocity(x - box, y, z)
    if mode == 'repeat':
        return value
    if mode == 'sign_flip':
        return tuple(-v for v in value)
    left = field.velocity(-box / 2, y, z)
    right = field.velocity(box / 2, y, z)
    return tuple(a + b - v for a, b, v in zip(left, right, value))


def joined_velocity(x, y, z, mode='matched_flip', field=FIELD, box=BOX):
    if x <= box / 2:
        return field.velocity(x, y, z)
    return neighboring_copy(x, y, z, mode, field, box)


def norm(vector):
    return math.sqrt(sum(v * v for v in vector))


def divergence(velocity, point, h=1e-5):
    total = 0.0
    for axis in range(3):
        plus, minus = list(point), list(point)
        plus[axis] += h
        minus[axis] -= h
        total += (velocity(*plus)[axis] - velocity(*minus)[axis]) / (2*h)
    return total


def face_jump(mode, field=FIELD, box=BOX, points=33):
    """Maximum vector trace mismatch over a stated sample of the x face."""
    coords = [-box/2 + box*i/(points-1) for i in range(points)]
    return max(norm(tuple(a-b for a, b in zip(
        field.velocity(box/2, y, z),
        neighboring_copy(box/2, y, z, mode, field, box))))
        for y in coords for z in coords)


def normal_slopes(mode, h, y=0.0, z=0.0, field=FIELD, box=BOX):
    """One-sided d(u_x)/dx. Keep distinct interface traces for discontinuities."""
    edge = box/2
    left = (field.velocity(edge, y, z)[0]
            - field.velocity(edge-h, y, z)[0])/h
    right = (neighboring_copy(edge+h, y, z, mode, field, box)[0]
             - neighboring_copy(edge, y, z, mode, field, box)[0])/h
    return left, right


def report():
    result = dict(
        question='Does an upside-down right-hand copy fix the wrap join?',
        interpretation='Vertical reflection of velocity values in a neighboring copy, not a coordinate reflection and not a flip inside the original box.',
        domain='Original x in [-3,3]; neighboring copy x in [3,9]. Only the x-face join is tested; y and z are unchanged.',
        parameters=dict(R=1.5, radial_scale=math.sqrt(3), axial_scale=1.0,
                        meridional_rate=1.0, swirl_rate=1.0),
        units='Nondimensional; velocity, derivative, and divergence are different quantities.',
        face_sampling='33 by 33 points including face boundaries',
        probe_point=[4.2, 0.0, 0.0],
        cases=[],
    )
    phi_edge = FIELD.envelope(3, 0, 0)
    # Here q=9, a^4=9, R^2=2.25: d_x u_x(3,0,0)=26*phi_edge.
    # Flipping the copy reverses that slope; divergence of the matched
    # offset is -2*d_x u_x(3,0,0)=-52*phi_edge throughout y=z=0.
    result['analytic_checkpoint'] = dict(
        original_x_slope_at_join=26*phi_edge,
        flipped_x_slope_at_join=-26*phi_edge,
        matched_flip_divergence_on_center_line=-52*phi_edge)
    for mode in ('repeat', 'sign_flip', 'matched_flip'):
        velocity = lambda x, y, z: neighboring_copy(x, y, z, mode)
        result['cases'].append(dict(
            mode=mode,
            sampled_face_value_jump=face_jump(mode),
            left_velocity_x_at_center=FIELD.velocity(3, 0, 0)[0],
            right_velocity_x_at_center=velocity(3, 0, 0)[0],
            refinement=[dict(h=h, left_slope=normal_slopes(mode, h)[0],
                             right_slope=normal_slopes(mode, h)[1],
                             divergence_inside_right_copy=divergence(velocity, (4.2, 0, 0), h))
                        for h in (1e-2, 1e-3, 1e-4, 1e-5)],
        ))
    result['reading'] = {
        'repeat': 'Original mismatch remains; each copy is divergence-free in its interior.',
        'sign_flip': 'The center-line endpoints meet and interior divergence stays zero, but velocity traces still mismatch elsewhere on the face. This does not establish global incompressibility across the join.',
        'matched_flip': 'All face values match by construction, but one-sided slopes differ and the endpoint offset creates nonzero interior divergence.',
        'scope': 'This tests these two direct velocity-value flips. It does not rule out other reflections or a separately constructed smooth periodic field. No Navier-Stokes time evolution or regularity conclusion follows.'}
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    result = report()
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    for case in result['cases']:
        last = case['refinement'][-1]
        print(f'{case["mode"]}: face jump={case["sampled_face_value_jump"]:.9g}; '
              f'left/right slope={last["left_slope"]:.9g}/{last["right_slope"]:.9g}; '
              f'interior divergence={last["divergence_inside_right_copy"]:.9g}')
