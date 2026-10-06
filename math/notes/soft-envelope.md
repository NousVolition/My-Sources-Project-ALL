# A soft envelope: hug, imprint, and pressure

[Math guide](../README.md) · [Mirror-family check](norm-entropy-mirror.md)

The contributor describes two mirrored halves that round inward like arms,
connect, breathe gently, and let go. They then distinguish **pressure-triggered
opening** from that timed release: an enclosure holds until a pressure limit
is reached, while retaining an imprint and remaining able to change.

## Two small prototypes

| Prototype | Implemented behavior | Checks |
| --- | --- | --- |
| [Timed hug](../hug_envelope.py) | Straighten, bend inward, join at both tips, gently stretch, release, and open on an illustrative eight-second schedule | Four checks: mirror symmetry, timing, ellipse geometry/eigenvectors, and smooth phase joins |
| [Pressure envelope](../pressure_envelope.py) | Hold closed below a chosen pressure-difference limit; retain a slowly relaxing shape memory; make small breathing changes; open when the limit is reached | Five checks: closure, pressure threshold, retained imprint, mirrored opening, and small area-preserving breathing |

The second prototype has no automatic timed release. Lower pressure does not
rejoin an envelope whose opening has already been triggered. A fresh envelope
starts a new trial. These are prescribed geometry/response experiments, not
calibrated clay mechanics, fracture calculations, fluid solvers, or transport
models. An enclosed curve is tested; physical trapping and permeability are
not simulated.

## What the equations mean

The closed hug is now connected to the original fluid work: [the hug boundary construction](hug-boundary.md) uses the actual closed pose to shape the existing smooth ring's stream function. This supplies smooth, divergence-free initial data with zero velocity near the cube faces. The timed motion and pressure-opening law remain separate; they are not added as forces to the fluid equation.

For the pressure example, `p = (inside pressure - outside pressure) / opening limit`.
This prototype accepts nonnegative differences only. Opening starts at `p >= 1`.
The limit is chosen for the experiment, not inferred from a material.

A stored imprint `m` follows `dm/dt = (p - m)/tau`. It responds with a time
constant of 0.8 seconds when loading and relaxes over 6 seconds when unloading.
Those choices make memory visible; they are not measurements. The step uses
the exact exponential update for a constant pressure over that step.

The horizontal scale is `s = 1 + m*(0.04 + 0.02*sin(2*pi*t/3))` and the vertical
scale is `1/s`. Their product is one, preserving the closed ellipse's area.
The horizontal and vertical directions are eigenvectors of this **axis-stretch
part**, with eigenvalues `s` and `1/s`. The complete bending motion is nonlinear.

Both arms reuse one geometry implementation. Here they are spatial mirrors
across `x = 0`. The earlier field experiment instead mirrored velocity values
in a neighboring box; the operations have different meanings and results.

## The graph's height and width

The right panel in the ring-profile figure uses `f(x)/f(0) - 1`, separately for
each curve, and zooms into `x` between -0.3 and 0.3. The center is therefore zero
and heights express relative increases. This is a change in the displayed
quantity. Multiplying `f(x)` by `A > 1` makes values taller; using `f(x/s)` with
`s > 1` stretches horizontally. Changing graph limits only changes the view.

## Response modes to explore next

The contributor proposes **fog, rain, ice, soft clay, hard clay, and broken**
as response names for points, ideas, or interactions. These are currently
proposed labels; no classifier or validated prediction is implemented.
The contributor's later sketch is `fire -> clay <- water`: clay as a candidate
adaptable intermediary, rather than an automatic pairing with the direct
opposite. This remains a brainstorm, not a fixed response rule. An accidental
"yes" response did not resolve the definition of pressure.
The latest proposed distinction is adaptation that preserves the boundary
versus repeated pressure followed by boundary failure. For conversation data,
candidate observations include the stated limit, subsequent repetitions,
whether the input changes in response, and the observed outcome. An outlier
is an unusual observation; that alone does not establish harmful behavior
or an intention to break a boundary. These labels concern observable exchanges,
not a diagnosis or a fixed category of person.
The meaning and observable measure of "pressure" for those data still need
to be specified. Save examples, their prior context, and the chosen response
before fitting thresholds. Material state cannot in general be predicted
from physical pressure alone.

```sh
python -m pytest math/tests/test_hug_envelope.py math/tests/test_pressure_envelope.py -q
```
