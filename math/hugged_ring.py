"""Use the existing closed hug to give SmoothRingField a smooth compact support.

The closed arm curve supplies the ellipse in an r,z section. Rotating that
ellipse about the z axis supplies an ellipsoidal envelope. A C-infinity spatial
gate multiplies the stream function before differentiation, preserving
incompressibility. This prepares initial data; it imposes no time-dependent
wall, pressure threshold, or external force during Navier--Stokes evolution.
"""
from dataclasses import dataclass, field
import math

from hug_envelope import state_at
from initial_field import SmoothRingField


def flat_gate(squared_radius, inner_fraction=0.7):
    """Return W(s), dW/ds; W=1 in the core and W=0 outside the ellipse.

    For a=inner_fraction**2 and t=(s-a)/(1-a), use
    W=b(1-t)/(b(t)+b(1-t)), where b(t)=exp(-1/t) for t>0, else 0.
    Every derivative matches at the two joins. Stable logistic evaluation
    avoids overflow. The older hug's quintic *time* interpolation is unchanged.
    """
    if not math.isfinite(squared_radius) or not 0 < inner_fraction < 1:
        raise ValueError("Require finite squared radius and 0 < inner_fraction < 1")
    a = inner_fraction**2
    if squared_radius <= a:
        return 1.0, 0.0
    if squared_radius >= 1:
        return 0.0, 0.0
    t = (squared_radius-a)/(1-a)
    if t <= 0:
        return 1.0, 0.0
    if t >= 1:
        return 0.0, 0.0
    log_ratio = -1/t+1/(1-t)
    e = math.exp(-abs(log_ratio))
    if e == 0.0:
        return (0.0, 0.0) if log_ratio > 0 else (1.0, 0.0)
    w = e/(1+e) if log_ratio >= 0 else 1/(1+e)
    derivative = -e/(1+e)**2 * (1/t**2+1/(1-t)**2)/(1-a)
    return w, derivative


@dataclass(frozen=True)
class HuggedRingField:
    base: SmoothRingField = field(default_factory=lambda: SmoothRingField(radial_scale=math.sqrt(3)))
    box_length: float = 6.0
    hug_seconds: float = 4.0
    inner_fraction: float = 0.7
    margin_fraction: float = 0.1
    radial_extent: float = field(init=False)
    axial_extent: float = field(init=False)
    arm_scale: float = field(init=False)

    def __post_init__(self):
        if not math.isfinite(self.box_length) or self.box_length <= 0:
            raise ValueError("Box length must be finite and positive")
        if not 0 < self.inner_fraction < 1 or not 0 < self.margin_fraction < 1:
            raise ValueError("Core and margin fractions must lie strictly between 0 and 1")
        pose = state_at(self.hug_seconds, points=3)
        if not 3.0 <= self.hug_seconds <= 5.0 or not pose["connected"]:
            raise ValueError("Use a closed hug pose (seconds 3 through 5) for boundary preparation")
        sx, sz = pose["eigenvalues"]
        scale = (self.box_length/2)*(1-self.margin_fraction)/max(sx, sz)
        object.__setattr__(self, "arm_scale", scale)
        object.__setattr__(self, "radial_extent", sx*scale)
        object.__setattr__(self, "axial_extent", sz*scale)

    def gate(self, x, y, z):
        s = (x*x+y*y)/self.radial_extent**2+z*z/self.axial_extent**2
        return flat_gate(s, self.inner_fraction)

    def stream_function(self, r, z):
        return self.base.stream_function(r, z)*self.gate(r, 0, z)[0]

    def velocity(self, x, y, z):
        w, dw = self.gate(x, y, z)
        if w == 0 and dw == 0:
            return 0.0, 0.0, 0.0
        original = self.base.velocity(x, y, z)
        if w == 1 and dw == 0:
            return original
        q = x*x+y*y
        factor = self.base.meridional_rate*self.base.envelope(x, y, z)*dw
        # Product-rule terms from psi_h=W*psi are essential. Merely multiplying
        # velocity by W would usually introduce divergence in the outer layer.
        radial_correction = -2*factor*z*z/self.axial_extent**2
        axial_correction = 2*factor*q*z/self.radial_extent**2
        return (w*original[0]+x*radial_correction,
                w*original[1]+y*radial_correction,
                w*original[2]+axial_correction)

    def velocity_cylindrical(self, r, z):
        if r < 0:
            raise ValueError("Cylindrical radius must be nonnegative")
        return self.velocity(r, 0, z)

    def arms(self, points=65):
        pose = state_at(self.hug_seconds, points=points)
        return {name: [[self.arm_scale*x, self.arm_scale*z] for x,z in pose[name]]
                for name in ("Norm", "Entropy")}

    def metadata(self):
        return dict(base="existing SmoothRingField", hug_source="hug_envelope.state_at",
                    hug_seconds=self.hug_seconds, box_length=self.box_length,
                    radial_extent=self.radial_extent, axial_extent=self.axial_extent,
                    inner_fraction=self.inner_fraction, margin_fraction=self.margin_fraction,
                    construction="C-infinity gate on stream function; swirl uses the same gate",
                    evolution="initial-data preparation only; f(x,t)=0 in the box runner")
