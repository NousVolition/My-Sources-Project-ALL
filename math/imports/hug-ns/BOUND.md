> **Review correction:** The conditional continuum argument below needs an actual solution and a bound on its vorticity. A finite grid estimate such as `86.1` does not certify continuum smoothness, including on the measured interval. The final paragraph overstates what those reported numbers establish. The numbers were supplied without logs. [Review and derivation](REVIEW.md#review-of-boundmd). The received text is preserved below.

<!-- END IMPORT REVIEW BANNER -->

# The bound, and where it stops

The equation after the start is

    d_t u + (u·grad) u = -grad p + nu * laplacian u
    div u = 0

The gate is not in this. It only built the start.

## What the equation does give

Let w = curl u, the spin. Taking the curl of the equation and dropping the pressure, because the curl of a gradient is zero,

    d_t w + (u·grad) w = (w·grad) u + nu * laplacian w

The first term on the right is the stretch. It can make the spin larger. The last term smooths it.

From that, the size of the spin over the whole cube satisfies

    d/dt ||w||_2^2 + 2 nu ||grad w||_2^2 <= 2 ||w||_∞ ||w||_2^2

Smoothing is on the left. Stretching is on the right, and it is paid for by the biggest spin, not by the average spin.

## The bound that would finish step 10

If the biggest spin stays integrable in time,

    integral from 0 to T of ||w(t)||_∞ dt  <  infinity

then the solution stays smooth up to T, and a little past it. That is the Beale–Kato–Majda bound. It is a real theorem. It is conditional.

Step 10 needs this integral to be finite for every finite T, for every smooth divergence-free start, with no extra stop put in by hand.

## Why the runs do not meet it

The inequality does not close. In three dimensions the biggest spin is not controlled by the average spin. A start can raise the biggest spin without raising the energy. The pulled tube did that, then faded. One faded climb does not bound every climb.

On the pulled tube, 48 cubed, nu 0.002, the biggest spin starts at 71.7 and is 29.3 at time 2. The time integral through time 2 is 86.1. Finite on that interval means this run stayed smooth on that interval. It was already known to stay smooth there. It is not the bound for all time.
