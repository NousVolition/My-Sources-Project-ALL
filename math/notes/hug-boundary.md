# The hug now prepares the ring's boundary

**Latest:** [The 384³ grid comparison at time 0.16](#finer-grid-at-time-016).

**Completed:** the existing closed hug shapes the existing smooth ring's initial field. Velocity becomes zero smoothly before every cube face, so opposite faces and all their derivatives match. The continuum field remains divergence-free. The Navier–Stokes equation and zero-external-force choice remain unchanged.

## Direct connection to the existing work

- Starting field: `box_experiment.SMOOTH`, the maintained `SmoothRingField` with ring radius 1.5, radial scale √3, axial scale 1, and both rates 1.
- Hug geometry: the actual `hug_envelope.state_at(4.0)` closed pose, including its reciprocal axis stretch.
- Mapping: interpret the hug's horizontal and vertical axes as radial distance and height in a meridional section; rotate the closed ellipse about the z axis.
- Scale: its largest semiaxis is 90% of the cube half-width. For the length-6 cube the radial and axial semiaxes are 2.7 and approximately 2.314815.
- The inner ellipse scaled by 0.7 is an unchanged core. The surrounding layer provides a smooth transition to rest. These spatial scale choices are explicit construction parameters.

The hug's four-second pose is sampled once to prepare fluid time `t=0`. The subsequent fluid motion is governed by Navier–Stokes; the earlier timed animation and imposed pressure-opening rule are not applied during evolution.

## The construction

For the two semiaxes `R_h` and `Z_h`, set

$$
s=\frac{x^2+y^2}{R_h^2}+\frac{z^2}{Z_h^2},\qquad a=0.7^2.
$$

Define `b(t)=exp(-1/t)` for `t>0`, and `b(t)=0` otherwise. Let `W(s)=1` for `s≤a`, `W(s)=0` for `s≥1`, and in between set

$$
t=\frac{s-a}{1-a},\qquad
W(s)=\frac{b(1-t)}{b(t)+b(1-t)}.
$$

All derivatives match at both joins: this is a C∞ spatial gate. The old quintic interpolation was a timing rule for the animation; it is not substituted for this all-orders spatial smoothness.

Apply the hug to the stream function before taking derivatives:

$$
\psi_h=W\psi_0,\qquad
u_r^h=-\frac1r\partial_z\psi_h,\qquad
u_z^h=\frac1r\partial_r\psi_h,\qquad
u_\theta^h=Wu_\theta^0.
$$

The existing ring uses `ψ₀=A r² z φ`. Writing `q=x²+y²` and `W_s=dW/ds`, the implementation evaluates the smooth Cartesian formula

$$
u_x^h=Wu_x^0-\frac{2A xz^2\phi W_s}{Z_h^2},\quad
u_y^h=Wu_y^0-\frac{2A yz^2\phi W_s}{Z_h^2},\quad
u_z^h=Wu_z^0+\frac{2A qz\phi W_s}{R_h^2}.
$$

It has no division by radius. The meridional divergence vanishes by equality of mixed stream-function derivatives, and axisymmetric swirl adds no divergence. This identity extends smoothly through the axis. Simply multiplying the old velocity by W would generally add divergence; a negative-control check detects that error.

The field is identically zero in an open strip next to every cube face. Consequently its periodic extension is C∞, including edges and corners. Viewed on all of space, it is also a smooth, compactly supported, finite-energy initial field. The whole-space and periodic evolutions remain distinct domain choices; the recorded runs use the periodic cube.

## Verification and saved runs

<!-- RESULTS_START -->
**Boundary mismatch:** 0.053709416 before the hug; **0 after**. The zero value follows from the field vanishing in a neighborhood of every face, as well as the sampled face check.

**Force:** `σ = 0`, `P_U = 0`; no external force. The existing solver and pressure projection were reused. No sine/cosine replacement was used.

Four short comparisons reached the same final time 0.08. All values are nondimensional. The raw finite-difference divergence of a sampled continuum field is not exactly zero; the existing projection corrects it before the two evolution branches begin.

| Grid | Time step | Raw sampled max divergence | Initial projection change (L2) | Final projected max divergence | Final projected max gradient | Final projected energy |
| --- | --- | --- | --- | --- | --- | --- |
| 16³ | 0.02 | 1.352074 | 13.2596% | 8.882e-16 | 5.553309 | 39.006958 |
| 32³ | 0.02 | 0.908968 | 3.3073% | 1.776e-15 | 8.078380 | 39.479609 |
| 64³ | 0.02 | 0.514969 | 0.7831% | 3.553e-15 | 9.768902 | 39.555986 |
| 32³ | 0.01 | 0.908968 | 3.3073% | 1.332e-15 | 8.068097 | 39.441322 |

The initial grid correction drops from 13.26% to 3.31% to 0.78% as the grid is refined. It modifies the sampled velocity, including the core; the exact unchanged-core statement applies to the analytic construction before this discrete correction. The gradient values still change with grid resolution, so evolution accuracy is not yet established. The zero seam and analytic smoothness do not rely on that evolution comparison.

![The actual closed hug around the prepared ring and its spatial weight](../figures/hug-ring-boundary.png)

[Saved JSON](../results/hug-boundary.json) · [All rows, including the branch without projection](../results/hug-boundary.csv)

<!-- RESULTS_END -->

The nine new checks cover reuse of the actual hug geometry, unchanged inner core, zero values in face neighborhoods, flat transitions, agreement with independent stream-function differences, local divergence convergence, the failure of naive velocity masking, rejection of open hug poses, and zero force during use of the existing box runner. Seven existing box-integration checks also passed.

## Run this construction

From the repository root:

```sh
python -m pytest math/tests/test_hugged_ring.py math/tests/test_box_integration.py -q
python math/box_experiment.py --profile hug --start projected --box 6 --points 16 32 64 --dt 0.02 --steps 4 --csv scratch/hug-boundary.csv
python math/box_experiment.py --profile hug --start projected --box 6 --points 32 --dt 0.01 --steps 8 --csv scratch/hug-half-step.csv
```

Implementation: [hugged_ring.py](../hugged_ring.py), using [hug_envelope.py](../hug_envelope.py), [initial_field.py](../initial_field.py), and the maintained [box_experiment.py](../box_experiment.py).

**Scope:** the boundary preparation is established analytically and checked numerically. The recorded short evolution runs are implementation checks; their gradients are not yet spatially converged. This supplies no long-time regularity or breakdown result.

<!-- REFINEMENT_START -->
## Finer-grid evolution check

The same hug-shaped ring was evolved to model time **0.08** on **32³, 64³, 128³, 256³** grids with time step **0.001**. The latest time-step check uses **256³ / 0.0005**; the earlier 128³ control is retained. All use the length-6 periodic cube, viscosity 0.01, zero external force (`sigma=0`, `P_U=0`), and the existing centered update with pressure projection.

**Closer agreement:** the final velocity difference is **3.1766% → 0.8613% → 0.2261%** between successive grids. The full-gradient difference is **17.0262% → 6.1472% → 1.7510%**. On the latest pair (128³ vs 256³), the differences are **0.2261%** and **1.7510%**, respectively. Halving the time step on the **256³ grid** changes velocity by **0.0140%** and gradient by **0.0582%**. These comparisons support refinement over this interval, with remaining spatial differences.

The 256³ half-step check is complete. Its starting array is identical to the 256³ full-step run. All previously completed runs and the earlier 128³ time-step control are preserved.

| Grid | Time step | Start max gradient | End max gradient | End energy | End max divergence |
| --- | --- | --- | --- | --- | --- |
| 32³ | 0.0010 | 8.496596 | 8.059470 | 39.407642 | 1.33e-15 |
| 64³ | 0.0010 | 10.710468 | 9.771791 | 39.473866 | 2.89e-15 |
| 128³ | 0.0010 | 11.516387 | 10.348072 | 39.486739 | 5.77e-15 |
| 128³ | 0.0005 | 11.516387 | 10.347793 | 39.484522 | 6.00e-15 |
| 256³ | 0.0010 | 11.715300 | 10.501479 | 39.489805 | 1.27e-14 |
| 256³ | 0.0005 | 11.715300 | 10.500958 | 39.487549 | 1.18e-14 |

![Refinement of the actual hugged ring](../figures/hug-refinement.png)

### Method and verification

- `navier.py` now has an optional NumPy step backend that evaluates the same centered stencils and explicit Euler update as its retained Python loop. The FFT projection uses the same centered derivative symbols. No PDE term or initial profile was replaced.
- **All 146 repository tests passed in the publication check.** These include the array-backend comparisons, independent projection checks, and physical-coordinate interpolation checks. The 256³ extension uses the same numerical source as the completed grid study.
- Each grid samples and projects the same analytic initial field. The projection change depends on the grid. Raw and projected starts are recorded separately.
- Even cell-centered grids are not nested. Fine fields are sampled at the coarse grid's physical coordinates with periodic cubic interpolation; its approximation error is not included in a rigorous error bound. A known smooth wave checks coordinate alignment and fourth-order interpolation refinement.
- Relative L2 differences divide the coarse-minus-fine norm by the interpolated fine norm. Gradient comparisons use all nine centered derivatives. The maximum-gradient plot is a separate pointwise diagnostic.
- Centered derivatives retain checkerboard null modes. Small discrete divergence does not by itself establish continuum accuracy. Positive viscosity and energy decay here are measured properties of these runs, not a general numerical stability guarantee.
- The largest sampled gradient decreases within every run, but the maxima still differ between grids. These results do not establish a blowup or all-time smoothness result.

### Reproduce without duplicating completed work

```sh
python -m pytest math/tests/test_array_step.py math/tests/test_refinement_comparison.py math/tests/test_projection.py math/tests/test_hugged_ring.py math/tests/test_box_integration.py math/tests/test_navier.py math/tests/test_legacy_parameters.py -q
python math/hug_refinement.py
python math/extend_hug_refinement.py --half-step
```

The runner retains full initial/final arrays and per-run measurements in ignored `scratch/hug-refinement/`. It reuses completed runs only when the source fingerprint and parameters match. Reviewed JSON/CSV are in `math/results/`; the earlier four-run boundary study is preserved separately. NumPy and SciPy are already in the repository requirements.

[Refinement JSON](../results/hug-refinement.json) · [Rows by time](../results/hug-refinement.csv)

**Completed next step:** [The continuation to time 0.16](#longer-run-to-model-time-016) is recorded below.

<!-- REFINEMENT_END -->

<!-- LONGER_START -->
## Longer run to model time 0.16

**Completed:** the 128³ / 0.001, 256³ / 0.001, and 256³ / 0.0005 runs now reach model time **0.16**. Each continued directly from its saved **0.08** velocity array. The earlier interval was reused, and the original six-run refinement record is preserved.

**The largest sampled gradient decreased at every recorded time in all three continuations.** On 256³ with the smaller time step, the largest sampled gradient changes from **10.500958** at 0.08 to **9.600445** at 0.16. Energy changes from **39.487549** to **38.702630**. The final centered divergence is **1.31e-14**.

| Comparison at time 0.16 | Velocity difference | Full-gradient difference |
| --- | --- | --- |
| 128³ vs 256³, both dt 0.001 | 0.3287% | 1.8051% |
| 256³, dt 0.001 vs 0.0005 | 0.0240% | 0.0819% |

The differences between calculations are larger at time 0.16 than at 0.08, even though the largest sampled gradient decreases. These are different measurements. The grid differences remain larger than the time-step differences, so spatial resolution accounts for the larger remaining discrepancy in these comparisons.

![Continued flow and numerical comparisons](../figures/hug-longer-time.png)

| Grid | Time step | Max gradient at 0.08 | Max gradient at 0.16 | Energy at 0.16 | Max divergence at 0.16 |
| --- | --- | --- | --- | --- | --- |
| 128³ | 0.0010 | 10.348072 | 9.490098 | 38.701445 | 6.22e-15 |
| 256³ | 0.0010 | 10.501479 | 9.601384 | 38.706818 | 1.31e-14 |
| 256³ | 0.0005 | 10.500958 | 9.600445 | 38.702630 | 1.31e-14 |

### What was preserved and checked

- The same smooth initial field, length-6 periodic cube, viscosity 0.01, centered update, and pressure projection are used. **External force remains zero** (`sigma=0`, `P_U=0`), and the scalar stays zero.
- The numerical source fingerprint matches the earlier runs. Every continuation records hashes of its source checkpoint and of the continuation code. The same grid and time step continue from that checkpoint; the clock is not reset.
- The saved velocity is restored exactly. No new initial projection is applied at the restart. Pressure is recomputed by the existing projection during each new time step.
- A small-grid regression check confirms that the restarted final array is **bit-for-bit equal** to uninterrupted evolution. A second check rejects a checkpoint from different numerical source code. Both checks passed.
- Comparisons reuse the earlier physical-coordinate interpolation and full-gradient measurement. These are differences between numerical runs, not exact-solution error bounds. The centered-grid null modes and finite resolution remain limitations.

The result describes this finite interval. It does not establish an all-time smoothness or breakdown proof.

### Continue or reproduce

```sh
python -m pytest math/tests/test_hug_continuation.py -q
python math/continue_hug_refinement.py
```

The continuation command uses the original checkpoint arrays under ignored `scratch/hug-refinement/`. On a fresh checkout, first follow the earlier refinement commands to generate those arrays. Completed matching continuations are cached and reused. The saved numerical summaries and figures can be read immediately without running the simulations.

[Continuation code](../continue_hug_refinement.py) · [Results and provenance](../results/hug-longer-time.json) · [Measurements by time](../results/hug-longer-time.csv)

**Completed next step:** [The 384³ grid comparison](#finer-grid-at-time-016) is recorded below.

<!-- LONGER_END -->

<!-- GRID384_START -->
## Finer grid at time 0.16

**Completed:** a new **384³** run reaches model time **0.16** with time step **0.001**. The saved **256³ / 0.001** result is reused. Both use the same initial-field construction, length-6 periodic box, viscosity 0.01, and **zero external force**.

**Both measured differences are smaller on the new grid pair.** The final velocity difference changes from **0.3287%** (128³ vs 256³) to **0.0619%** (256³ vs 384³). The full-gradient difference changes from **1.8051%** to **0.3509%**.

| Comparison at time 0.16 | Velocity difference | Full-gradient difference |
| --- | --- | --- |
| 128³ vs 256³, dt 0.001 | 0.3287% | 1.8051% |
| 256³ vs 384³, dt 0.001 | 0.0619% | 0.3509% |
| Earlier 256³ time-step control | 0.0240% | 0.0819% |

The grid changes have different refinement ratios: **2** and **1.5**. A smaller difference on the new pair therefore does not, by itself, establish a convergence order. The time-step control remains at 256³; a 384³ half-step run has not been performed.

![Comparison including the 384 cubed grid](../figures/hug-grid-384.png)

| Grid | Initial max gradient | Final max gradient | Final energy | Final max divergence |
| --- | --- | --- | --- | --- |
| 128³ | 11.516387 | 9.490098 | 38.701445 | 6.22e-15 |
| 256³ | 11.715300 | 9.601384 | 38.706818 | 1.31e-14 |
| 384³ | 11.751008 | 9.620783 | 38.707794 | 1.95e-14 |

The new 384³ run starts with a sampled maximum gradient of **11.751008** and ends at **9.620783**. This maximum is a different measurement from the relative L2 difference of the full gradient tensor.

### Verification and reuse

- The numerical source fingerprint matches the previous runs. The existing centered update, pressure projection, analytic starting field, and force settings were retained.
- The 384³ field was sampled on its own grid from the same analytic construction. It was not created by stretching a saved 256³ result.
- Checkpoints at **0.04, 0.08, 0.12, and 0.16** preserve the new work. Restarts use the saved velocity exactly, without a new initial projection; the scalar remains zero.
- Six targeted checks passed: four field-comparison checks, including a new **3:2 grid-alignment check**, and two checkpoint-restart checks. The new transfer ratio is tested at matching physical coordinates.
- The comparison uses the existing periodic cubic interpolation and centered-gradient tensor. The interpolation error is not enclosed by a rigorous bound. Centered-grid null modes and finite resolution remain limitations.

These are finite-time numerical measurements. They do not establish an all-time smoothness or breakdown proof.

### Reproduce this comparison

```sh
python -m pytest math/tests/test_refinement_comparison.py math/tests/test_hug_continuation.py -q
python math/extend_hug_refinement.py --finer-grid 384
```

This command extends the saved longer-time study. On a fresh checkout, first follow the preceding refinement and continuation commands to create the earlier checkpoint arrays. Completed matching runs and checkpoints are reused; raw arrays remain under ignored `scratch/hug-refinement/`.

[Measured results](../results/hug-grid-384.json) · [Rows by time](../results/hug-grid-384.csv) · [Shared extension runner](../extend_hug_refinement.py)

**Next numerical check:** compare the 384³ result with half the time step at the same final time before extending the interval further.

<!-- GRID384_END -->
