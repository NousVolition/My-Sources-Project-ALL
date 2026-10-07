# Initial fields and the periodic cube

[← Math guide](README.md) · [Saved results](results/box-experiment-results.csv)

This folder makes the submitted box experiments reproducible. It keeps the original radial Gaussian for comparison and provides the corrected smooth field separately. The original advection and diffusion step is retained.

## Review of all submitted code

| Code group | Finding |
| --- | --- |
| Original stream function and velocity | The derivatives are correct for positive radius. The Cartesian axial velocity has an axis cusp when the ring radius is positive and the height is nonzero. |
| First seven tests | All pass. They establish the selected numerical checks, not smoothness or evolution. |
| Seven geometry tests | All pass for the original bump. Symmetry in radius must become symmetry in squared radius for the corrected profile. Finite radial velocity does not prove smoothness. |
| Six comparisons with the cube | All pass for the original step. Tests of a function's attributes or a value differing from speed do not mathematically establish which PDE it solves. The nonzero-divergence experiment does demonstrate that this update lacks a projection. |
| Projection and step wrapper | Fixed the missing dx squared, damped the oscillating iteration, added residual stopping and an optional Fourier solve. The potential must be multiplied by density/dt to interpret it as the time-step pressure approximation. |
| Eight-point, dx=0.75 experiment | This is the same length-6, eight-point case included in the later script. Reproduced with the corrected projection. |
| Length-18 grids | Reproduced at 16 and 32 points. Their spacings differ from the corresponding length-6 grids. Added a 48-point control to hold spacing fixed. |
| Length-6 grids | Reproduced at 8, 16 and 32 points. On the finer grids the largest raw divergence occurs at the periodic seam. |

## Run

From the repository root:

    python -m pip install -r requirements.txt
    python -m pytest
    python math/box_experiment.py

The default experiment is the submitted length-6 box with 8, 16, and 32 points in each direction, four steps of size 0.02, and the original field sampled before either branch evolves. The projection is corrected. Labels use “projection” because the stored correction potential is not pressure itself.

To run the comparisons recorded in the accompanying CSV:

    python math/box_experiment.py --suite --csv scratch/box-results.csv

To select a common smooth start:

    python math/box_experiment.py --profile smooth --start projected --box 18 --points 48

The full suite includes the requested length-18 grids with 16 and 32 points, plus a 48-point grid that matches the length-6, 16-point spacing. It also runs a half-time-step check for the smooth length-6, 32-point case to the same final time, 0.08.

## What was corrected

1. The submitted potential iteration omitted the squared grid spacing. For centered divergence and centered gradient, the Poisson stencil reaches two cells away, and the right-hand-side contribution is **minus 4 times dx squared times divergence**.
2. Undamped Jacobi can oscillate. The standard-library projection uses damping 2/3 and stops according to the equation residual. It raises if the solve fails. The optional NumPy Fourier backend solves the same discrete equation, using sine symbols of centered differences.
3. Sampling now uses the spacing stored on the simulation. The smooth Cartesian field handles the axis exactly; it needs no small-radius replacement.
4. Both experiment branches use zero scalar, zero scalar injection and sigma=0. The last setting removes the tiny constant tail in the original tanh force and isolates fluid dynamics.

The class in math/navier.py keeps its original step method. The new step_with_pressure method calls that update and then project. project defaults to standard-library damped Jacobi; the experiment explicitly selects the faster Fourier backend.

## Results for the original bump, length 6

The initial maximum divergence is measured before any time step or projection.

| Points per direction | Spacing | Initial max divergence | Initial interior max divergence | Max divergence at time 0.08, no projection | Max divergence at time 0.08, projected |
| --- | --- | --- | --- | --- | --- |
| 8 | 0.75 | 0.46351 | 0.46351 | 0.48463 | 2.22e-16 |
| 16 | 0.375 | 0.91375 | 0.23399 | 0.62311 | 6.66e-16 |
| 32 | 0.1875 | 1.71323 | 0.07559 | 1.23564 | 1.78e-15 |

“Interior” excludes the outermost cell layer where the difference stencil wraps across the periodic seam. For the 16- and 32-point grids, the largest initial divergence is at that seam. A finer grid amplifies a fixed boundary jump even while interior difference errors shrink.

At the same spacing 0.375, changing length 6 with 16 points to length 18 with 48 points lowers initial maximum divergence from 0.91375 to 0.23399. The sampled opposing-face velocity mismatch falls from 0.81910 to 8.58e-24. The remaining divergence comes from sampling and differencing a continuum field.

At fixed length 18, the 16- and 32-point raw maxima are 0.22016 and 0.40763: these coarse samples do not show monotone convergence. A bigger box at the same point count also gives a coarser grid.

## Results for the smooth field

The smooth runs use R=1.5, a=sqrt(3), b=1, A=B=1. This matches the original bump's local radial width at its peak; the complete profiles differ. The corrected envelope is symmetric in squared radius, not in radius. Its Cartesian velocity is smooth through the axis.

The smooth comparison projects the sampled start once and copies that common field to both branches before evolution. This avoids mixing the initial sampling correction with the first time-step correction.

At length 6 and 32 points, the initial correction is 0.9636% of the sampled velocity in grid L2 norm. At time 0.08:

| Time step | Max divergence, no projection | Max divergence, projected | Max gradient, projected | Grid energy, projected |
| --- | --- | --- | --- | --- |
| 0.02 | 1.05212 | 1.33e-15 | 4.08213 | 38.80037 |
| 0.01 | 1.03447 | 1.33e-15 | 4.07891 | 38.77337 |

The gradient diagnostic is the largest Frobenius norm of the centered velocity gradient. Grid energy is half the sum of squared speeds times dx cubed, per unit density. All values are nondimensional.

The length-6 smooth field still has a nonzero boundary mismatch. Projection does not repair the underlying choice of periodic extension, and it need not decrease the maximum gradient. Centered differences also miss checkerboard modes. These short experiments check the implementation; they establish neither continuum convergence nor long-time stability.

## Verification

The consolidated repository has 80 passing tests: the 54 earlier local checks, seven initial-field checks, and 19 distinct recovered Copilot checks. Duplicate function implementations and repeated uniform-field checks were combined. The original 20 still pass despite the old velocity's axis cusp. New tests show that cusp and the matching derivatives of the smooth replacement.

Projection tests cover non-unit spacing, odd and even grids, the oscillating Jacobi mode, failure reporting, agreement of both solvers, preservation of the mean, idempotence, orthogonality, and non-increase of energy during projection alone. One test explicitly demonstrates the centered operator's checkerboard blind spot.

The verified environment used Python 3.12.14, NumPy 2.3.5, and pytest 9.1.1. The seven original standard-library initial-field tests now live in math/tests alongside the other checks.

The next accuracy study should choose genuinely periodic initial data, or specify a controlled approximation to the whole-space problem, then refine space and time with the boundary treatment fixed.

The editable derivation is [smooth-initial-field.tex](notes/smooth-initial-field.tex). For the historical projection approach, see [Chorin's 1968 paper](https://math.berkeley.edu/~chorin/chorin68.pdf); the stencil used here follows from composing this code's two centered first derivatives.
