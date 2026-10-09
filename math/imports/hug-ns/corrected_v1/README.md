# Separately versioned hug-ns repair

October 9, 2026 · Version `hug-ns-corrected-v1`

[Math index](../../../README.md) · [Repository review](../../../../REPOSITORY-REVIEW.md) · [Historical defect evidence](../REVIEW.md) · [Saved verification](verification.json)

**The three known implementation defects now have a separate corrected version and passing regression checks.** The original source, provenance, supplied tables and three strict expected failures remain historical evidence. The active matched-stretch and vortex studies continue using their own frozen sources and starting fields.

## What changed

| Defect | Repair | Verification |
| --- | --- | --- |
| Corner in the max-coordinate gate | Product of three smooth coordinate gates | At the old coordinate tie, the derivative mismatch falls from 4.44444444 to 1.11e-10 with a 1e-6 probe spacing; reflected ties, a three-coordinate tie, plateau and endpoints also pass |
| Retained Fourier cutoff alias | Keep integer modes only when `3*abs(k) < N` on every axis, for state and quadratic products | Independent explicit convolution checks at N=12, 16, 17, 18 and 48 agree within 2e-15; the old N/3 endpoint is excluded |
| Rounding above the estimated step and using only the initial speed | Recompute the limit every step, shorten to output times, and halve trials that violate predictor/result limits or become nonfinite | The 1.49-step reproducer takes steps 0.075 and 0.03675, ending exactly at 0.11175; a manufactured accelerating state exercises rejection and smaller later limits |

The copied solver retains the domain length 6, the four stream-function formulas, their viscosities (0.01 for rings, 0.002 for tubes), Fourier pressure projection, explicit two-stage Heun method and zero force. It does not subtract or impose a new mean. The pressure diagnostic retains the earlier corrected sign and filtered same-point calculation. The zero-field checker diagnostic now returns zero instead of dividing by zero. Invalid grid/collar/time parameters and nonfinite states are rejected in the applicable entry points.

## The smooth-gate choice changes the starting field

Let `F` be the original scalar flat ramp, zero at/below 0 and one at/above 1. The new gate is

```text
w(x,y,z) = product over a in {x,y,z} of F((L/2 - abs(a) - eps)/eps),
0 < eps < L/4.
```

Each factor is constant one near `a=0`, so the absolute-value corner lies inside a constant plateau. At the two transition endpoints, every derivative of the flat ramp vanishes. Each factor is therefore smooth, as is their product; the zero collar also gives a smooth periodic extension. The plateau remains the cube `abs(a) <= L/2-2*eps` on all axes, and the gate remains zero when any `abs(a) >= L/2-eps`. Values inside overlapping transition bands change. This is an explicitly chosen replacement gate, not a reconstruction of the original field.

On the verification grid N=16, the old and new cutoffs are identical. The measured initial velocity changes below therefore isolate the new gate within the retained initial-construction code:

| Case | Relative L2 change from historical initial velocity |
| --- | ---: |
| hug | 0.103640% |
| sharper | 0.125879% |
| tube | 5.397355% |
| tubes | 0.069819% |

These are changes of initial data, not discretization error estimates. The saved verification includes initial energies, means and collar velocities for both versions. For example, the repaired tube's maximum absolute velocity component in the zero-gate collar is still **0.507414** at N=16. Global Fourier differentiation, the sampled product-rule construction and projection do not enforce a zero-velocity wall. Initial-construction convergence remains unverified.

## Cutoff and time-step scope

For a retained component index, `|k| <= K = floor((N-1)/3)`. An alias from two retained inputs into another retained mode would require `k+l-m = ±N`. This is impossible because `|k+l-m| <= 3K < N`. Applying this independently to each axis removes aliases into retained modes for the quadratic evolution products. It does not make arbitrary, initially sampled products exact or supply a spatial error bound.

The time-step estimate is inherited: the minimum of `0.15*dx/max_speed` and `0.08*dx²/(6*nu)`, with no viscous restriction when `nu=0`. Accepted trials satisfy that estimate at the initial state, Euler predictor and final Heun state. Every accepted step records its actual size, all three limits and the number of rejected trials. Samples include divergence and the signed same-point energy-rate split. Runs carry the version and source hashes, and the command-line tools refuse to overwrite an existing output file.

This estimate is **not an accuracy controller or a proof of stability**. In particular, explicit Heun has no nonzero stability interval on the imaginary axis for pure linear advection. Zero viscosity is supported as an input but does not come with a stability guarantee. Trial halving is a guard against the chosen speed/diffusion thresholds; grid and time refinement are still required for trajectory claims. The runner is for separate checks and does not provide a checkpoint restart system.

## Verified scope

The new version has **28 passing regression tests**. The full repository suite passed **318 tests**, with the original **three strict expected failures** confined to the preserved historical implementation.

- Four N=16 starts, each through **t=0.01**, with ten output intervals. Every accepted state/predictor was finite; recorded divergence stayed below **2.48e-15**, and energy decreased at every recorded interval.
- An exact viscously decaying shear tested with 4, 8 and 16 Heun steps to t=0.2. Maximum errors were **1.25665e-5**, **3.10944e-6** and **7.73371e-7**, consistent with second-order convergence on this analytic solution.
- A separate retained-high-mode field satisfies the instantaneous viscous energy identity and the pressure split check. Manufactured scheduler inputs test acceleration and nonfinite-trial rejection; those inputs are tests of scheduling, not fluid trajectories.
- The new regression tests also check output preservation and that the saved verification matches its source hashes. The repository's original defect tests remain strict expected failures for the historical implementation.

Python 3.12.14 and NumPy 2.3.5 were used. [Verification script](verify.py) · [Regression tests](../../../tests/test_hug_ns_corrected_v1.py). There is no new long-running matrix, physical concentration claim, or comparison pretending that the two versions start from the same field.

From the repository root, using the pinned project environment:

```sh
python -m pytest -q math/tests/test_hug_ns_corrected_v1.py math/tests/test_hug_ns_submission.py -rx
python math/imports/hug-ns/corrected_v1/verify.py --out scratch/hug-v1-verification.json
python math/imports/hug-ns/corrected_v1/run.py hug --n 16 --time 0.01 --samples 10 --out scratch/hug-v1-short.json
```

Create `scratch/` first and choose unused output filenames. The next unresolved step for this version is a separately specified grid/time refinement comparison with a common, documented initial-field construction. These short repair checks do not resolve the existing studies' spatial/time-convergence failures, nonsmooth raw-strain joins, or initial maxima outside the central tube.
