# The hug: repaired contact and harder reciprocal Lorenz tests

[Open the report and recorded molding](report.html). [Earlier material model and failures](../clay-hug/README.md).

This version addresses the object's overlap under strong forcing. The hug remains a filled, yielding 2-D material around a fixed rigid object. Its bond rest lengths change under load and relax afterward, producing a fading imprint. “Breathable” describes easy molding; this material has no latched opening rule.

## What changed

The earlier quadratic contact penalty permits penetration at finite force. Increasing its stiffness alone cannot impose a strict boundary. This version replaces it with a localized logarithmic barrier, uses constraints on **every material edge**, adds barriers at the square walls, and rejects invalid integration trials. Positions are never clipped after solving.

For a positive clearance lower bound `g`, activation distance `d = 0.04`, and barrier strength `k`, the contact energy is `b(g) = -k (g-d)^2 log(g/d)` when `0 < g < d`, and zero otherwise. Nonpositive gaps are outside the solver domain. The energy derivative supplies the force. This creates a small artificial stand-off distance, rather than exact zero-gap contact. Default `k = 10000` is a dimensionless modeling choice. The former penalty's strength and the barrier's strength are not physically interchangeable.

For circles and ellipses, each edge's bound comes from minimizing the squared transformed radius analytically along the complete segment. For an ellipse, `g = min(a,b) (sqrt(q_min)-1)` is a conservative distance lower bound. Its positive sign certifies that the segment stays outside the ellipse. Its derivative uses the minimizing point; at endpoints the clamped minimizer is used.

Each edge barrier is weighted by the sum of one third of its adjacent reference triangle areas. Lobed subdivision divides this weight equally across the subsegments. Wall energies use nodal reference-area weights. These are chosen discretization weights, not measured clay contact parameters.

For the three-lobed object, the radial gap has gradient norm at most `sqrt(1 + (0.594/0.702)^2) < 1.32` outside its minimum-radius disk. Each edge is divided into `s` subsegments and each midpoint gap is reduced by `1.32 * edge_length / (2*s)`. An exact segment-to-center minimum-radius guard establishes the domain of that bound. This conservative offset decreases as `s` increases; the subdivision sweep measures its effect. A first-crossing argument also establishes the bound along accepted motion: a path entering the minimum-radius disk must first reach nonpositive object gap, which the positive swept certificate excludes.

## Time integration and geometric certificate

The solver compares a full RK4 trial with two half trials and accepts the two half steps only when its local error estimate passes. All RK stages must stay in the barrier domain with positive rest lengths and noninverted cells. Rejected trials are retried with a smaller step. Recording times are reached directly, without interpolating the saved states.

For each accepted half step, whole-edge clearance at its two endpoints is reduced by a Lipschitz bound times half the maximum endpoint displacement. A positive result certifies every straight interpolation between those states. Circle/ellipse clearance is 1-Lipschitz in endpoint positions; the lobed certificate uses `1.32*(1+1/s)`. Signed triangle area during interpolation is quadratic in time; its minimum is checked analytically. Square-wall constraints are convex, so positive endpoint wall clearances imply positive interpolated clearances.

The certified numerical path is the piecewise straight path between accepted states. This is not a guarantee for an arbitrary RK dense-output polynomial, an exact continuous solution, other shapes, or unlimited forcing. Trial stages can fail and are discarded. Maximum attempted steps are not the actual accepted steps. Acceptance counts, domain rejections, minimum steps, and minimum swept clearance are saved per configuration.

## Both directions of influence

The Lorenz drive applies `-(alpha/10) X grad(A)` to the material, with `A = sum(reference_node_area * radius)`. Material velocity gives `A_dot`. In a two-way run, the reaction changes the Lorenz equation by `X_dot += (alpha/10) A_dot / C`. One-way controls use the identical material forcing rule but omit this reaction. Default reservoir capacity `C = 0.01` is a chosen model parameter and is swept.

Port power into the material is `P = -(alpha/10) X A_dot`. The reservoir receives `-P` in two-way mode, so exchange cancels from the combined energy balance. Material viscous and internal losses remain dissipative. Lorenz native forcing and dissipation remain present: this is not a closed physical fluid system.

The bookkeeping reservoir energy is `E_L = C*(X^2+Y^2+Z^2)/2`; its native power is `C*((10+rho)*X*Y - 10*X^2 - Y^2 - (8/3)*Z^2)`. `X` is treated as an abstract drive coordinate. The coupling is not a derived physical back-reaction law for a convection experiment.

The integrated net port work is a state variable computed with RK4. Separate positive and negative gross work estimates use trapezoidal integration at every accepted half step; their refinement errors are reported. A chaotic trajectory difference alone is not evidence of a physical mechanism: zero-coupling controls, explicit reaction work, matched initial states, short solver comparisons, mesh checks, and longer blocks provide separate checks.

## Actually completed results

<!-- RESULTS -->
Run `build_report.py` after experiments and analysis to populate this section from actual results.
<!-- /RESULTS -->

![Old penetration and new certified clearances](figures/contact.png)

![Reciprocal effects across longer time windows](figures/reciprocal.png)

## Reproduce

Python 3.12 was used. Exact NumPy/SciPy/Matplotlib versions are recorded in `requirements-lock.txt`; compatible ranges are in `requirements.txt`.

Run these commands from this directory (`reports/clay-contact-stress` in the repository).

```text
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-lock.txt
.venv\Scripts\python validate_contact.py
.venv\Scripts\python run_experiments.py
.venv\Scripts\python run_sensitivity.py
.venv\Scripts\python independent_checks.py
.venv\Scripts\python analyze_results.py
.venv\Scripts\python build_report.py
```

The runner resumes already completed named configurations. To rerun independently, use a fresh copy with an empty `data` directory, then run validation first. Four local worker processes are used. `protocol.json` freezes configuration settings and source hashes before computation; each run also records the actual model hash. Full states and geometry are in compressed NPZ files, with times in model units. Independent solver recordings are separate. `checksums.json` records SHA-256 for all deliverable files except itself. `test_playback.cjs` requires Node, Playwright and Chrome; set `PLAYWRIGHT_MODULE` to your Playwright installation and optionally pass a sandboxed preview wrapper.

The version-1 `clay_hug.py` is retained as the material base class, and `baseline_metrics.py` supplies shared measurements. **Run the commands above**, rather than the old support module's standalone runner. `contact_hug.py` overrides its contact law and integration. `data/initial-attempts` preserves the recording-time roundoff failure and its source. `data/low-barrier-attempts` preserves completed trials with barrier strength 100, including an energy-error miss. A too-large finite-difference displacement in the stronger-barrier work check was also retained; the corrected check chooses displacement by the derivative magnitude.

## Limits and uncertainty

This is a dimensionless 2-D annular prototype already surrounding an object. It is not a slab wrapping itself around a new object, calibrated wet clay, a 3-D material, or Navier–Stokes flow. There is no friction, fracture, moving rigid body, pore fluid, or quantum calculation. Positive contact clearance establishes numerical exclusion in this model; it does not validate real clay's constitutive behavior.

The stronger barrier also changes the modeled contact response, so old/new trajectory differences are not purely a solver improvement. Mesh refinement, barrier activation width and strength, and lobed subdivision affect the result. Self-contact forces are not implemented: post-run checks look for proper crossings of nonadjacent edges at every saved frame, while the swept certificates enforce obstacle exclusion and local cell orientation. They are not a full continuous self-collision guarantee. Tiny cell area guards can stop a computation; a stop must remain a failed run.

The independent DOP853 checks test the same equations and therefore validate numerical agreement, not physical accuracy. Three Lorenz initial states and finite block ranges do not establish a statistical confidence interval, an invariant distribution, universal stabilization, or global existence. A finest grid is still an approximation. No claim about the Clay Navier–Stokes problem follows from these tests.

## Primary references

- [Li et al., Incremental Potential Contact (2020)](https://ipc-sim.github.io/file/IPC-paper-fullRes.pdf): barrier contact and geometric feasibility as a general numerical approach. This custom overdamped solver does **not** implement the full IPC algorithm, its Newton solver, friction law or general continuous collision detection.
- [SciPy solve_ivp documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html): the independent DOP853 method and tolerances.

No experimental clay measurements were used in this repair. The comparisons here are against the earlier numerical contact failures, analytic gradient/work identities, independent integration, and refinement—not published measurements of physical clay.
