> **Archived intake review.** This text records the earlier supplied files and the checks available then. Read the [current overview](../../README.md) and [diagnostic status](../../DIAGNOSTIC-STATUS.md) for later repairs and remaining evidence gaps. Numerical claims below retain their original scope. Only relative navigation links were adjusted when archiving.

# Review of the supplied hug-ns files

[Project home](../../../../../README.md) · [Math index](../../../../README.md) · [Experiment overview](../../README.md) · [Recorded checks](../../review-results.json)

**Pressure reporting is corrected. The focused checks record 32 passes and three known numerical failures. Sixteen result tables are included.**

The latest supplied archive is `hug-ns (5).zip`, containing five code/document files and eight result JSON files. Its contents match `hug-ns (4).zip` exactly. Its `solver.py` matches the separately attached file, and its `BOUND.md` matches the separate note. Earlier receipts added the pressure function and bound note; this update adds the result tables. One current copy of each is published. See [the results review](../../RESULTS-REVIEW.md).

[Source hashes](../../provenance.json) and [an exact patch](../../review-changes.patch) record received and published versions. The gate, initial fields, projection, fluid right-hand side, time integrator, cutoff, and time-step scheduler are retained. Code changes affect pressure reporting and runner diagnostics. The documentation has also been revised for clear descriptions of results and limits. Earlier project code and saved runs remain intact.

### What was corrected in the same-point calculation

The uploaded solver includes `pressure_at_peak`. For a divergence-free flow, pressure satisfies `Delta p = -div(advection)`, so `p_hat = div_hat / k²`. The supplied function used the negative of that Fourier coefficient. Its reported pressure contribution therefore had the wrong sign. The main `rhs` projection had the correct sign; this defect was in the diagnostic.

An independent low-mode check found a reported sum of `-1.4224256592`, while the solver's actual local energy rate was `+0.6803376304`. Reversing the pressure sign closes that difference. [Before-correction evidence](../../received-pressure-check.json).

The corrected diagnostic also uses the same projected state and filtered advection as `rhs`, so its three terms describe the implemented update at the same grid point. It reports that point's index and position, signed carrying, pressure, smoothing, and their sum. Low-mode, retained-high-mode, and analytic-shear checks agree with the solver to within `3.34e-16` in this environment.

`run.py` now calls this diagnostic, saves it as `peak_split`, prints the three signed terms, and records `dt` and step count. This fixes the missing console comparison. These additions do not change the simulated velocity arrays.

Each term is a rate of **speed squared divided by two**, not speed itself. Positive means increasing local kinetic energy and negative means decreasing it. The selected fastest grid point may change between samples. These results do not track one particle or one fixed location across time.

### Results supplied during the review

The user quoted Grok reporting that smoothing slows the fastest point in a coarse run, carrying and pressure switch sign, speed rises to `5.22` and later falls to `2.19`, and a finer run is at `6.06` at time `0.53` while still climbing. That message also states the gate is used only at the start.

The supplied JSON files record these runs. Arithmetic checks pass, and the two original starting fields were rebuilt. The later velocity arrays and producing scripts were not supplied. The pressure interpretations need recalculation with the corrected diagnostic. A pressure-diagnostic error does not by itself invalidate velocity results from `rhs`. The newer tables now include the finer run through time 1.6. No new long run was launched for this import.

## The invisible boundary

The boundary is implemented in `solver.py:gate`. It can affect the calculation without being drawn. The scalar gate is one in the interior and zero near the cube faces. In `build`, it shapes the initial velocity through a stream-function construction. It has been retained in this import.

After initialization, `rhs` and `advance` evolve velocity with advection, viscosity, and a pressure projection. They contain no separately evolving boundary, imprint variable, contact rule, or measurement of the separation between two arms. The velocity can change shape; this code has not demonstrated the requested breathing, enveloping motion with a gap maintained between the halves.

There are three distinct objects:

| Object | What this upload implements |
| --- | --- |
| Periodic cube | A computational domain of side 6; opposite faces wrap |
| Initial gate | A scalar function used to prepare the start; no later gate update |
| Two enveloping halves that stay apart | No explicit curves, tracked boundary, or gap check in these files |

The gate being invisible is not an error. The smoothness and support issues below are measurable properties of its implementation.

## How this code differs from the existing project

The uploaded step approximates the unforced periodic incompressible Navier–Stokes equation with Fourier derivatives, a spectral pressure projection, and two-stage explicit Heun stepping. It has positive viscosity (`0.01` for the two ring cases and `0.002` for the tube cases). Pressure is implicit in the evolution projection. The corrected diagnostic reconstructs its gradient at the fastest point; it does not report a pressure difference across a tracked hug boundary.

It uses new starting data. Before numerical projection its velocity is

```text
u = (D_y(w)*p0 + w*D_y(p0), -D_x(w)*p0 - w*D_x(p0), 0).
```

For the `hug` case, `p0 = exp(-((x²+y²-1.5²)²+z²))`. This is a Cartesian stream-function construction. It is different from the existing axisymmetric `psi = r²*z*phi`, its meridional velocity, its separately specified swirl, and its ellipsoidal hug gate in [hugged_ring.py](../../../../hugged_ring.py). Numerical projection can add a third velocity component. The tube functions also define new initial fields.

The earlier project uses finite differences and different starting data. Matching the name `hug`, box length, viscosity, or final time does not make the two trajectories identical. No saved 256³/384³ checkpoints were reused or replaced by this import.

## Findings that need correction before relying on larger runs

### 1. The underlying gate has corners in its transition region

`flat` is smooth at its scalar endpoints, but `gate` composes it with `max(abs(x), abs(y), abs(z))`. That maximum has a corner where two coordinates tie inside the transition band.

At `(2.325, 2.325, 0)`, the measured left x-derivative is `0`; the right derivative approaches `-4.44444444` as the probe spacing decreases from `1e-4` to `1e-6`. Thus the written gate is not differentiable everywhere. Calling the scalar ramp flat does not make this entire spatial gate smooth.

A finite Fourier reconstruction of the arrays is a smooth periodic trigonometric polynomial. This does not repair the smoothness claim about the underlying gate formula, or prove convergence as the grid is refined. A repair would require choosing and documenting a smooth spatial gate while retaining the boundary's intended role; this review does not substitute one.

### 2. The Fourier cutoff retains an alias at N divisible by three

In `grids`, `abs(index) <= n//3` retains the endpoint `N/3`. At N=48, multiplying a retained mode 16 by itself produces mode 32, which aliases onto retained mode -16. The direct cosine-square check leaves a false oscillation of amplitude `0.5` after filtering. The same issue is reproduced at N=12.

The README's 48³ examples encounter this cutoff. A standard repair for this unpadded quadratic product is the strict integer condition `3*abs(index) < N` on each axis, applied consistently to state and products, with a regression check. Changing that filter changes subsequent numerical results, so the supplied solver is preserved and the defect is marked explicitly.

### 3. The runner can exceed its own step estimate

`evolve` uses `round(t_end / dt)` before replacing `dt` with `t_end / steps`. A request of `1.49` estimated steps becomes one step that is 49% larger than the estimate. It also calculates that estimate only once, although the README describes speed increasing in some cases.

A correction should use a ceiling for fixed subdivision, validate finite positive input parameters, and monitor or recompute the step restriction when speed changes. The existing restriction is a numerical estimate, not a guarantee of accuracy or stability. No long-run stability claim is verified here.

The original missing console comparison has been corrected as described above. `compare` still selects the `tube` case and saves both the old global magnitudes and the new same-point split.

## Additional limits in the measurements

- **Zero gate does not mean zero stored velocity there.** Fourier differentiation is global, and the truncated derivative does not obey an exact product rule on arbitrary sampled inputs. Projection/filtering also acts globally. At 16³, the hug's maximum absolute velocity component where the gate is zero is `0.0000215173` before projection and `0.201849` after projection. For the pulled tube, those values are `0.568410` and `0.522915`. These are measurements of this coarse review grid, not error bounds or grid-convergence results. A nonzero collar is compatible with a periodic flow; it means this preparation does not enforce an impermeable zero-velocity wall.
- **“Steepening” is an advection magnitude.** `terms` measures `max |(u·grad)u|` and `max |nu*laplacian(u)|`. Neither directly measures growth of a velocity gradient. The maxima may occur in different places, and the pressure contribution is omitted. Their relative sizes do not decide whether speed or gradients will grow.
- **“Slope” has a specific definition.** `slope_max` is the largest absolute individual partial derivative. Earlier reports use a different gradient norm. Their values cannot be compared directly as the same diagnostic.
- **The checkerboard diagnostic is narrow.** `checker` measures one all-axis alternating mode, which this solver's filter removes; it does not test all short waves. It divides by zero for a zero-energy field. The tested starts have nonzero energy.
- **No restart evidence was supplied.** The latest archive includes result tables, but no velocity checkpoints or exact producing scripts/step records. The reviewed runner now records `dt` and step count, but it does not save evolving divergence or per-step finite checks. The separate review does record those quantities for its small runs.

## Review of BOUND.md

[BOUND.md](../../BOUND.md) includes a correct periodic enstrophy estimate, assuming a smooth divergence-free solution and Euclidean vector norms. The estimate follows from

```text
abs(integral omega · (grad u) omega)
  <= ||omega||_infinity ||grad u||_2 ||omega||_2
   = ||omega||_infinity ||omega||_2^2.
```

The last equality uses periodicity and incompressibility: `||grad u||_2 = ||curl u||_2`. Combining it with the vorticity energy identity gives the inequality in the supplied note. It remains conditional on controlling the time integral of the continuum vorticity maximum. BKM-type continuation criteria for Navier–Stokes are established results; see [Chemin and Zhang, introduction, p. 133](https://smf.emath.fr/sites/default/files/2024-05/ens_ann-sc_49_131-167__sample.pdf).

The conclusion about the finite numerical integral is too strong. The newer `bound-check.json` supplies sample values corresponding to `71.7`, `29.3`, and `86.1`, but no producing vorticity routine, time-quadrature method, or error bounds. The initial maximum was independently reproduced. Integrating the saved samples does not reproduce the recorded integral; see [the results review](../../RESULTS-REVIEW.md). Even verified finite grid measurements would not establish smoothness of the continuous PDE up to time 2. BOUND.md now states the estimate and this remaining requirement directly. Original source hashes and the documentation changes are retained in the provenance record and patch.

## Checks performed on the reviewed version

Environment: Python 3.12.14, NumPy 2.3.5. The repository's pinned requirements provide these numerical dependencies; the uploaded `requirements.txt` itself specifies only `numpy`.

The code import passed the full local repository suite: **243 passed and 3 expected failures**. After adding the sixteen result files, the focused code-and-data suite records **32 passed and 3 expected failures**. These totals come from different test scopes; they are not added together. Each expected failure reproduces one of the three open numerical defects. Syntax and undefined-name checks passed.

| Check | Observation |
| --- | --- |
| Projection on a known periodic field | Divergence `7.22e-16`; repeat projection error `3.33e-16`; energy did not increase |
| Single Fourier shear with exact viscous decay | One Heun step agrees with its known update within `4.45e-16` |
| Low-mode instantaneous energy identity | Rate `-3.89651981755`, matching the independently computed viscous dissipation |
| All four starts | Each advanced 10 steps of `0.001` on 16³, ending at time `0.01`; all fields finite at every step |
| Divergence across those steps | Maximum `1.83e-15` |
| Energy across those steps | Decreased in every case; no step increase detected |
| Corrected same-point split | Agrees with the solver at the selected point on three independent fields |
| Gate join, cutoff alias, step rounding | Three reproduced defects; marked strict expected failures in tests |

These checks verify the listed code behavior on short runs. The longer trajectories described in the original upload, including runs to times 2, 4 and 40, were not replayed in this review. Their accuracy over those intervals needs separate resolution and time-step evidence.

To reproduce just this review from the repository root:

```sh
python -m pip install -r requirements.txt
python -m pytest -q math/tests/test_hug_ns_submission.py -rx
python math/imports/hug-ns/review_checks.py --out scratch/hug-ns-review.json
```

Create `scratch/` first if it does not exist. Fifteen code checks should pass and three should be reported as `XFAIL` (known defects). The separate result-file checks are not included in that command. Unexpected success is also an error so a later repair must update the review.

## Included files

This folder contains the reviewed solver, runner, mathematical notes, sixteen result tables, check scripts and source hashes. The patch records the diagnostic fixes and documentation edits. Main project and math-index links make the work discoverable. The three remaining numerical fixes are recommendations, not implemented changes. The gate is retained, and the solver core before `pressure_at_peak` is byte-identical to the received source.
