# Mirror-fluid tests: reproduced results

[Hug-ns overview](../README.md) · [Charts and report](../../../../reports/files/fluid-tests-completed.html) · [Verification](verification.json) · [Model and test coverage](MODEL.md)

**All four supplied tests were executed and their current reference results reproduced.** Grid: 33³. Domain side: 6. Viscosity: 0.01. The fluid solver uses Fourier derivatives, projection and two-stage Heun stepping with zero external force.

![Completed fluid tests](../../../../reports/files/fluid-tests-completed.png)

| Test | End time | Result |
| --- | --- | --- |
| [Small break](results/fluid-break-to-2.json) | 2.00 | Peak spin falls from 14.789406 to 5.070971. Accumulated maximum spin reaches 17.646129. |
| [Five source configurations](results/fluid-sources.json) | 0.50 | Symmetric starts preserve symmetry. Reflected source cases have opposite signed D and equal unsigned E. |
| [Mirror-averaging intervention](results/fluid-pitchfork.json) | 0.24 | At 0.12 the code replaces the field by its mirror average, then evolves for another 0.12. |
| [Signed q record](results/fluid-correction.json) | 0.16 | Positive and negative inputs give opposite q. The zero-break control remains at numerical roundoff. |

## Measurements

```text
M[u] = velocity reflected across x, including reversal of its x component
E(t) = norm[u(t) - M[u(t)]] / norm[u(t)]
D(t) = <u(t) - M[u(t)], phi>
q_next = q + dt * (-0.2*q + 0.8*D)
W(t) = max_x |curl u(x,t)|
I(t) = integral from 0 to t of W(s) ds
```

The inner product is the spatial mean of component products. The antisymmetric template phi has unit norm under that inner product. E is nonnegative and dimensionless. D and q retain direction; q is a filtered record and never enters the fluid equation. W and I use the simulation's length and time units. Saved times are rounded to three decimal places; I is accumulated at every step.

Signed inputs test **sign preservation**. Zero-break inputs test **preserved symmetry**. Averaging with the mirror is an imposed reset of the velocity field. That test does not demonstrate spontaneous restoration or a pitchfork bifurcation.

The unperturbed pair is independent of z. The added localized perturbations depend on z. These starts differ from both the earlier hugged ring and the matched-stretch field.

## Reproduction and the time-window correction

The first correction script stopped at 0.12, while its reference JSON came from 0.16. Both windows were run. The latest uploaded script uses 0.16 and matches the tested executable code; the only difference from the local diagnostic copy is its module description. [Version check](version-2-check.json). The 0.12 result is retained under an explicit filename.

Reference comparisons pass with absolute tolerance 1e-12 plus relative tolerance 1e-10. The largest absolute difference among current reference values is 1.32e-13. These are reproducibility checks at the supplied grid; they do not establish grid convergence.

From this folder:

```sh
python -m pip install -r requirements.txt
python run_all.py
```

New calculations go into `rerun/`. The saved `results/` and supplied `references/` remain available for comparison.
