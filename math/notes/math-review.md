# Math review: Navier-Stokes working framework

This review concerns the original research draft. Its singularity conclusion needs the corrections or proofs listed below. The current [hug-shaped ring study](hug-boundary.md) uses the corrected smooth field and zero external force; its completed numerical work is recorded separately.

1. **Initial data at the axis.** With \(r=\sqrt{x^2+y^2}\), the factor \(\exp(-((r-r_0)^2+z^2)/\alpha^2)\) is not differentiable at \(r=0\) when \(r_0>0\). The claim that the displayed fields are smooth Schwartz data is therefore unsupported. A smooth ring-centered replacement is \(\phi=\exp(-(r^2-r_0^2)^2/\alpha^4-z^2/\alpha^2)\), because it depends on \(r^2=x^2+y^2\). Full velocity formulas are in [the smooth-field note](smooth-initial-field.md).
2. **Forcing and the Clay formulation.** Clay's Statement C allows a prescribed smooth force \(f(x,t)\). The draft's force depends on the evolving scalar \(S\), so the coupled model is a different system unless a reduction to the prescribed-force formulation is proved. Also, smooth \(S\) does not in general make \(|\nabla S|\) smooth at zeros of \(\nabla S\); a tanh switch outside that norm does not fix the cusp.
3. **Energy estimate.** A finite \(L^2\) energy bound does not prevent a pointwise or derivative blow-up. The asserted uniform bound on the scalar gradient needs an independent estimate; bounded scalar injection alone does not supply the stated constant \(K\) for the coupled problem.
4. **Pressure.** A pressure Poisson equation follows formally by taking divergence. It does not give an independent global smoothness guarantee while smoothness of \(u\) and the force remains unproved. The index and transform formulas should be rederived carefully.
5. **CKN theorem.** Partial regularity constrains the size of a possible singular set for suitable weak solutions. It does not force a Type-I profile or prove an isolated singular point. The draft needs a construction and estimates for any such conclusion.
6. **Discrete experiment.** The original step in math/navier.py omits pressure. The maintained implementation adds a compatible centered projection, with residual stopping and an optional Fourier solve. The publicly saved Copilot branches contain three other pressure variants: one omits dx squared, and two use a nearest-neighbor Poisson stencil incompatible with their centered gradient and divergence. The new numerical checks do not establish continuum convergence or a singularity theorem.

## Productive next step

Reuse the corrected initial field and its existing proof. The [experiment record](../EXPERIMENTS.md) preserves the earlier comparisons. The [hug boundary and refinement study](hug-boundary.md) now supplies the chosen smooth boundary preparation and the current spatial/time-step measurements. Extend that maintained construction and its runners for follow-up work.

## Primary references

- [Clay Mathematics Institute, official problem description](https://www.claymath.org/library/monographs/MPPc.pdf), especially Statements A-D and the discussion of partial regularity.
- [Caffarelli, Kohn, and Nirenberg (1982), original partial-regularity paper](https://doi.org/10.1002/cpa.3160350604).
