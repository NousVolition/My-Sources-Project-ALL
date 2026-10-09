# Test design: local stability, nonlinear outcomes, and reversibility

| Stage | Measurement | Acceptance/control |
| --- | --- | --- |
| 1. Calibrate RK4 | Eight supplied linear systems; two fixed initial vectors each; compare at T=1 with exp(At). | Steps 0.05, 0.025, 0.0125; require finest relative error <10⁻⁵ and observed order 3.5–4.5. |
| 2. Test local approximation | Supplied cubic system near (0,0) and (1,0); shrink perturbations fourfold in stages. | Use its exact solution. Check that linearization error decreases with the expected leading nonlinear term. |
| 3. Challenge center inference | Three explicitly chosen radial systems with the same Jacobian and eigenvalues ±i. | Compare radius with the exact formula; test stable, neutral and unstable cases separately. |
| 4. Predict population outcomes | Supplied rabbit–sheep model on a 47×47 initial-condition grid. | Compare the saddle tangent with full dynamics; halve RK4 step; independently check two trajectories near the basin boundary. |
| 5. Test reversibility | Evolve, apply the stated reversal R, evolve forward again, then apply R again. | Compare undamped and damped mechanics, RK4 and symmetric Verlet; also test both supplied reversible exercises. |

A failure of a numerical acceptance check stops the scripts before publication. No passing result should be interpreted as proof for an untested system. Exact mathematical properties provide the ground truth for these benchmarks. The recorded protocol specifies all starting states, steps and observation times.
