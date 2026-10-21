# Working equation: an active, asymmetric guard

**Working draft — completed comparisons cover the [local two-SIMS example](sims-active-guard/README.md) and a [bounded Navier–Stokes realization](sims-navier-stokes-guard/refinement/README.md). The coupled laws remain provisional. The 64² refinement passes the unchanged detailed-field tolerance; both older failures remain preserved.**  
[Connected tests](README.md) · [Sources and earlier equation audit](sources-and-details.md)

## The question

**How do participants change through an interaction when the boundary regulating that interaction also changes, and the influence can differ in each direction?**

The guard belongs to the relationship. It remains a state of the model throughout the evolution. It can regulate what is transmitted, resisted or allowed, and the interaction can change that regulation. It need not be a separate participant or a physical wall.

The earlier document supplies a starting idea. Its terminology and formulas remain provisional. This revision does not impose equal weighting, identical responses, maximization of a particular objective, or permanent operation at a threshold.

## The essential correction

For participant i, with j acting on i:

$$
\boxed{
\frac{d x_i}{dt}
=
f_i(x_i,t)
+
\sum_{j\ne i} R_{i\leftarrow j}
\bigl(x_i,x_j,g_{i\leftarrow j},t\bigr)
}
$$

$$
\boxed{
\frac{d g_{i\leftarrow j}}{dt}
=
H_{i\leftarrow j}
\bigl(x_i,x_j,g_{i\leftarrow j},g_{j\leftarrow i},t\bigr)
}
$$

| Symbol | Meaning |
| --- | --- |
| `x_i` | Participant i's current state; it may have several components. |
| `f_i` | Its own response or tendency, apart from this pair interaction. |
| `g_(i←j)` | The guard state governing the effect of j on i. |
| `R_(i←j)` | What the interaction actually does to i, given the current guard. |
| `H_(i←j)` | How that directional guard responds to the current participants and relationship. |

**Both directions are represented separately:**

$$
g_{A\leftarrow B}\ \text{need not equal}\ g_{B\leftarrow A}.
$$

Their response times, strengths and rules may differ. Compare the response laws after exchanging the participants: opposite-signed changes alone do not establish asymmetry. A larger participant is not automatically assigned more influence; that requires an explicit model choice or evidence.

At every evolution step, the current guard enters the interaction, and the current interaction state enters the guard update. A guard may settle to a steady value; remaining part of the equations does not require it to change at every instant.

This is the **structure** of the revision. The functions R and H still require a stated choice before it becomes an executable model. They must have genuine dependence on the guard and participants respectively; simply listing an unused argument does not create feedback.

<details>
<summary><strong>One explicit, provisional example for two abstract participants</strong></summary>

The following smooth example makes the structure concrete. It is an assistant-proposed modeling choice, not a formula derived uniquely from the original document.

Let x_A and x_B be dimensionless scalar response states, and let g_AB mean g_(A←B), with g_BA meaning g_(B←A). Choose:

$$
\dot x_A=\lambda_A(b_A-x_A)+k_{AB}g_{AB}(x_B-x_A)+s_A(t),
$$

$$
\dot x_B=\lambda_B(b_B-x_B)+k_{BA}g_{BA}(x_A-x_B)+s_B(t),
$$

$$
\tau_{AB}\dot g_{AB}
=
\sigma\!\left(\theta_{AB}+\beta_{AB}(x_B-x_A)\right)-g_{AB},
$$

$$
\tau_{BA}\dot g_{BA}
=
\sigma\!\left(\theta_{BA}+\beta_{BA}(x_A-x_B)\right)-g_{BA},
\qquad
\sigma(z)=\frac{1}{1+e^{-z}}.
$$

- b_i is an individual baseline, not a measured personality or autonomy score.
- λ_i > 0 is the rate of return toward that baseline.
- k_ij ≥ 0 sets how strongly an admitted interaction changes i.
- 0 ≤ g_ij ≤ 1 is the fraction of that interaction transmitted. In this example, a smaller value means stronger restriction.
- τ_ij > 0 is the guard's response time. Different values allow one direction to adapt more slowly.
- θ_ij sets the baseline transmission bias. β_ij sets how a directional state difference changes transmission; its sign is an explicit assumption.
- s_i(t) is a specified disturbance to the abstract state.
- With physical time units, λ_i and k_ij have units of inverse time, τ_ij has time units and s_i has inverse-time units. Here the time scale remains uncalibrated.

For a guard initially in [0,1], its derivative points inward at both endpoints, so this interval is invariant for the continuous equations. No clipping or imposed opening event is part of this model.

The two directional gains k_AB g_AB and k_BA g_BA can differ even when the participants follow the same equation template. Distinct baselines allow persistent differences: in the illustrative fixed-guard case λ_A = λ_B = 1, s_A = s_B = 0, set c_AB = k_AB g_AB and c_BA = k_BA g_BA. At equilibrium,

$$
x_A-x_B=\frac{b_A-b_B}{1+c_{AB}+c_{BA}}.
$$

A shared rule therefore does not force identical states. This identity checks that particular fixed-guard limit; it is not a result from a new simulation.

**Limits of this example:** it models graded transmission and attraction toward another state. It does not yet implement a strict forbidden-action set, repulsion, irreversible memory, an optimizer, or a physical rupture. If the intended guard requires any of those features, R and H must be revised explicitly. There is no programmed “break” to mistake for a discovery.

</details>

<details>
<summary><strong>How the guard remains inside a Navier–Stokes version</strong></summary>

For a constant-density incompressible continuum, retain momentum balance and pressure:

$$
\rho(\partial_t u+u\cdot\nabla u)
=
-\nabla p+\nabla\cdot\left[2\mu D(u)+\Sigma_g\right]+f_{\rm ext},
\qquad
\nabla\cdot u=0,
$$

$$
D(u)=\frac{\nabla u+\nabla u^T}{2}.
$$

The standard fluid part uses stress divergence. This form retains stress contributions when material properties vary; replacing a varying-viscosity stress by only μΔu would omit terms. See the [conservation-law reference](https://www.comsol.com/multiphysics/fluid-flow-conservation-of-momentum-mass-and-energy).

Here Σ_g is a **proposed guard-dependent stress**, with stress units. It must be given a constitutive rule. The abstract relationship does not determine that rule by itself.

A field realization also evolves the participant states and directed guards:

$$
(\partial_t+u\cdot\nabla)x_i
=
\kappa_i\Delta x_i+f_i(x_i,t)
+\sum_{j\ne i}R_{i\leftarrow j}(x_i,x_j,g_{i\leftarrow j},t),
$$

$$
(\partial_t+u\cdot\nabla)g_{i\leftarrow j}
=
\kappa^g_{ij}\Delta g_{i\leftarrow j}
+
H_{i\leftarrow j}(x_i,x_j,g_{i\leftarrow j},g_{j\leftarrow i},D(u),t).
$$

In these equations x_i is an internal state field for population i at each position, not its spatial coordinate or mass fraction. All fields use one carrier velocity u. Diffusivities κ are nonnegative. Modeling separate material velocities or actual chemical concentrations would require additional conservation equations and a new mapping.

### One possible closure, kept explicitly optional

To demonstrate a complete feedback path, use the two-participant example above for the local reactions, add χ_ij times the squared strain magnitude to each sigmoid argument, and choose

$$
\Sigma_g=2\eta_g D(u),
\qquad
\eta_g=a_{AB}g_{AB}+a_{BA}g_{BA},
\qquad
a_{AB},a_{BA}\ge 0,\quad \mu>0.
$$

Thus the sigmoid argument becomes θ_ij + β_ij(x_j-x_i) + χ_ij [D(u):D(u)]. In dimensional form a_ij has viscosity units and χ_ij has time-squared units. At least one nonzero χ and a connected nonzero a are needed for this direct strain→guard→stress loop; the local interaction also needs nonzero k and dependence of its guard on the state.

This supplies an explicit continuous loop:

**flow and participant states → directed guards → interaction and stress → later flow and states.**

For this optional closure, the guard changes viscous resistance. The sign of χ determines whether strain increases or decreases transmission. The assignment of transmission to viscosity is an added modeling choice; it is not an identification supplied by the original idea.

The asymmetry here is between populations A and B. A scalar viscosity does not represent a preferred spatial direction, a moving solid surface, or a membrane with independently specified permeability.

### Conservation and the scope of the example

On a periodic domain, with smooth fields, no external force and the positive-viscosity closure above, integration by parts gives

$$
\frac{d}{dt}\int \frac{\rho}{2}|u|^2\,dV
=
-\int 2(\mu+\eta_g)\,D(u):D(u)\,dV
\le 0.
$$

The spatial integral of momentum is also constant in this setting. A changing guard can change the response and the dissipation rate, but this particular closure supplies no mechanical energy. An active or elastic guard would need its own stress and energy accounting.

Unequal abstract influence is not permission to assign unmatched internal mechanical forces. For a closed physical pair, interaction forces balance; unequal masses can still give unequal accelerations. The directed x_i equations above describe abstract internal responses and do not claim to conserve mass or momentum by summing those states.

The [completed Navier–Stokes SIMS comparison](sims-navier-stokes-guard/refinement/README.md) selects a two-dimensional periodic domain, one smooth starting pattern and fixed parameters for this optional closure. It compares shared and directional evolving guards, equal and unequal weights, frozen guards and controls disabling individual feedback paths. The original sampled averaged responses were stable under grid and time-step refinement. The 64² refinement passes the unchanged detailed-field tolerance; both older failures remain preserved. The optional closure is not a calibrated model of water, smoke, isotopes or people.

</details>

## What to check before calling a future result a framework test

- The guard actually affects the response throughout the simulated interval.
- The participants or flow actually feed back into the guard.
- A control with the guard frozen separates adaptation from a fixed relationship.
- A symmetric control and a directional swap distinguish assumed asymmetry from an indexing error.
- Distinct starting states and individual baselines are retained and reported.
- Any event called a “break” has an observable definition independent of a solver failure or an imposed trigger.

## Status and provenance

This draft follows the user's clarification that the guard stays inside the interaction, relationships may be asymmetric, and earlier labels are placeholders. The initial PDF remains an unchanged source. The prior [equation audit](framework-audit.md) remains available with neutral wording.

**Completed:** structural revision, algebraic limit checks and conservation accounting for the optional fluid closure, followed by the [bounded two-SIMS comparison](sims-active-guard/README.md). Its active-guard runs passed step-refinement and independent-solver checks for the chosen local equations.  
**Also completed:** the [bounded Navier–Stokes comparison](sims-navier-stokes-guard/refinement/README.md), with the guard continuously coupled to transport, strain, influence and viscous stress. Conservation and averaged-response screens passed in the original comparison. The 64² refinement passes the unchanged detailed-field tolerance; both older failures remain preserved.  
**Not yet established:** calibrated response laws, physical time scales, continuum convergence of the spatial guard fields, or validation of the whole framework.

No frozen solver, existing experiment, saved result or active study was changed.

