# Tests connected to the original framework

**Starting question:** How do participants change through an interaction when the boundary regulating that interaction also changes, and influence can differ in each direction?

The original document is a provisional starting point. The guard remains part of the evolving relationship, and each direction may have its own response. [Working equation and assumptions](interaction-equation.md) remain provisional; the local and Navier–Stokes SIMS examples now have completed comparisons below.

**These completed studies test parts of that question. None implements or validates the whole framework.** The connection is assessed from the actual rules and controls, not from names such as water, smoke or hug.

| Part of the question | Finished study | What it actually tests |
| --- | --- | --- |
| An evolving guard inside fluid motion | [Navier–Stokes SIMS: equal and unequal influence](sims-navier-stokes-guard/refinement/README.md) | Two-dimensional flow with continuous guard feedback. The 64² refinement passes the unchanged detailed-field tolerance; both older failures remain preserved. |
| An evolving guard and unequal influence | [Equal and unequal influence with an active guard](sims-active-guard/README.md) | Two abstract SIMS with the guard updating throughout. Compares equal settings, unequal settings, reversed influence and frozen controls. |
| Shared rules and different local responses | [Coordination with freedom to refuse](../models/decentralized_coordination/README.md) | Hierarchy, local coordination, refusal and position. A common update rule can leave disagreement. |
| Classification and changing group boundaries | [Abstract SIMS populations](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/social-organization/README.md) | Cooperation, group transitions, memory and assigned classification penalties. Material names are symbolic labels. |
| Unequal access under constraints | [Community power](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/community-power/README.md) | How network access, budgets and priorities change fulfillment. Its objective is collective fulfillment, not the PDF's individual autonomy. |
| Roles, expectations and room to refuse | [Adaptable roles](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/adaptable-roles/README.md) | Assigned roles, learned expectations, budgets and handoffs; changing a label and changing access have different consequences. |
| A boundary that retains an imprint | [Moldable hug refinement](hug-molding-refinement/README.md) | Deformation and retained rest-length changes. Individual zigzags remain mesh sensitive; this is a mechanical component model. |
| Reciprocal influence | [Material changes the flow that transports it](water-biology-study/coupled-feedback/README.md) | A controlled two-way feedback mechanism under an assumed viscosity law. It has no autonomy objective or enforcement gate. |

The SIMS links open the separate Nous-Volition repository. Its files were inspected and linked, not reorganized.

[Sources, equation audit and archive](sources-and-details.md)
