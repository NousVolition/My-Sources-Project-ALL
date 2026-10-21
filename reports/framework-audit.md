# Returning to the original framework

**Audit date:** October 20, 2026, America/Los_Angeles.  
[Connected tests](README.md) · [Sources and archive](sources-and-details.md)

## What happened

The work gradually changed from a question about imposed weighting, autonomy and relational boundaries into questions about pressure thresholds, material deformation, fluid growth and numerical accuracy. Several useful component experiments remain. Their results do not establish the original framework.

The clearest break is documented in the hug-to-fluid construction: the hug's closed shape was sampled **once** to prepare the initial fluid field. Its timed motion and opening rule were absent from the subsequent fluid evolution. Later refinement made that separate fluid calculation more accurate without restoring the missing boundary relationship. [Construction record](https://github.com/NousVolition/My-Sources-Project-ALL/blob/19b6da5d8496173e9e2068842f32eaa9d9126fb9/math/notes/hug-boundary.md#direct-connection-to-the-existing-work).

## 1. The starting point to preserve

**Current clarification:** the earlier terminology and formulas are exploratory. The active guard and potentially asymmetric relationship are retained in the [working equation](interaction-equation.md). Equal weighting and operation at a prescribed boundary distance are not imposed as universal requirements.

Source: the supplied *True_B_Lister_Mathematical_Framework.pdf*, pages 1–3. The entries below describe provisional ideas in the document; they are not final names, required constitutive laws or established empirical results.

| Original element | Meaning to preserve | What an implementation must identify |
| --- | --- | --- |
| Provisional equal-weighting notation: `W(i, α) = W(j, α)` | Imposed equality of weighting across participants | The act α, whose weighting W is imposed, and how enforcement works. Equal weighting does not by itself specify identical movement or identical outcomes. |
| Provisional constrained-choice notation: maximize an objective while remaining near a constraint | Preserve autonomous intent near the enforcement boundary | An actual objective U, available actions, the boundary and a defined distance or margin. |
| The boundary advances or recedes with context or assigned position | The feasible room for action can differ by relationship and position | Which interaction changes the boundary, whose action is constrained, and the direction of influence. |
| Your subsequent “guard” description | The role can be reciprocal and asymmetric: each side may constrain the other differently | Separate directional relationships; reciprocity alone does not imply equal influence. |

The document gives an optimization schema. It does not supply a Navier–Stokes initial field, a material law, a time evolution rule or measured coefficients.

Two mathematical definitions need clarification when implementing it: ordinary distance cannot take the negative values used for an “outside” state in the PDF, so crossing would need a signed margin or an explicit side label. Also, a boundary having no intersection with an action set does not, by itself, show that no feasible action exists. These are definition gaps to resolve while preserving the intended question.

## 2. Which studies actually connect

“Component connection” below means a rule or controlled comparison addresses an identifiable part of the framework. It does not mean the whole framework was implemented, or that a toy result establishes a fact about people.

### Population, access and autonomy components

| Study and evidence | Actual connection | Missing connection or limit |
| --- | --- | --- |
| [Coordination with freedom to refuse](../models/decentralized_coordination/README.md), [methods](../models/decentralized_coordination/METHODS.md) | Compares imposed hierarchy, shared local rules and probabilistic refusal/initiative. Crossed identity and seat comparisons address unequal influence by position. The deterministic common-rule condition has 152 unanimous and 360 frozen-disagreement outcomes among 512 starts. | A common update rule is not a direct implementation of equal W. No PDF autonomy objective, changing boundary or distance constraint. |
| [Abstract SIMS population](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/social-organization/README.md), [model](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/social-organization/model.py) | Cooperation, group membership, rigidity, memory and classification penalties. Renaming controls distinguish label names from assigned consequences. | Ice/water/fog/mist/gas/smoke are symbolic categories. Penalties and adaptation rules are assumptions. No complete weighting/objective/boundary mapping. |
| [Community power](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/community-power/README.md), [methods](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/community-power/methods.md) | Explicit optimization under access, effort and priority constraints. Shows how unequal access or weighting changes who benefits. | Maximizes weighted collective fulfillment, not an individual's autonomous intent near a boundary. The allocation model is static. |
| [Adaptable roles](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/adaptable-roles/README.md), [methods](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/adaptable-roles/METHODS.md) | Assigned roles, expectations, learning, budgets, refusal and handoffs. Separates a label change from a change in practical constraints. | Explicitly chosen update rules; no measured population or calibrated autonomy metric. |
| [One resisting position](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/one-resisting-position/README.md) | Compares a resisting hub, leaf and ring position; distinguishes first agreement from lasting agreement. | Opposition is programmed. Loss of unanimity in that condition is not discovery of an unprogrammed “break.” |

### Boundary memory and reciprocal influence components

| Study and evidence | Actual connection | Missing connection or limit |
| --- | --- | --- |
| [Pressure envelope](../math/notes/soft-envelope.md), [source](../math/pressure_envelope.py) | A chosen opening threshold and an imprint that loads and relaxes at different rates. | Prescribed pressure; opening is latched at normalized pressure 1. It is a prototype, not a discovered boundary law or a flow calculation. |
| [Hug dynamics](hug-dynamics/README.md) and [two-way stress tests](hug-stress/README.md) | Adds deformation/memory dynamics and reciprocal coupling between drive and boundary response. | The threshold is still imposed. Opening does not automatically vent the driving pressure. Missed-pulse and stiffness failures remain part of the record. |
| [Moldable hug refinement](hug-molding-refinement/README.md) and [material revision](clay-gap-refinement/README.md) | Retained deformation, changing rest lengths, contact and reciprocal loading. | A mechanical surrogate for a yielding boundary. No U or relationship-specific constraints. The material revision changes the constitutive model; it is not simply a more accurate rerun. Mesh-sensitive imprints and zigzags remain visible. |
| [Coupled material–fluid feedback](water-biology-study/coupled-feedback/README.md), [source](water-biology-study/coupled-feedback/coupled_model.py) | Material is transported by the flow and changes its viscous stress. Feedback-off controls isolate this specified reciprocal mechanism. | The viscosity law is assumed. No guard, autonomy objective or changing permission boundary. No biological or isotope calibration. |
| [Historical HUG SIMS](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/dynamics-fractals-sims/HUG_SIMS_METHODS.md) and [CLAY SIMS](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/dynamics-fractals-sims/CLAY_SIMS_METHODS.md) | Local states, retained memory and reciprocal network coupling. | The graph coupling is symmetric and loading is shared. In HUG SIMS, the recorded opening event does not feed back into the dynamics. These are earlier component models, retained behind the current reports. |

### Supporting results and illustrations

- [Fluid memory and erasure](fluid-memory-erasure/README.md) and its [extension](fluid-memory-extension/README.md): history-dependent response through the present physical state. Useful supporting comparisons, but no autonomy, classification constraint or acting SIMS. The first study's disturbance-response prediction did not improve. The extension's report records incomplete raw-data publication.
- [Organization pilot](fluid-organization-pilot/README.md): comparisons of prediction and disturbance-specific response. Ordinary deformation was more predictable in some conditions; that does not establish the claimed disturbance-specific benefit.
- [Molecular positional influence](water-molecule-influence/README.md): shared physical laws can produce unequal temporary local responses. This supports an analogy; it does not implement enforcement or autonomous intent.
- [SIMS recovery and entrainment](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/sims-recovery-entrainment/README.md): model-dependent recovery and position effects, with negative long-time results preserved.
- [The records/classifier demonstration](https://github.com/NousVolition/Nous-Volition/blob/0ae1701cfe7570fb54c29d689a69650172c55927/studies/sims-recovery-entrainment/animations/records.md) is conceptually close: a record changes access, access changes state, and classification changes later access. It is a programmed illustration, not an empirical test or an independently discovered mechanism.

## 3. Where the Navier–Stokes equation changed the question

The preserved [initial PDE draft](../archive/navier-stokes-framework.txt) introduced the following force into fluid momentum:

```text
f(S) = -(σ/ρ)(|∇S| - Γcrit)
       × [1/2 + 1/2 tanh((|∇S| - Γcrit)/ε)] e_k

S_t + u·∇S = D ΔS + ρ_ε(x - x₀) P(U)
```

This force is also visible in [the implementation](https://github.com/NousVolition/My-Sources-Project-ALL/blob/19b6da5d8496173e9e2068842f32eaa9d9126fb9/math/navier.py#L39-L80).

| Change | Why it loses or narrows the original connection |
| --- | --- |
| Boundary proximity becomes scalar-gradient magnitude `|∇S|` | No derivation identifies the PDF's distance to an enforcement boundary with a spatial concentration gradient. They are different quantities. |
| The guard becomes a fixed `Γcrit`, resistance scale and direction `e_k` | There is no rule for a context-dependent boundary or separate influence in each direction of a relationship. |
| `P(U)` becomes a prescribed bounded source | Sharing the letter U does not implement maximizing autonomous intent. The code receives a source value; it does not solve the PDF's constrained choice problem. |
| One scalar field replaces participant/act-specific weighting | The draft supplies no representation of `W(i, α)`, its imposed equality, or distinct feasible actions for each participant. |
| ε is reused as source/switch smoothing | The PDF uses ε for proximity to a boundary. A documented change of variables and units would be needed to identify those uses. |
| The question becomes finite-time fluid singularity | A large fluid gradient, threshold crossing, a numerical stop and a change in permitted action are different events. No equivalence was established. |

This is the conceptual drift. It is separate from whether a fluid solver accurately integrates its own stated equations.

### A second break: the hug stops acting during fluid evolution

The [hug-boundary record](https://github.com/NousVolition/My-Sources-Project-ALL/blob/19b6da5d8496173e9e2068842f32eaa9d9126fb9/math/notes/hug-boundary.md#L5-L15) states that the closed pose at animation time 4 is used only to construct fluid time zero. Later runs set external force and scalar source to zero. The hug then supplies an initial shape, not an evolving boundary.

The [vortex matrix protocol](study/protocol.json) likewise specifies an unforced periodic fluid and chosen initial vortices. Numerical acceptance checks are accuracy controls, not guards in the original sense. Its 0–0.40 interval is model time, not a calibrated fraction of a second in a room.

The [hidden-flow note](../math/notes/hidden-flow-evolution.md) already says that its separate trigonometric example was chosen by the assistant and that a mapping from the contributor A/B/C idea was not established. It belongs in the unmapped archive.

### Mathematical claims that also needed correction

The existing [mathematical review](../math/notes/math-review.md) records additional problems in the initial draft:

- The ring profile is not smooth at the axis as originally claimed.
- `|∇S|` need not be smooth at its zeros; applying a tanh outside it does not remove that issue.
- Bounded scalar injection does not supply the asserted gradient bound.
- A finite total-energy bound and a formal pressure equation do not prove global smoothness.
- Partial regularity does not construct or force an isolated Type-I singularity.
- The original time update omitted pressure; a compatible pressure projection was added later.

Those repairs and numerical checks are useful. They neither prove the singularity claims nor supply the missing mapping to the original framework.

## 4. What has been set aside

[The archive](archive/not-yet-connected/README.md) collects the unforced vortex-growth, matched-stretch, hidden-flow, mirror/reduced-model, isotope, quantum and unrelated benchmark branches. “Not yet connected” means no demonstrated mapping to this framework was found in the reviewed records. It does not mean every result is false or that a failed check should disappear.

The old long navigation page has been relocated into that archive. Original scientific files, failed runs, figures, release downloads and checkpoint paths retain their locations. This preserves reproducibility and ongoing work while removing the disconnected entries from the main results page.

## 5. How to keep the starting point from getting lost again

A test should appear as a framework test only when its description identifies:

1. The participant state, its own response, the directional guard state and the ongoing feedback it implements. If weighting, choice or an autonomy objective is used, define it explicitly.
2. Which features are omitted or deliberately held fixed.
3. The rule for each direction of influence, including any assumed asymmetry.
4. The observable called a “break,” and whether the model prescribes it or the experiment tests for it.
5. A comparison that could show the proposed connection is absent.

These are mapping requirements. The subsequent [working equation](interaction-equation.md) separates the revised structure from optional response laws; neither the audit nor that draft introduced a new simulation. The PDF's broader social and clinical claims have not been validated by these numerical studies.

## Evidence and preservation

This is an audit of source, protocols, reports and existing verification records. It does **not** claim a fresh rerun or independent recomputation of every saved field.

- Main repository inspected at [19b6da5d](https://github.com/NousVolition/My-Sources-Project-ALL/tree/19b6da5d8496173e9e2068842f32eaa9d9126fb9).
- Companion repository inspected at [0ae1701c](https://github.com/NousVolition/Nous-Volition/tree/0ae1701cfe7570fb54c29d689a69650172c55927).
- [Archived earlier navigation](archive/not-yet-connected/previous-results-index.md) preserves the report links and historical publication statuses.
- No numerical source, data, checkpoint, experiment matrix or scientific acceptance threshold is changed by this reorganization.

