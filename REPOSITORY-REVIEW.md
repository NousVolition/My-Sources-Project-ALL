# Review of both repositories

Review date: October 9, 2026. [Project home](README.md) · [Current status](STATUS.md) · [Streams and Rocks](https://github.com/NousVolition/Nous-Volition)

## What was checked and repaired

The review began from My-Sources-Project-ALL commit `068046b739d69dfd4288552e20b798eb08b0339c` and Nous-Volition commit `7f0b200658d162450346b92d218cc6a5bbe591d0`. The initial 2,261 tracked files were retrieved or matched to their Git blob hashes. During review, Nous-Volition advanced to `0efb14238301188d4babbe253c92bc3767fcc947`; its 19 added/changed files were also retrieved and hash-verified, giving 2,271 files in the combined publication base. The cleanup preserves that new switching/input-history addendum. My-Sources then advanced to `8f8ec40399c3d007d86231c47c5ee3839d034f0b`, adding the [separate 39-run smooth-start stress suite](reports/navier-stokes-stress/README.md). Its 143 added/changed files were retrieved and hash-verified, bringing the combined publication base to 2,413 files. This later suite was preserved and included in syntax/link checks; its trajectories were not independently replayed in this review. The review checked Python and JSON parsing, local document file links, study manifests, landing pages, published methods/results, and the existing unit tests. It is not a line-by-line scientific validation of every experiment, a check of every external website, or a replay of the completed ensembles.

| Finding | Repair or disposition |
| --- | --- |
| Nous-Volition's CI only printed placeholder messages | Replaced with the existing 84 unit tests and checks of source syntax, JSON, local Markdown file links and all 196 entries in the three packaged study manifests |
| Project status and math navigation described earlier imported tables as the latest work | Added dated current status and links to the two active study indexes; retained the earlier history and test totals in their original scope |
| The report index called its updated matched-stretch count an October 8 snapshot | Corrected the description of the October 9 completion counts |
| California-paper link did not match the actual stored filename | Corrected the link; original document bytes and filename preserved |
| Broken character encoding in two Markdown guides and the motion viewer's embedded selector labels | Repaired the affected punctuation and superscript characters |
| Local Nous-Volition mirror was older than public main | Refreshed the mirror from the hash-verified public files while preserving local extras |

The My-Sources consolidated suite passed **290 tests**, with **three strict expected failures**. Those expected failures reproduce known defects in the preserved supplied hug-ns implementation: the spatial gate's corner, a retained cutoff alias when the grid size is divisible by three, and rounding that can enlarge the estimated time step. They are not three passing numerical controls. [Detailed defect evidence and repair requirements](math/imports/hug-ns/REVIEW.md#findings-that-need-correction-before-relying-on-larger-runs).

Nous-Volition passed **59 dynamics/fractal tests**, **12 social-organization tests**, and **13 recovery/entrainment tests**. Its published study manifests matched all 196 referenced files. The resisting-position source and saved records were included in the parsing/link review; its full 24,000-round reproduction was not rerun. These scopes are explicit so a green badge is not mistaken for verification of every result.

The isotope study's `manifest_sha256.json` describes its frozen original delivery. The 792 absent NPZ files belong to the separately distributed raw-data archive, and the GitHub README was added after that manifest. These are documented publication differences, not deleted trajectories. Use the [release and raw-data instructions](reports/molecular/isotope-mass/README.md), rather than interpreting that historical manifest as a current whole-repository inventory.

## Connections supported by the existing work

The studies share ways to ask and test questions. Their variables, equations, driving, and evidence remain specific to each model.

### 1. Identity, position and interaction rules

The [nine-SIMS study](https://github.com/NousVolition/Nous-Volition/blob/main/studies/one-resisting-position/README.md) distinguishes a held color from active opposition and compares hub and leaf positions. The [recovery/entrainment study](https://github.com/NousVolition/Nous-Volition/blob/main/studies/sims-recovery-entrainment/README.md) separately swaps occupants and changes where the imposed rhythm enters. The [molecular isotope experiment](reports/molecular/isotope-mass/README.md) uses placement and label controls; the [fluid marker work](reports/matched-stretch/marker-stress/separation-roles.md) follows which label has the largest measured separation.

The common question is whether an effect follows an identity, a location, or a specified interaction. A largest value, copied color, or highly connected position alone does not identify a permanent leader. Each study measures a different response.

### 2. Geometry matters through the rule that reads it

In [the SIMS fractal comparison](https://github.com/NousVolition/Nous-Volition/blob/main/studies/dynamics-fractals-sims/README.md#sims-experiment-geometry-versus-connections), the Koch boundary and matched ring have identical adjacency lists and identical paired outcomes. Their different drawings do not enter the update rule. Menger changes the graph as well as the geometry.

In [the fluid-marker relationship comparison](reports/matched-stretch/relationship-fractal-study/README.md), distances enter edge weights and an added attraction/repulsion drift. Geometry can therefore affect that model's motion. Its states and added motion do not feed back into the saved fluid. The two experiments expose different mechanisms; neither isolates a universal benefit of fractal dimension.

### 3. Symmetry, signed variables and imposed changes

The [C/D source controls](math/imports/hug-ns/mirror-fluid/source-attribution/README.md) compare reflected sources, equal and unequal weights, and each source alone. Equal signed contributions cancel in the tested setting. The [mirror-fluid definitions](math/imports/hug-ns/mirror-fluid/MODEL.md) distinguish unsigned E, signed projections, the passive q record and the explicit mirror-averaging intervention.

This connects to label-permutation controls in SIMS: a transformation that the rules ignore should preserve the outcome, while a changed rule or input is a separate intervention. Mirror averaging actively replaces a field. An assigned opposing seat actively changes its update rule. Neither intervention, by itself, establishes spontaneous symmetry breaking.

The [dynamic-q model](math/imports/hug-ns/mirror-fluid/dynamic-q/README.md#branch-selection) also supplies a precise link to the [bifurcation examples](https://github.com/NousVolition/Nous-Volition/blob/main/studies/dynamics-fractals-sims/README.md#batch-2-tipping-points-and-multiple-stable-states). For constant inputs and nonnegative gamma and g, the fitted drift `-gamma*q - g*q^3 + input` is nonincreasing and cannot have two isolated stable branches. The illustrated pitchfork `x' = r*x - x^3` can have two stable branches when r is positive. The different sign of the linear term changes the mechanism; a cubic term alone does not make the two models equivalent.

### 4. First arrival is different from persistence

The resisting-position experiment records first consensus separately from endpoint unanimity. Recovery experiments separate threshold attainment from sustained recovery; oscillator examples distinguish a center, a decaying transient and an attracting cycle. In the fluid [central-tube analysis](reports/matched-stretch/central-response/README.md), the identity of the maximum changes and stretching changes sign. A rise followed by a fall does not yet establish repeated breathing or a stable cycle. The mathematical examples help specify the measurement needed; they do not supply missing fluid evidence.

### 5. A measured input is different from a forecast

The [dynamic-q tests](math/imports/hug-ns/mirror-fluid/dynamic-q/README.md) improve conditional predictions when future measured source histories are supplied, but their frozen-input forecast is worse than simple decay; the cubic term adds no benefit. The [17 driven-fluid comparisons](math/imports/hug-ns/mirror-fluid/model-layers/FLUID.md) reject the two tested scalar-to-fluid mappings under the declared prediction screen. The [marker-history pilot](reports/matched-stretch/history-prediction/README.md) uses whole-run holdouts and causal inputs, but its small prediction gain is not established across grid refinement.

Together these results identify the missing step: a proposed connection must specify its measured variable, inputs available at prediction time, independent held-out runs, numerical accuracy and a baseline it improves upon. Fitted coefficients from these experiments are not established physical constants.

## Remaining work, in order

1. **Finish and report the authorized controls.** At this review's publication snapshot, matched-stretch has 28/48 verified completed results and the separate vortex matrix has 46/48. Follow their [matched-stretch](reports/matched-stretch/README.md) and [vortex](reports/study/README.md) indexes for subsequent batches. Keep the prescribed horizon of 0.40 and the existing starting fields, methods and per-run parameters.
2. **Resolve the numerical interpretation before claiming physical concentration or instability.** The [128-grid timestep control](reports/matched-stretch/completed-timestep128/README.md) fails its curve screen from t=0.31. The original raw strain is nonsmooth across periodic joins and initial maxima lie outside the central tube. The [160-grid departure comparison](reports/study/completed-departure160-base/README.md) still shows material global-peak grid dependence despite much closer central-region measurements. Small energy-budget residuals or seed-to-seed agreement do not establish spatial convergence. Report finite-time perturbation rates as such.
3. **Version any repaired solver separately from its historical results.** The three imported hug-ns defects above remain explicit. Repairing the gate requires a documented choice of a smooth construction; correcting the cutoff and step scheduler changes the numerical problem being replayed. These repairs need a distinct version and their existing regression checks before new results can be compared. Active numerical sources and archived intake evidence were preserved during this review.
4. **Use the completed prediction failures to constrain the model.** A quantitative mapping from scalar state to fluid observables remains unestablished. Repeating a failed mapping or fitting the passive q record to its own defining equation would not supply new evidence. Any new equation or experiment requires a stated prediction and comparison design before execution.
5. **Recover missing inputs where they block exact work.** The [SIMS atlas](https://github.com/NousVolition/Nous-Volition/blob/main/studies/sims-recovery-entrainment/README.md#files-and-scope) still lacks the heteroclinic exercise's governing equations. The [dynamics study](https://github.com/NousVolition/Nous-Volition/blob/main/studies/dynamics-fractals-sims/README.md#batch-3-cycles-exclusion-arguments-and-nonlinear-oscillators) lacks the annulus vector field and exact equations for the cropped Figure 8.1.7. The [paper guide](papers/README.md) identifies missing supporting documents. Labelled illustrations are already separate from reconstructions of those missing sources.

This review adds repository checks, repairs navigation/text, and connects already published evidence. It adds no numerical matrix, longer trajectory, new physical coefficient, or new empirical claim.
