# Reports by topic

[Project home](../README.md) · [Current status](../STATUS.md) · [Mathematics and code](../math/README.md)

Start with a study's **README** for its question, results, charts and limitations. Open its data and verification files when you want the supporting detail.

## Fluid calculations

Published status: October 10, 2026 (UTC). Each study has its own starting field and protocol.

| Study | Published results | Start reading |
| --- | --- | --- |
| Matched-stretch controls | **29/48 complete**, through 0.40; 19 remain | [Current study page](matched-stretch/README.md) · [Latest 128-grid 0.5% seed](matched-stretch/completed-n128-halfpercent-seed101/README.md) |
| Three vortex surroundings | **48/48 complete**, through 0.40 | [Final matrix and charts](study/final-matrix/README.md) · [Study page](study/README.md) |
| Separate smooth-start stress suite | **39/39 complete** | [Results and numerical limits](navier-stokes-stress/README.md) |
| Earlier hugged ring | Saved evolution, grid and timestep comparisons | [Construction and results](../math/notes/hug-boundary.md) |

**Numerical limits:** the matched-stretch field has nonsmooth periodic joins and initial maxima outside the central tube; its spatial and timestep screens remain inadequate. The completed vortex matrix also lacks established spatial convergence. These computations do not resolve the Navier–Stokes existence and smoothness problem.

<details>
<summary>Earlier fluid controls and starting-field comparisons</summary>

- [Original 48/64/80 replay: fixed and moving cutoff measurements](support/fixed-cutoff-run/README.md).
- [Three completed 0.1% perturbation seeds](matched-stretch/completed-seeds/README.md).
- [Three completed 0.5% perturbation seeds](matched-stretch/completed-halfpercent/README.md).
- [Separate periodic starting field: five initial grids](matched-stretch/periodic-start-check/README.md). The 256-grid width check passes; the large change in starting energy is reported.
- [First completed 128-grid 0.1% seed](matched-stretch/completed-n128-seed101/README.md): saved-field checks, baseline separation and resolution limits.

</details>

## Water, biological particles and dynamics tests

| Study | Completed work and limits |
| --- | --- |
| [Water/biology three-question study](water-biology-study/README.md) | Solute and particulate transport, buoyancy, controlled material properties and stochastic freezing; code, figures, full raw outputs and uncertainty. Model predictions and published evidence are distinguished. |
| [Dynamics reference benchmarks](dynamics-reference-tests/README.md) | Exact solutions, Euler/Heun/RK4 convergence, fixed points, Lorenz dynamics and waterwheel conservation. |
| [Expanded dynamics and pond-vibration tests](dynamics-reference-tests/stress/README.md) | Long-run Lorenz statistics, aliasing, amplitude/phase errors, adaptive and alternating steps, limit cycles and a moving-bed pond pilot. Ground inputs are synthetic; no real footstep, rock or concert strength has been calibrated. |

For the illustrated offline views, download the repository and open [water/biology](water-biology-study/report.html), [baseline dynamics](dynamics-reference-tests/report.html), or [expanded stress results](dynamics-reference-tests/stress/report.html). GitHub itself displays HTML source. [Checksummed publication inventory](water-dynamics-publication.json).

## Organization, transport and response

| Question | Report |
| --- | --- |
| How does a disturbance change organized fluid? | [Fluid organization pilot](fluid-organization-pilot/README.md): paired disturbances, history targets and numerical controls |
| Does past marker geometry predict later deformation? | [Marker-history prediction](matched-stretch/history-prediction/README.md): whole-run holdouts; finer-grid comparison remains inconclusive |
| Which moving labels separate fastest? | [Marker stress and changing roles](matched-stretch/marker-stress/README.md) |
| How do peaks compare with a changing reference? | [Peak and reference viewer](visuals/peak-and-reference.html) |

## Molecular water

These studies use separate molecular or toy models. Their own reports state the applicable methods and limitations.

| Study | Report |
| --- | --- |
| Quantum water and observable fluid behavior | [Vibrational predictions, molecular controls and measured-property fluid tests](quantum-water-pilot/README.md): completed raw-data package; bulk quantum transport remains unvalidated |
| One heavy water molecule | [Isotope placement, identity controls and raw-data provenance](molecular/isotope-mass/README.md) |
| Identical water molecules | [Positional influence and follow-up experiments](water-molecule-influence/README.md) |

## Mirror-fluid and reduced models

[Four reproduced fluid tests](../math/imports/hug-ns/mirror-fluid/README.md) · [Signed calibration](../math/imports/hug-ns/mirror-calibration/README.md) · [Definitions and controls](../math/imports/hug-ns/mirror-fluid/methods-controls/README.md)

The supplied q is a passive signed record, E is unsigned, and mirror averaging is an explicit intervention. The later fitted q and cubic-feedback experiments are separate model tests; their coefficients are not established physical constants.

<details>
<summary>Mirror-fluid reports, source controls and earlier completed comparisons</summary>

### Completed mirror tests — October 8

- **[Separate C and D contributions](../math/imports/hug-ns/mirror-fluid/source-attribution/README.md):** ten new fluid controls and [chart](../math/imports/hug-ns/mirror-fluid/source-attribution/source-effects.png).

- **[Four reproduced fluid tests](../math/imports/hug-ns/mirror-fluid/README.md):** [charts](files/fluid-tests-completed.png), [HTML report](files/fluid-tests-completed.html), code and result tables.
- **[Signed mirror calibration](../math/imports/hug-ns/mirror-calibration/README.md):** 20 completed runs, [charts](files/mirror-calibration.png) and [HTML report](files/mirror-calibration.html).
- [Files and hashes for this addition](verified-mirror-publication.json).

![Reproduced fluid tests](files/fluid-tests-completed.png)

### Completed comparisons

- [Dynamic source weights and cubic term: 26 runs, four unseen starts](../math/imports/hug-ns/mirror-fluid/dynamic-q/README.md).
- [Four additional 65-grid controls for the dynamic q model](../math/imports/hug-ns/mirror-fluid/dynamic-q/GRID65.md).

- [Breathing rerun: six reproduced cases with prescribed per-step driving](../math/imports/hug-ns/mirror-fluid/breathing-rerun/README.md).

- [Breathing run followed through every step: actual-fluid animations and the changing E/D measurements](../math/imports/hug-ns/mirror-fluid/breathing-rerun/timeline/README.md).

- [Reproduced one-lean ball and breathing-fluid series](../math/imports/hug-ns/mirror-fluid/breathing-rerun/one-lean/README.md).

### Measurements and definitions

- [Corrected fluid dashboard, explicit definitions, 516 saved-field checks and completed controls](../math/imports/hug-ns/mirror-fluid/methods-controls/README.md).
- [Norms, signs, reflection, numerical methods and evidence categories](../math/imports/hug-ns/mirror-fluid/methods-controls/METHODS.md).

### Additional controls

- [Four reduced-model layers: 46 completed runs and plots](../math/imports/hug-ns/mirror-fluid/model-layers/README.md). All 17 driven fluid comparisons are complete: [prediction scores and numerical controls](../math/imports/hug-ns/mirror-fluid/model-layers/FLUID.md). The tested scalar mappings fail the unseen-run prediction check.
- [All 22 numerical controls completed](../math/imports/hug-ns/mirror-fluid/methods-controls/README.md), including grid, timestep, translation, rotation and initial-shape checks.

</details>

## Visuals and conversation studies

[Conversation study and source records](../language/conversation-study/README.md) · [Conversation report guide](../language/README.md)

**Interactive controls:** [download the project](https://github.com/NousVolition/My-Sources-Project-ALL/archive/refs/heads/main.zip), extract it, and open `reports/index.html`. GitHub's file viewer displays HTML source. Older viewers retain the snapshot date shown in their report.

<details>
<summary>Browse the interactive viewers and earlier illustrations</summary>

Earlier interactive-viewer snapshot: 2026-10-08T20:40:59.597392+00:00. That display retains its historical data. The [final matrix report](study/final-matrix/README.md) contains all 48 completed records and current numerical conclusions.

![Vortex viewer](files/adversarial-vortex-study-visual-preview.png)

| Report | Purpose |
| --- | --- |
| [Vortex shapes and movement](files/adversarial-vortex-study.html) | Rotate the 3D field and play saved times. Aligned, reverse and departing surroundings. |
| [Flow results and controls](files/change-influence-result.html) | Imported results, fixed thresholds, and what each numerical check measured. |
| [Hugged-ring calculations](files/hug-boundary-results.html) | The earlier ring calculation, grid comparisons, energy and time-step checks. |
| [Peak and reference](visuals/peak-and-reference.html) | Compare the largest spin with its changing reference. |
| [Norm, Entropy and the hug](visuals/norm-entropy-mirror.html) | The earlier mirror and timed-envelope illustration; this is a separate prototype. |
| [Flip and join](visuals/flow-edge-seam.html) | The illustrated join between flipped copies. |
| [Why the corner matters](visuals/why-the-corner-matters.html) | Explore the difference between a tall smooth peak and a corner. |
| [People in context](visuals/people-in-context.html) | The A/B/C conversation illustration. |
| [Conversation without words](visuals/conversation-without-words.html) | Inspect categories, reply order and lengths. |
| [Conversation study](../language/conversation-study/README.md) | Source-linked observations and conversation patterns. |

## Earlier experiments

- [Earlier test using a different starting flow](files/hidden-flow-results.html)

</details>

## Data, verification and history

| File or folder | How to use it |
| --- | --- |
| `README.md` | Read the formatted explanation, figures and limitations |
| `result.json`, `measurements.json`, `comparisons.json` | Inspect saved measurements and comparison values |
| `verify.py`, `plot.py`, `compare.py` | Inspect or reproduce the checks and charts |
| `verification.json`, manifests and hash lists | Check the stated verification scope and file identity |
| `completed-…` folders | Read a dated batch; use the study's main README for subsequent completions |

- [Original report inventory and source hashes](manifest.json) and [verification of that saved snapshot](numerical-progress-verification.json) describe their dated publication scope.
- [Vortex protocol](study/protocol.json), [saved comparisons](study/comparisons.json), and [run measurements](study/runs) provide the supporting records.
- [Adaptive peak-spin audit](matched-stretch/adaptive-peak/README.md) records the failed initial six-cell width screen.
- [Project history](../STATUS.md) and [review of both repositories](../REPOSITORY-REVIEW.md) retain earlier findings and dates.
- [Earlier math notes and code](../math/README.md).

Full restart arrays for the two 48-run fluid matrices remain in the local calculation workspace; the published JSON files are measurement records. Other studies state their own raw-data release arrangements.

**Seeing orange/red lines and minus signs?** You are viewing a commit or pull request diff. Those lines show the previous version. Return to the **Code** tab and open the README to read the current report.

