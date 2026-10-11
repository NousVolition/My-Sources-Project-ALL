# Navier–Stokes research: hug-ns

**Working toward a solution to the Clay Navier–Stokes existence and smoothness problem.**

The central question is whether smooth three-dimensional fluid motion must remain smooth, or whether it can develop a singularity in finite time. This project develops flow constructions, follows their evolution, and checks the mechanisms that amplify or disperse vorticity.

The **hug** is the starting idea: a form that surrounds, yields, and changes shape. Its mathematical implementations and their limits are recorded in the studies below. The current fluid calculations use periodic boundaries and no external force after initialization.

## Start here

| I want to… | Open |
| --- | --- |
| Read the results and charts | [Reports by topic](reports/README.md) |
| Check what is complete and what remains | [Current status](STATUS.md#current-numerical-status--october-10-2026-utc) |
| Understand the equations or reproduce a calculation | [Mathematics and code](math/README.md) |
| Trace connections between studies | [Review of both repositories](REPOSITORY-REVIEW.md) |

**Reading on GitHub:** open a `README.md` in the **Code** tab for its formatted report. A commit or pull request’s **Files changed** tab shows editing history: orange/red lines with `−` are the old version, and green lines with `+` are the replacement.

## Start with the mathematics

| Study | What it contains |
| --- | --- |
| **[Navier–Stokes refinement stress suite](reports/navier-stokes-stress/README.md)** | 39 completed runs up to 256³, smooth vortex starts, perturbation controls, spatial and time-step checks, and full raw-data release. No continuum breakdown established. |
| **[Completed mirror-fluid tests](math/imports/hug-ns/mirror-fluid/README.md)** | Four reproduced tests, peak-spin and symmetry charts, actual results and runnable scripts. |
| **[Signed mirror calibration](math/imports/hug-ns/mirror-calibration/README.md)** | 20 completed runs with reflected starts, two grids and half-timestep controls. |
| **[Vortex stress tests: final matrix](reports/study/final-matrix/README.md)** | One central vortex with aligned, reverse, and departing surroundings; saved shapes, measurements, grid comparisons and time-step controls. |
| **[Hug-ns code and results](math/imports/hug-ns/README.md)** | The supplied solver, matched-stretch calculations, result tables and reviewed numerical controls. |
| **[Hugged-ring construction](math/notes/hug-boundary.md)** | The earlier smooth ring, its enclosing starting shape, and the calculation history through model time 0.40. |
| **[Math notes and code](math/README.md)** | Definitions, derivations, earlier experiments and instructions for using the code. |

The new vortex stress tests, the imported hug-ns calculations, and the earlier ring study have different starting fields. Their records identify which construction each result belongs to.

## See the flow

**[Reports and visuals →](reports/README.md)**

For the interactive controls, [download the project](https://github.com/NousVolition/My-Sources-Project-ALL/archive/refs/heads/main.zip), extract it, and open `reports/index.html`. The same report collection is saved locally and in this repository. GitHub’s file viewer shows HTML source; the downloaded pages run the controls in a browser.

## What the work needs to establish

The objective is a mathematical result about smoothness or breakdown for the Navier–Stokes equations. The numerical studies examine peak vorticity, accumulated maximum vorticity, core width, stretching, energy, and sensitivity to the calculation settings.

**Current stage:** saved numerical experiments and checks. A proof resolving the Clay problem has not been established. A finite recorded peak is an observation over the calculated interval, not a proof about all later times or all allowed starting flows.

- [New stress-test settings](reports/study/protocol.json)
- [Saved grid, time-step and method comparisons](reports/study/comparisons.json)
- [Hug-ns results review](math/imports/hug-ns/RESULTS-REVIEW.md)
- [Hug-ns code review](math/imports/hug-ns/REVIEW.md)
- [Report inventory and publication snapshot](reports/manifest.json)

## Water, particles and numerical validation

**[Water, biological particles and freezing](reports/water-biology-study/README.md)** — the complete three-question computational study: miscible water, separate solute and particle transport, density/temperature stratification, controlled material-property comparisons, and probabilistic ice nucleation. Includes runnable code, controls, raw outputs, uncertainty, sources and verification. No new biological measurements or cultures.

**[Dynamics reference tests](reports/dynamics-reference-tests/README.md)** and **[expanded stress tests](reports/dynamics-reference-tests/stress/README.md)** — exact-solution benchmarks, Lorenz dynamics, waterwheel conservation, spatial aliasing, oscillation amplitude and timing, alternating/adaptive steps, limit cycles, and a synthetic pond response to ground motion. The pond inputs are assumed, not calibrated footsteps or concert measurements. These studies do not establish a new fluid law or a solution of the Clay problem.

Download the repository and open `reports/dynamics-reference-tests/stress/report.html` or `reports/water-biology-study/report.html` for the illustrated offline reports. All study files are stored individually, including raw arrays; duplicate ZIP bundles and temporary environments are omitted. [Publication inventory](reports/water-dynamics-publication.json).

## Material response and molecular studies

**[Moldable hug with reciprocal Lorenz drive](reports/clay-hug/README.md)** — a filled 2-D yielding material that reshapes around rigid objects and retains a fading imprint. Includes matched one-way/two-way drive tests, work balances, numerical stress limits, refinement, raw trajectories and an interactive offline report. This is an abstract material prototype, not calibrated physical clay.


**[Quantum water and observable fluid behavior](reports/quantum-water-pilot/README.md)** — H₂O/D₂O vibrational predictions, 24 molecular runs and 131 fluid simulations/check runs, with raw data and stress tests. The finer 100 mm/s whole-flow comparisons pass the 1% grid screen; quantum prediction of bulk transport remains unvalidated.

**[Fluid organization under stress: completed paired-disturbance pilot](reports/fluid-organization-pilot/README.md)** — code, full recordings and tests for the [Streams and Rocks study](https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization). Marker-history prediction, water/air regimes, dilute fog/smoke controls, energy budgets and explicit numerical limits.

**Separate molecular study:** [One heavy water molecule](reports/molecular/isotope-mass/README.md) records 544 explicit-water MD trajectories, isotope-placement and identity controls, short-time response tests, raw-data provenance and limitations. This confined classical-water experiment is distinct from the Navier–Stokes work; it establishes neither a persistent molecular leader nor new fluid physics.

[Identical water molecules: positional influence](reports/water-molecule-influence/README.md) collects its own experiments and controls. Each study states its model, evidence and limits.

## Background: where the ideas came from

The human–AI conversations, mirrored categories, and pressure-and-memory illustrations document the development of the ideas. They provide context for the mathematical project.

- [Pressure, imprint and the soft envelope](math/notes/soft-envelope.md)
- [Norm and Entropy: mirrored halves](math/notes/norm-entropy-mirror.md)
- [Conversation study and source records](language/conversation-study/README.md)
- [Conversation categories](language/exploration/README.md)
- [Coordination with the freedom to refuse](models/decentralized_coordination/README.md): a separate nine-agent simulation and optional volunteer protocol exploring local organization, refusal, and identity versus network position. No human data or water-physics claims.

## Connections and remaining work

[Review of both repositories](REPOSITORY-REVIEW.md): related symmetry, identity/position, persistence and prediction tests; current numerical limitations; and the next unresolved decisions. [Streams and Rocks](https://github.com/NousVolition/Nous-Volition) contains the SIMS, recovery, oscillator and finite-fractal studies.

## Source, history and reproduction

- [Math source and instructions](math/README.md#use-the-code)
- [Project history](STATUS.md)
- [Working papers](papers/README.md)
- [Earlier archive](archive/README.md)
- [Original stream-function repository](https://github.com/NousVolition/Nous-Volition)

The report collection is a dated copy of saved local work. Ongoing calculations can produce newer results locally before the next publication; the snapshot date is recorded in its [index](reports/README.md).


