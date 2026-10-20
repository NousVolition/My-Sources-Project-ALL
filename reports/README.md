# Finished tests

A clear view of the completed test reports. Each **Sources and details** link opens the full explanation, charts, supporting files and checks.

**Finished does not mean every check passed.** The limits below are part of the results; these numerical experiments do not resolve the Navier–Stokes smoothness problem.

## Fluid tests

| Test | Result and limits | Supporting information |
| --- | --- | --- |
| Vortex surroundings | 48 runs completed. The surroundings change the response; peak estimates remain sensitive to grid size. | [Sources and details](study/final-matrix/README.md) |
| Smooth-start stress tests | 39 runs completed. Grid and time-step checks retain the unresolved late-time results. | [Sources and details](navier-stokes-stress/README.md) |
| Completed matched-stretch controls | Finished parameter and disturbance tests through time 0.40. Spatial accuracy remains unresolved. | [Sources and details](sources-and-details.md#completed-matched-stretch-controls) |
| Fresh-flow memory and erasure | 135 configurations completed. History helps predict ordinary deformation; the disturbance-response prediction does not improve. | [Sources and details](fluid-memory-erasure/README.md) |
| Marker-history prediction | Completed prediction and refinement checks. A reliable benefit across resolutions remains unproven. | [Sources and details](sources-and-details.md#organization-transport-and-response) |
| Individual fluid markers | Different starting positions produce different paths. Detailed paths fail grid and recording-time checks. | [Sources and details](matched-stretch/marker-stress/README.md) |
| Fluid organization | Completed paired-disturbance pilot. Results and numerical limits are specific to the tested models. | [Sources and details](fluid-organization-pilot/README.md) |

## Water and molecular tests

| Test | Result and limits | Supporting information |
| --- | --- | --- |
| Water, particles and freezing | Completed transport, stratification and freezing simulations. These are model studies, without new laboratory measurements. | [Sources and details](water-biology-study/README.md) |
| Material feedback into fluid | 30 paired comparisons completed. Material changes the modeled fluid response; the assumed material law still needs calibration. | [Sources and details](water-biology-study/coupled-feedback/README.md) |
| Quantum water and fluid behavior | Completed molecular controls and fluid tests. Prediction of bulk transport from quantum behavior remains unvalidated. | [Sources and details](quantum-water-pilot/README.md) |
| One heavy water molecule | 544 trajectories completed. Isotope placement and identity controls do not establish a persistent molecular leader. | [Sources and details](molecular/isotope-mass/README.md) |
| Identical water molecules | Completed positional-influence experiments and follow-ups. Conclusions depend on the stated model and observation time. | [Sources and details](water-molecule-influence/README.md) |
| Gas return selection | Completed forecast comparisons. Short-horizon improvements coexist with simple-baseline and longer-horizon failures. | [Sources and details](water-molecule-influence/gas-return-selection-extension/README.md) |

## Material and dynamics tests

| Test | Result and limits | Supporting information |
| --- | --- | --- |
| Moldable hug: material length and longer forcing | Completed original/revised-model comparisons and reciprocal forcing tests. Gap agreement improves; some imprint checks remain unresolved. | [Sources and details](clay-gap-refinement/README.md) |
| Hug unloading and zigzags | Four controls completed. Overall gap agreement improves; individual zigzags remain sensitive to the mesh. | [Sources and details](hug-molding-refinement/README.md) |
| Heating: decay to growth | Completed convection benchmark. A disturbance decays below the tested threshold and grows above it. | [Sources and details](dynamics-reference-tests/stress/convection-transition/README.md) |
| Dynamics reference tests | Completed exact-solution, accuracy and conservation benchmarks, including Lorenz and waterwheel models. | [Sources and details](dynamics-reference-tests/README.md) |
| Dynamics stress and pond motion | Completed numerical stress tests. Pond forcing is synthetic and is not calibrated to real footsteps or concerts. | [Sources and details](dynamics-reference-tests/stress/README.md) |

## Mirror-fluid tests

| Test | Result and limits | Supporting information |
| --- | --- | --- |
| Four reproduced fluid tests | Four tests reproduced, with measured spin and mirror symmetry. The supplied signed record is passive. | [Sources and details](../math/imports/hug-ns/mirror-fluid/README.md) |
| Signed mirror calibration | 20 runs completed, with reflected starts, two grids and smaller time steps. | [Sources and details](../math/imports/hug-ns/mirror-calibration/README.md) |
| Source contributions and reduced models | Completed source, driving and reduced-model comparisons. Fitted terms are separate model tests; failed predictions remain recorded. | [Sources and details](sources-and-details.md#mirror-fluid-and-reduced-models) |

## More information

[All supporting records, earlier tests and progress history](sources-and-details.md) · [Project home](../README.md)
