# History prediction at increased Reynolds number

**Acquisition in progress. There is no completed independent prediction result yet.**

This experiment lowers viscosity from 0.025 to 0.01 while keeping the forcing power, initial speed, localized disturbances and prediction comparisons used in the [completed stationary study](../stationary-history/README.md). It targets a measured Taylor Reynolds number above 50; the achieved value will be reported from the new recordings. The earlier stationary experiment measured about 29.2 and did not establish a history advantage.

The prospective protocol uses 128 cubed grids, with 160 cubed controls, and fresh seeds for 24 training, 8 validation and 32 independent test flows. Four times and two disturbance strengths supply eight dependent targets per flow. The uncertainty calculation treats each independent flow as one sampling unit.

Two separate pilot initial fields are integrated through time 44 on both grids. A matched impulse pilot compares the original and halved step, then the finer grid at the same smaller step. Numerical and pilot gates must pass before calibration. Calibration stationarity must pass before fitting and freezing models. The test jobs then consume the immutable frozen artifact and verify its source, calibration and pilot hashes before integration. GitHub dependencies establish the order; no held-out targets enter model selection.

The declared test-resolution controls use four preselected test flows at the final intervention, time 44. Both restart at time 42 from the same base field, recomputing the history and response. This checks short-time spatial and time-step sensitivity, alongside separate full-duration grid pilots; it is not a continuum extrapolation of every realization.

History is compared with richer and global current measurements, wrong-run histories, and training-only averages by strength and intervention time. Support requires improvement against every comparator with a paired 95 percent run-bootstrap interval and passing all numerical and stationarity gates. A benefit from partial-observation history would not add a new memory force to Navier–Stokes.

## Runnable acquisition

Use Python 3.12 and install [requirements.txt](requirements.txt). Run `python -m pytest test_solver.py test_design.py -q`. The [protocol](protocol.json) fixes the settings, independent splits, numerical gates and prediction design. [cloud_stage.py](cloud_stage.py) drives the archived pilot, calibration, freeze, confirmation, refinement and evaluation stages. The GitHub workflow runs those stages in order and stops dependent stages if a gate fails. Completed observations, plots, uncertainty and provenance will be published in [Streams and Rocks](https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization), linked to this source folder. Pending integrations and proposed experiments remain distinct from completed results.

## Physical basis

The solver integrates three dimensional incompressible Navier–Stokes in a periodic cube, with a strict two-thirds Fourier cutoff, double precision RK4 and constant-power large-scale negative damping. The difference-energy budget includes shear production, viscous dissipation and the change in forcing work between copies. Reynolds number and cutoff times Kolmogorov length are measurements, not proof of an inertial range. These diagnostics and independent resolution checks follow the concerns discussed in [forced isotropic turbulence](https://arxiv.org/abs/1406.6317) and [finite spatial and temporal resolution](https://doi.org/10.1103/PhysRevFluids.3.064603).
