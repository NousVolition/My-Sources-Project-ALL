FLUID ORGANIZATION UNDER STRESS — COMPLETED PILOT

Open report.html for the complete report and all added dynamical-systems references.

Completed: 93 paired fluid configurations, 7 particle tracking configurations,
sinusoidal and triangular phase-locking controls, six linear-system controls,
seven nonlinear pendulum period controls, paired pendulum impulses,
four fast-slow relaxation controls, weak Van der Pol / Duffing controls,
parametric swing and bifurcation controls, and 22 automated tests.
The circuit has a conditional model and algebraic tests;
no component-specific circuit run is claimed. Ice is a proposed separate extension.

Primary result: past marker geometry improves overall future neighborhood
deformation prediction, but does not improve the additional deformation caused
by the specified fluid disturbance. The finer-grid checks preserve this distinction.
The fast-water stress case fails the marker-path resolution screen (0.0348 RMS
versus a 0.005 limit), although its energy-gain change is within the 2% screen.
See results/stress_extension.json for the subsequent 38-to-50-to-62 grid
comparison. Its result does not retroactively validate the original 26-grid paths.

REPRODUCE (Python 3.12 recommended; run inside this folder)
  python -m pip install -r requirements.txt
  python -m pytest test_pilot.py -q --junitxml=results/tests.xml
  python run_experiments.py --stage all
  python particles.py
  python prediction.py
  python phase_locking.py
  python dynamics.py
  python pendulum.py
  python relaxation.py
  python weak_nonlinear.py
  python parametric.py
  python bifurcations.py
  python extend_stress.py
  python analyze.py
  python verify.py

Fluid runs resume only if source/configuration/data hashes agree. Modified code
requires a fresh data directory: run_experiments.py --data new_data --stage all.
Prediction and particle scripts also accept --data. analyze.py and verify.py
rebuild/check the delivered default data/ layout. Retain LF source bytes when
using the exact delivered cached run hashes.

FILES
  protocol.json                 Original exploratory fluid/prediction protocol
  core.py                       Paired 3D NS, markers and energy budgets
  vendor/numerics.py             Pinned upstream Fourier operators
  run_experiments.py             Screening, regimes, independent seeds, refinement
  prediction.py                 Past-only features and seed-separated evaluation
  particles.py                  Inertial fog/smoke approximation and gravity control
  phase_locking.py               Sinusoidal phase-locking experiment
  dynamics.py                    Triangular response and linear-system controls
  pendulum.py                    Nonlinear periods and controlled impulse response
  relaxation.py                  Fast-slow Van der Pol cycles and phase response
  weak_nonlinear.py              Weak Van der Pol attraction and conservative Duffing
  parametric.py                 Periodic swing forcing, rest and seeded amplification
  bifurcations.py               Saddle-node ghost passage and Hopf cycle controls
  extend_stress.py               Targeted 50-grid and 62-grid stress continuation
  test_pilot.py                  Analytic, invariance and numerical tests
  data/                         Saved states, trajectories and diagnostic records
  results/parameter_regimes.csv  Regime means and run-bootstrap intervals
  results/all_runs.csv           Every fluid configuration and quantitative results
  results/numerical_controls.csv Spatial and temporal sensitivity measurements
  results/particle_results.csv   Particle trajectories summarized through time
  results/prediction*.json       Held-out scores and numerical robustness
  results/figures/               Rebuildable PNG and SVG scientific plots
  results/verification.json      Checks of actual saved numerical records
  manifest.json                 SHA-256 inventory
  references/                   Additional excerpts supplied by the user

The HTML embeds its generated scientific plots; the reference image appendix
and links to data use local companion files. Keep the full folder together.
No remote resources are required to read the report.
