# Stationary turbulence history study design

Historical planning memo retained from the initial screening stage. The operative specification is protocol.json; model-freeze receipts identify the version used for confirmation. This memo is not a current completion report.

Original status: pilot screening in progress. No confirmatory test runs have started.

The next experiment uses independently seeded periodic three dimensional
incompressible Navier–Stokes fields with sustained power input in Fourier modes
0 < |m| <= 2. Each disturbed copy uses the same constant power control law,
evaluated on its own state. The perturbation budget includes the resulting
forcing difference. This is an idealized forced fluid, in dimensionless units.

Pilot candidates are being screened for stationarity, Reynolds number and
resolution before fixing viscosity and mesh. The earlier transient study had
an average isotropic-approximation Taylor Reynolds number near 20; a higher
number alone will not establish developed high Reynolds number turbulence.

The planned prediction target remains the time integral of normalized difference
energy over 1.5 time units after an impulse. Candidate strengths are 0.25 and
0.75 times baseline RMS velocity. Several interventions can be launched from
separate copies along one baseline run; all observations from a run will stay
together in the same split and uncertainty will resample whole runs.

Compare global current measurements, richer current measurements, matched causal
histories, wrong-run histories, and a training-only mean. Models, preprocessing,
feature definitions, seed splits, numerical gates, hypothesis tests and analysis
code will be frozen before independent test integrations. Pilot and calibration
outcomes may guide this design; test outcomes may not.

Stationarity will be assessed using successive windows of energy, dissipation,
spectra and component statistics, accumulated turnover times and ensemble
uncertainty across independent runs. Large-scale forcing must balance mean
dissipation and a nonlinear transfer to higher modes must remain present.
Numerical checks include exact/analytical solver tests, energy budgets, a halved
step, and a finer spatial grid. Short matched restarts and long-run statistical
comparisons answer different questions in a chaotic flow and will be labeled.

This design does not assume that Navier–Stokes has extra physical memory beyond
its full current state. A history benefit with partial measurements may reveal
unobserved present information. Eulerian field history is the focus of this
extension; passive marker history and other media are separate questions.

References: [constant-power forcing and stationary DNS](https://arxiv.org/abs/1406.6317);
[low-Reynolds self-organization under this forcing](https://doi.org/10.1088/1751-8113/48/25/25FT01).
