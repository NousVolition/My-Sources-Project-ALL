# Expanded dynamics stress study

This is an extension of the neighboring baseline package. Open `report.html` for measured results. Read `stress_plan.json` for the parameter values and numerical acceptance criteria fixed before the run. Simulated time is dimensionless, not real elapsed seconds or a calibrated waterwheel clock.

## Completed scope

The separate [heated-fluid transition study](convection-transition/README.md) now tests whether changing maintained heating switches a small disturbance from decay to growth. It includes full-inertia Boussinesq calculations, numerical refinement, rejected controls and an independent comparison with the inertia-free limit. [Illustrated results](convection-transition/report.html). This is a different, thermally driven boundary-value problem from the original unforced fluid and fitted dynamic-q studies.

The complete configuration consists of 24 long Lorenz runs (eight independently seeded starts × three paired step sizes), 120 finite parameter-scan trajectories (ten rho values × six starts × two steps), 16 prescribed-rotation transport cases (four grids × four steps), three Fourier-mode time-stability controls and an independent adaptive-solver trajectory comparison. A follow-up accuracy study adds 144 discrete-map method/step/duration combinations and four directly stepped audits. There are 22 independent implementation/benchmark tests in the combined extension. Actual completion, failures and elapsed time are recorded in `data/summary.json`, `data/checks.json`, `data/oscillation-summary.json` and the test XML; this configuration description alone is not an execution record.

The first independent numerical derivative check failed marginally; its original XML is retained. The diagnostic interval was refined with the tolerance unchanged. The refined test verifies the expected reduction in derivative error. This was a reference-check resolution issue, not a change to the transport equation or a relaxed acceptance criterion.

## How to decide how long to run

### Added: outdoor pond response to ground movement

The user specified ground vibrations acting on an outdoor puddle or pond, and clarified interest in both many footsteps adding together and many tiny organisms/particles being affected. These are separate roles: vibration sources supply energy through the ground, while suspended particles may be transported by the resulting water motion. The added model does not represent biological sensing, damage or collective agency.

The concert connection is empirical context: [Caltech's 2024 account](https://www.caltech.edu/about/news/swifities-shake-it-off-and-help-seismologists-solve-mystery-of-how-concertgoers-shake-things-up) describes crowd-generated seismic motion at SoFi. [PNSN](https://pnsn.org/blog/beast-quake-taylor-s-version-from-the-vault) says the Seattle concert did not cause or trigger an earthquake. The claim about an earthquake alarm is not verified. [SCEDC publishes concert waveforms](https://local.scedc.caltech.edu/data/concert-tremor.html); these have not been downloaded or processed into this pond simulation. [USGS seismic seiches](https://www.usgs.gov/programs/earthquake-hazards/seismic-seiches) and [linear moving-bottom wave theory](https://arxiv.org/abs/physics/0611182) support the general physical coupling, not a numerical footstep calibration.

`pond_vibration.py` solves one exact spatial mode of linear shallow-water motion in an assumed basin L=2 m, H=0.1 m, density 1000 kg/m³, gravity 9.81 m/s² and linear drag d=0.15 s⁻¹. These dimensions and damping are hypothetical. The equations eta_t+H*u_x=b_t and u_t+g*eta_x=-d*u use vertical bed displacement b=B(t)cos(kx), k=pi/L. With eta=A(t)cos(kx), u=U(t)sin(kx), the temporal equations are A'=B'-HkU and U'=gkA-dU. The bed is a zero-mean deformation mode; the idealized ends have zero horizontal flux. Localized impacts and other modes are omitted.

B(t) is amplitude*sin(pi*(t-start)/width)^6 during a pulse and zero outside. Start time is 1.037 s, avoiding automatic alignment with fixed steps. Assumed amplitudes are 10 and 100 micrometres; durations are 0.2 and 2 s. These are stress-test inputs, not typical measured footstep or rock displacements. RK4 steps 0.1,0.05,0.025,0.0125,0.00625 s give 20 cases over 20 s, plus a zero-motion control. The reference integrates the closed-form damped-wave Green function against B' with 64-point Gauss quadrature. A 32-point check and independent DOP853 integration with pulse-resolving maximum steps verify the reference.

Energy per unit width E=rho*L*(g*A²+H*U²)/4 satisfies E'=rho*L*g*A*B'/2 - rho*L*H*d*U²/2. Work and drag loss are integrated as additional state variables, with the residual reported independently from surface error. Volume change integrates eta-b across the basin and is zero analytically for this mode; numerical quadrature confirms roundoff closure, not general spatial-grid convergence. Both reference and solver states/diagnostics are retained. The 0.1% relative surface-error and 0.5% relative energy-residual budgets pass at the finest step in all four inputs. Five benchmark checks and eight independent tests passed. The sweep has not been fully repeated.

Ten times the assumed bed amplitude gives ten times the linear response. At 100 micrometres input, the after-pulse surface-mode peak is 9.03 micrometres for the 0.2 s pulse and 84.23 micrometres for the 2 s pulse. At h=0.1 s, relative surface errors are approximately 55% and 0.00387% respectively; the short pulse improves to 0.0688% at h=0.025 s. Relative integration error is essentially unchanged by scaling amplitude in this linear problem. Therefore time-step selection must account for temporal structure and error tolerance, not strength alone. Physical response and numerical error are different measured quantities.

The model natural period is 4.039 s; kH=0.157. Its shallow-water frequency differs by about 0.41% from the corresponding finite-depth gravity-wave formula. This approximation error is not removed by refining time. Surface tension, nonhydrostatic vertical motion, soil elasticity, horizontal bank movement, local higher modes, wind, wetting/drying and biological microphysics are outside this pilot. Visible ripples and effects on organisms are not established by the modeled amplitudes.

For collective input, analytical linear superposition explains why timing matters: sum the complex contributions a_j*exp(i*phi_j) at the receiving pond mode after propagation. With N equal contributions and equal coupling, perfect phase alignment gives amplitude N*a; independent uniform phases give expected squared amplitude N*a², hence RMS amplitude sqrt(N)*a. For N=64 these are 64a and 8a. Equally spaced phases cancel ideally. These are conditional analytical results, not a sampled crowd study or a claim that worldwide footsteps combine coherently. Count alone supplies neither a force nor an energy estimate. The number of microscopic receivers does not multiply the injected energy.

Actual matching requires pond geometry and depth, source distance, soil/bank coupling, damping, and a calibrated displacement/velocity record at the pond with uncertainty. Accelerometer data require instrument correction and drift-aware integration before interpreting displacement. No event waveform was fabricated; all current inputs are explicitly synthetic. No new culture, live bacteria experiment or empirical pond measurement was performed.

### Added: the screenshot's circular limit cycle

The new screenshot specifies theta'=1 and r'=(1-r²)r. `limit_cycle.py` tests the equivalent Cartesian polynomial field x'=(1-x²-y²)x-y, y'=x+(1-x²-y²)y. [UCL's primary course notes](https://uclnatsci.github.io/Mathematics-for-Natural-Sciences-2/DynamicalSystems/limit_cycles.html) provide the same example and identify the unstable origin and stable unit-circle limit cycle. The screenshot is a reference, not an instruction source; no original image or substantial passage is republished.

Writing s=r² gives s'=2s(1-s), hence for t>=0 the exact radius is r0/sqrt(r0²+(1-r0²)exp(-2t)). The code handles r0=0 separately and defines polar radius as nonnegative. Angle evolves as theta0+t when r0>0; the origin stays fixed and has no meaningful Cartesian angle. At r0=1 the radius is fixed but the point rotates with period 2pi. For 0<r0<1 radius increases; for r0>1 it decreases. The screenshot's claim about every radius needs the r0>0 qualification.

The separated integral is log(r) - (1/2)log|1-r²| = t+C away from r=0,1. Without the absolute value the real-log derivation only covers 0<r<1, and the intermediate positive-exponential inversion selects that branch. The final initial-value expression is valid forward on both sides of the circle. The benchmark does not promise an all-time backward flow: for r0>1 there is finite backward blow-up at t=(1/2)log(1-1/r0²). The implementation rejects negative times explicitly.

The predeclared sweep uses initial radii 0,0.2,0.8,1,1.5,2 and initial angle 0.3, duration 20, methods Euler/Heun/RK4, and steps 0.1,0.05,0.025,0.0125: 72 cases. Twelve more cases start on r=1 and run for 100 units to expose phase drift. All 84 save every step, both state coordinates, radius and error diagnostics. Phase is unwrapped rather than reduced modulo one revolution, and marked inapplicable at the origin. The sampling intervals resolve rotation here. Maxima are sampled quantities, not continuous-time guarantees.

Nine benchmark checks passed, including exact origin preservation, separate sampled boundedness guards, finest-step RK4 error against the exact trajectory, and a full-state repeat of one fine RK4 100-unit trajectory. Twelve independent tests validate the initial-value formula against DOP853, the polar/Cartesian identities, radial differential residual and composition law, the origin exception and cycle period, domain handling, and fourth-order convergence. The full 84-case set was not independently repeated.

Accuracy budgets are maximum radial error 0.001 and phase error 1 degree, separately from a sampled radius guard of 3. At r0=1 these are 0.1% radial error and one-degree angular error. At h=0.05 over 100 units, Euler reached radius 1.01243 with phase error 2.35459 degrees; Heun had radius error 0.06490% yet phase error 2.38555 degrees, passing the radius budget and failing phase. RK4 had radius error about 1.98e-7 and phase error 0.0002971 degrees. Stable-looking motion alone does not establish trajectory accuracy.

This is a dissipative state-space oscillator, not a physical incompressible water velocity field: divergence of the Cartesian field is 2-4r². Radius changes away from the cycle, so no radius or physical mass conservation is claimed. It tests mathematical integration; it introduces no bacteria, wet-lab result or new Navier-Stokes law. There are no random inputs, so no seed ensemble or statistical confidence interval is applicable.

### Added: the exact paired-digit sequence

The clarified input `12 12 12 22 11 22 11 22 32 23` expands to `[1,2,1,2,1,2,2,2,1,1,2,2,1,1,2,2,3,2,2,3]`. `sequence_transport.py` repeats that 20-step cycle using each digit as a multiplier of a small base interval. This choice of relative units was stated before running; literal steps of one to three model units were not tested. The full reverse and uniform steps of 1.75 base intervals are equal-work controls.

The sum of the cycle is 35. At C cycles over T=8, base interval = 8/(35C), step count = 20C, and RK4 field evaluations = 80C. The uniform interval is 8/(20C). Five cycle counts (80,160,320,640,1280) were tested on grid 64 in both prescribed speed profiles; all three schedules were also tested on grids 32 and 128 at C=1280. These are 42 completed cases. The original accuracy, conservation and sampled-boundedness budgets are unchanged.

At C=1280, the user-pattern maximum sampled profile error is about 6.026e-7 and the uniform error about 2.341e-7. Uniform steps were approximately 2.6 times more accurate at equal work in these flows. Reversing the repeating pattern made very little difference. Refinement recovered approximately fourth-order time convergence, and the user pattern met the strict accuracy budgets at C=640 and C=1280. Grid 32 retained error around 1.99 despite refined time steps. None of these outcomes claims an optimal pattern across other flows.

There are 26 passing benchmark checks and 11 independent tests covering exact decoding, invalid inputs, step sizes and counts, an independent constant-frequency RK4 amplification-product audit, and fourth-order convergence. Raw NPZ files retain every accepted time and diagnostic array, with full spatial states at each cycle boundary. Two selected finest grid-64 pattern runs repeated every intermediate state exactly. The full repeat is now completed: all 336 retained numerical arrays across 42 case files matched exactly in the same environment. The full reproduction took 41.11 seconds, including regenerating the cases, selected intermediate-state repeats, plotting and comparison. `data/sequence/reproduction.json` records the result and both source hashes. The earlier usage-limit interruption is retained as history in `data/sequence/reproduction-status.json`. `verify_sequence_reproduction.py` reruns the complete comparison. This establishes reproducibility of saved numerical contents, not empirical validity or a proof of accuracy between sampled times.

The original executed source is preserved as `data/sequence/executed_sequence_transport.py.txt`, matching the original source hash in the reproduction record. The repeated run used the current source, which differs only in plot tick formatting. The displayed figure was regenerated with clearer axis labels. All numerical arrays matched the original run exactly; equations and parameters are unchanged. Both original and repeated source hashes are recorded.

### Added: adaptive intervals and opposite orders

`adaptive_transport.py` adds 16 completed, deterministic transport comparisons and ten independent tests in `test_adaptive_transport.py`. The results, source hashes, dependency versions and exact accepted times are retained under `data/adaptive/`; the new report section links to the table and figure. These are separate from the original 39 checks and 22 independent tests.

The model is q_t + omega(t) q_theta = 0 on the periodic ring, with q(theta,0)=cos(17 theta), and m=60+q. For T=8 the two prescribed angular-speed profiles are 0.2+3.8 sin²(pi t/T) and 4-3.8 sin²(pi t/T). Both have the same net rotation, 16.8 radians, and reverse the order of slow and fast motion. The exact solution is cos(17[theta-Phi(t)]), where Phi is the integrated speed. With B(t)=t/2-T sin(2 pi t/T)/(4 pi), Phi=0.2t+3.8B or 4t-3.8B respectively. No forcing, leakage, diffusion or biological mechanism is included. This is a prescribed incompressible transport benchmark, not a full Navier-Stokes solution.

The spatial derivative uses an FFT on 32, 64 or 128 points. DOP853 is tested at relative tolerances 1e-6 and 1e-9, absolute tolerance 0.001 times the relative tolerance, maximum step 0.1. It evolves the fluctuation q to avoid letting the large mean weaken error control. Native accepted states/times and dense output every 0.005 are saved. The independent tests compare the integral to quadrature, the derivative to analytic sine, nonautonomous RK4 to an unrelated exact ODE, and both tight grid-64 adaptive runs to the analytic reference and a repeated integration. See [SciPy's primary documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) for local error controls; they are not global accuracy guarantees.

Two manual RK4 schedules use h=0.008 outside t=2..6 and h=0.001 inside, or the opposite. Both have 4,500 steps and 18,000 field evaluations. Stage times account for the changing velocity. Putting the small intervals in the fast part gave profile error about 0.000975 versus 0.0972 for the opposite placement. Reversing the speed order reverses the better schedule. Neither passes the strict 1e-5 profile target. The manual diagnostics use every accepted step; maxima are sampled, not guaranteed continuous bounds.

The tight grid-64 adaptive runs each used 690 accepted steps and reached maximum sampled profile error around 1.8e-8. Grid 128 gave comparable error. Tightening tolerance did not repair grid-32 aliasing: error remained around 1.99. All 16 cases passed separate mass and sampled boundedness guards; four tight resolved adaptive cases passed all accuracy budgets. Loose resolved adaptive cases missed the profile target despite solver success. This is an accuracy comparison, not proof of an optimal schedule or general computational speed advantage.

The plan specifies profile/amplitude errors at most 1e-5, accumulated phase error at most 0.01 degrees, mass error at most 1e-10, and a sampled fluctuation bound of 2. The latter is a finite-run guard, not a general stability proof. Phase uses the unwrapped sampled Fourier coefficient; both diagnostic sampling schedules resolve its phase increments in this test. The mass diagnostic uses the mean of q, equivalent to the mass change in 60+q without subtracting two large means. There are no random inputs and therefore no seed ensemble or statistical confidence interval. Two selected adaptive cases repeat exactly in this environment; all 16 cases were not independently rerun.

Changing calculation intervals does not change the physical velocity. Changing spatial resolution is another operation. Adaptive mesh refinement is not implemented here, and later refinement cannot uniquely reconstruct a Fourier mode already lost to aliasing. No arbitrary numeric sequence has been declared optimal. A proposed sequence needs defined units and the same error/work budgets for a meaningful comparison.

### Duration remains a separate question

There is no universal wall-clock duration that validates a dynamical model. Separate four questions:

1. **Is the method accurate?** Compare with analytic solutions and independent integration, and refine the step. Stable integration can still damp real oscillations into numerical invisibility.
2. **Is the spatial grid adequate?** Refine the grid and resolve the forcing. Mass conservation alone does not reveal aliasing.
3. **Are sampled statistics stable with duration?** Compare paired prefixes and independent starts. Here we inspect 100, 250, 500 and 1,000 model units after a 50-unit transient.
4. **Is the physical/model timescale itself long?** Close to a stability boundary, an eigenvalue's real part approaches zero. A long laboratory wait or computation may show only a transient. A local e-folding estimate is 1/|Re(lambda)|, with growth above a threshold and decay below it.

At rho_H±0.01 for the classic Lorenz parameters, the computed local e-folding time is about 3,300 units; five e-folds are about 16,500 units. This is a local linear timescale prediction, not a nonlinear experiment run for that duration. The 200-unit parameter scan and 1,000-unit rho=28 ensemble cannot certify settling at those near-critical parameters. The factor five is an illustrative observation window, not a proof of convergence or a universal stopping rule.

For the rho=28 ensemble, declared stopping criteria compare the two finest steps and the 500-unit and 1,000-unit prefixes. Growth-rate mean tolerance is 0.02 across steps and 0.015 across duration; mean-z tolerance is 0.15; positive-x occupancy tolerance across steps is 0.04. These are chosen engineering tolerances for a pilot, not universal standards. Any failure remains in the report and calls for targeted investigation instead of automatically widening tolerance or claiming convergence.

Integration cost is approximately proportional to duration × number of trajectories / step for this explicit method. A longer run with unchanged resolution does not correct a biased or underresolved calculation. The report provides measured wall-clock time from this machine; scaling estimates can change with load and implementation.

## Long ensemble and uncertainty

The classic Lorenz field uses sigma=10, beta=8/3, rho=28. It reuses the independently tested baseline `core.py`. The seeds are 101,202,303,404,505,606,707,808 and initial points are sampled independently from [-0.5,0.5]×[-0.5,0.5]×[20,30]. Every seed is paired across steps 0.005,0.0025,0.00125. Each run advances the state through a 50-unit transient, initializes the tangent basis to identity, then measures for 1,000 units. QR reorthogonalization occurs every 0.1 unit and tangent rates are aggregated in 50-unit blocks.

There is no separate discarded tangent-basis transient, so the finite-time exponents include a basis-alignment contribution of order inverse duration. The near-zero exponent and duration prefixes should be interpreted with that in mind. The sum of the measured exponents is checked against the exact divergence, -41/3. State bounds use states sampled at QR times; they are finite-computation guards, not continuous-time global boundedness proofs.

Approximate 95% Student t intervals describe between-start spread in eight finite-time estimates. The two steps and overlapping prefixes are paired, not independent samples. The 50-unit blocks are retained but not treated as independent replicates. The intervals do not account for all integration bias, tangent-basis alignment bias, distributional assumptions or global ergodic questions. Positive growth at rho=28 supports the expected sensitive behavior but is not a new mathematical theorem.

Raw NPZ files retain 10,000 states sampled at 0.1-unit intervals, 20 tangent-rate blocks, the spectrum, the initial point and provenance. They do not retain every RK4 substep. Prefix statistics are computed from those saved arrays. One selected finest-step run is repeated in `audit_stress.py`; exact array comparison is recorded separately. All 24 extended runs are not repeated as part of that audit.

## Parameter scan

At each rho, two starts are within 0.001 of the two nonzero equilibria when they exist (or near the origin below rho=1). Four further starts are drawn from a broad box with an explicitly fixed seed. Each case runs for 200 units at both steps. The code is vectorized across trajectories, and the batch equations and selected trajectories are independently checked against single-model formulas and an adaptive DOP853 solution.

A sampled terminal window t=180…200 is classified as near an equilibrium if every sample lies within absolute distance 0.001 of one chosen equilibrium. Sampling is every 0.5 unit. This is a finite-window observational criterion: it does not establish asymptotic stability, classify all other trajectories as chaotic, or locate a basin boundary. The last 20 time units and threshold are fixed in the plan. The maximum state guard during this scan is checked at each full RK4 step.

Local Jacobian eigenvalues and their e-folding times are saved in a separate CSV. No unstable periodic orbit, homoclinic orbit, continuation diagram or global coexistence boundary was computed. Approximate reference thresholds from the original screenshots have not become empirically validated thresholds through this scan.

## Transport stress and negative controls

To isolate resolution from nonlinear feedback, this extension uses prescribed angular speed omega=2, not a dynamically coupled wheel. The periodic equation is

    m_t = Q(theta) - m - 2*m_theta
    Q = 60 + 28*cos(theta) + 12*cos(9*theta) + 6*cos(17*theta)
    m(theta,0) = 60.

The source has a global lower bound 60-28-12-6=14, so it is nonnegative. The mean remains 60 exactly. Each positive Fourier coefficient obeys c_k'=q_k-(1+2ik)c_k, c_k(0)=0; therefore c_k=q_k[1-exp(-(1+2ik)t)]/(1+2ik), where q_k is half the cosine amplitude. The real mass density is 60+2Re(sum c_k exp(ik theta)). This independently derived exact solution is the profile benchmark, and a finite-difference PDE residual checks it without reusing the FFT derivative.

Mode 17 aliases to mode 1 on grid 16 and mode -15 on grid 32. The sampled forcing consequently drives the wrong semidiscrete equation. Smaller time steps converge to that wrong spatial solution. Grids 64 and 128 resolve all forced modes and recover time convergence. This is still a smooth manufactured example; it does not validate discontinuous inflow, real cup geometry, evaporation or biological microphysics. The original nonlinear wheel-coupling tests remain in the baseline package.

For temporal stability, a separate complex scalar Fourier mode y'=-376 i y has exact amplitude one. The RK4 amplification polynomial is R(z)=1+z+z²/2+z³/6+z⁴/24. We compare direct integration to |R(-376 i h)|^(1/h) at h=0.002,0.004,0.008. The first two steps avoid growth yet significantly damp the true oscillation; the last is unstable. Detecting those failures is the intended negative-control outcome, not validation of the inadequate steps.

## Reproduction, provenance and interrupted runs

### Follow-up: amplitude and phase are separate requirements

In response to the question about stable integration erasing a real oscillation, the package now explicitly accepts or rejects oscillatory settings against two accuracy budgets: relative amplitude error at most 0.001 (0.1%), and absolute accumulated phase error at most one degree. Numerical stability and at least eight samples per true period are also required. These are illustrative, predeclared pilot thresholds for this mode, not calibrated experimental requirements or a universal sampling guarantee. A setting that preserves amplitude can still be rejected for phase drift.

`oscillation_accuracy.py` compares RK4, implicit midpoint, fourth-order Gauss–Legendre and an exact-rotation oracle for y'=-376 i y. It evaluates 12 steps and three horizons (1,10,100), or 144 combinations, from each integrator's exact discrete amplification map. This avoids pretending that all 144 cases were explicitly time-stepped. Four representative cases are actually advanced for 40,000 steps each; the Gauss–Legendre stages are solved independently as a two-stage linear system. Their amplitude and accumulated phase agree with the map predictions, including direct roundoff effects.

The implicit midpoint stability function is (1+z/2)/(1-z/2). The Gauss–Legendre order-four function is (1+z/2+z²/12)/(1-z/2+z²/12). For imaginary z their exact-arithmetic norms are one. The RK4 norm is computed through log1p of the identity |R(-iw)|²=1-w⁶/72+w⁸/576 to avoid rounding very small attenuation to zero. Phase is accumulated rather than reduced modulo one cycle. Coarse cases with inadequate sampling are rejected without inventing a uniquely resolved phase. The exact-rotation oracle has no discretization error for this special constant-frequency mode, but its output is still subject to the chosen sampling rule.

At h=0.00025 and T=10, implicit midpoint preserves amplitude yet has about 158.42° phase error; Gauss–Legendre order four has about 0.02335° phase error. For T=100, the largest tested steps meeting both budgets are 0.000125 for RK4, 0.00025 for Gauss–Legendre order four, and 0.00000390625 for implicit midpoint. They are the largest accepted values in a discrete tested set, not exact optimal step bounds. At RK4's selected step, predicted amplitude error is about 0.00599% and phase error about 0.08753°.

The change is to validation and accept/reject rules. No blanket integrator replacement is applied to Lorenz or the coupled wheel. An implicit method may cost more per step, and an exact mode update is not a general nonlinear solver. A physical extension needs its own resolved frequencies, desired error budget, conservation properties and validated coupling. Extra run duration alone cannot repair a time step that fails these budgets.

### Commands

From the parent package, with its pinned dependencies installed:

```text
python stress/run_stress.py
python stress/oscillation_accuracy.py
python stress/adaptive_transport.py
python stress/sequence_transport.py
python stress/limit_cycle.py
python stress/pond_vibration.py
python -m pytest stress/test_pond_vibration.py -q --junitxml=stress/data/pond-vibration/independent-tests.xml
python -m pytest stress/test_limit_cycle.py -q --junitxml=stress/data/limit-cycle/independent-tests.xml
python stress/verify_sequence_reproduction.py
python -m pytest stress/test_sequence_transport.py -q --junitxml=stress/data/sequence/independent-tests.xml
python -m pytest stress/test_adaptive_transport.py -q --junitxml=stress/data/adaptive/independent-tests.xml
python -m pytest stress/test_stress.py stress/test_oscillation_accuracy.py -q --junitxml=stress/data/independent-tests.xml
python stress/audit_stress.py
python stress/make_stress_report.py
python make_report.py
python verify_package.py
```

The initial failure XML is an archival record supplied with this package; a clean successful rerun does not recreate an intentionally obsolete test. `audit_stress.py` recreates the numerical derivative diagnostic at all three intervals as well as the selected long-run repeat.

`python stress/run_stress.py --resume` reuses completed long-run checkpoints only when their source/plan fingerprint matches the current `core.py`, `run_stress.py` and `stress_plan.json`. It reruns the inexpensive diagnostics and rebuilds summary tables. Report elapsed time for a resumed run is the time of that invocation, including any checkpoint loading; per-case reused flags are recorded. The initial complete run is not to be confused with a faster resumed read.

`data/progress.json` is a convenience snapshot. Its simple remaining-time estimate averages completed runs and therefore underestimates work before the finer steps begin; the completed elapsed time and per-run costs are the authoritative records. Output files are overwritten on reproduction; keep a copy for comparisons. The package manifest is generated at parent level after the stress report is finished.

This work is local mathematical testing. It adds no empirical biology evidence, culturing procedures, invented microbial force or new Navier–Stokes claim. Physical mixing, particulate transport and ice nucleation remain the separate questions addressed in the earlier water/biology study. Primary scientific sources and textbook provenance are linked in the parent README.
