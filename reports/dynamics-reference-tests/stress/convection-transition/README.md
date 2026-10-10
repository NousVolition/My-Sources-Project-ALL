# Heating can switch a fluid disturbance from decay to growth

**Completed synthetic fluid benchmark, with full inertia.** Increasing the maintained bottom-to-top temperature difference changes the outcome in this ideal layer. At 95% of the theoretical critical heating, the initial velocity mode falls to **0.1042 times** its starting amplitude after six thermal diffusion times; at 105%, it reaches **9.078 times** its starting amplitude. These are numerical model results, not water or pond measurements.

![Measured sign change and numerical control](figures/transition.png)

## What was completed

79 full nonlinear 2D Boussinesq runs, including controls and selected repeats; 12 separate linear modal full-inertia/Stokes comparisons; **20/20 declared final checks** and **15/15 independent tests** passed. Final matrix elapsed time: 72.48 seconds on the recorded machine; this excludes the rejected initial attempt, test execution and report generation. The complete 79-run matrix was not independently repeated. Selected N=32, Ra/Ra_c=1.05 trajectories were repeated exactly at two steps; the refined comparison includes all recorded arrays.

| Heating / critical heating | Rayleigh number | Refined measured growth rate | Analytic infinitesimal rate | Final / initial mode amplitude |
| --- | ---: | ---: | ---: | ---: |
| 0.80 | 526.009 | -1.56294480 | -1.56294279 | 8.01269e-05 |
| 0.95 | 624.636 | -0.37485902 | -0.37485595 | 0.104153 |
| 1.00 | 657.511 | -0.00000348 | +0.00000000 | 0.999979 |
| 1.05 | 690.387 | +0.36559203 | +0.36559595 | 9.07773 |
| 1.20 | 789.014 | +1.41300283 | +1.41300829 | 5037.35 |

Growth is the slope of log amplitude of the critical vertical-velocity Fourier mode over t=1…6. Energy would grow at twice the amplitude rate for a pure growing mode; those rates are not interchangeable. Each refined row averages four independently seeded mixtures with matched seed across heating conditions. The four horizontal phases nearly collapse onto the same slow linear mode after the transient. Their t intervals in the raw summary are therefore extremely narrow, mainly roundoff-scale, and **are not useful bounds on numerical error or real-world uncertainty**. Step/grid changes, analytic discrepancies and modeling limitations are reported separately.

## Physical mechanism and the actual equations

Rising fluid carries temperature perturbations; buoyancy changes its vertical acceleration. With sufficient heating from below, that feedback beats diffusion and viscous loss. This is classical Rayleigh–Bénard convection. [Fabre's university lecture notes](https://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md) describe this mechanism and the stress-free neutral curve; [Ó Náraigh's course notes](https://maths.ucd.ie/~onaraigh/acm40740/acm_40740_jan2016_v1.pdf) provide further theory. These sources supply theoretical context, not experimental data analyzed here.

The model is u_t + u·grad(u) = −grad(p) + Pr Lap(u) + Pr Ra theta e_z; theta_t + u·grad(theta) = Lap(theta) + w; div(u)=0. The conductive state is motionless with temperature decreasing linearly upward. Here Ra=g alpha DeltaT H^3/(nu kappa), Pr=nu/kappa=1. Length is scaled by H, time by H^2/kappa, velocity by kappa/H, and temperature perturbation by each run's imposed DeltaT. Changing Ra represents changing DeltaT while holding geometry, coefficients and gravity fixed. No actual water properties, temperature in degrees or pond dimensions were calibrated.

The horizontal period is 2 sqrt(2) H. At z=0,H the layer is impermeable, stress-free and held at fixed temperatures: w=0, d_z u=0, theta=0. These ideal walls differ from a no-slip pond bed or a free water surface. Uniform horizontal velocity is constrained to zero; no claim of decay of an unconstrained uniform translation is made.

The same dimensionless velocity perturbation is used at every Ra, so the initial dimensional velocity matches when kappa and H are fixed. The primary thermal perturbation is the same *fraction* of DeltaT, not the same absolute temperature at different DeltaT. Two additional controls divide initial theta by Ra/Ra_c to keep its absolute size relative to the critical DeltaT fixed; they recover the same growth-rate conclusion. Halving the overall seed amplitude also preserves the linear-regime result. Four weak higher spatial modes accompany the critical roll, with independently seeded horizontal phases. Zero disturbance remains exactly zero even above the threshold in these deterministic equations.

For horizontal wavenumber k and vertical mode sin(pi z), let q=k^2+pi^2 and a=k^2/q. The independently evaluated modal matrix on (W,Theta) is [[−Pr q, Pr Ra a],[1,−q]]. Its determinant changes sign at Ra=q^3/k^2. Minimization gives k_c=pi/sqrt(2), Ra_c=27 pi^4/4=657.511364480. At Pr=1, the leading eigenvalue is −q+sqrt(Ra a). The full PDE runs retain nonlinear advection and inertia; they are compared with this infinitesimal prediction, rather than generated by that formula. This study brackets the observed sign change and checks the predicted neutral point; it does not infer an exact experimental threshold from five sample settings.

## Resolution, budgets and rejected calculations

Fourier differentiation uses N=16,24,32 on the horizontally periodic domain and a reflected vertical extension of length 2H. Thus the physical half-domain has N/2 intervals in depth; N is not the number of independent vertical physical cells. Strict 2/3 filtering removes quadratic aliases. Horizontal velocity is even under vertical reflection; vertical velocity and temperature are odd. Both reflection parity and the divergence-free projection are enforced at every stage. Time advancement uses exact half diffusion steps around RK4 for advection and coupling. The split scheme is generally second order; the equal-diffusivity linear Pr=1 benchmark has commuting diffusion/coupling operators and approaches fourth order. That special result is not claimed for arbitrary nonlinear runs or Pr.

The accepted primary step is 0.005, with 0.01 and selected 0.0025 comparisons. Maximum 0.01-to-0.005 growth-rate change is 7.62e-05; maximum refined N=24-to-32 change is 7.22e-16. Separate larger-amplitude runs at theta amplitude 0.02 compare full final fields using Fourier interpolation: N=16→24 difference 1.98e-09; N=24→32 difference 9.49e-13. These smoothly resolved roll cases do not establish resolution for turbulence or different starting fields.

Kinetic energy K=mean(|u|^2)/2 obeys K'=Pr Ra mean(w theta)−Pr mean(|grad u|^2); temperature variance V=mean(theta^2)/2 obeys V'=mean(w theta)−mean(|grad theta|^2). The code integrates production with RK stages and diffusion loss through its exact substeps. It records their separate accumulated residuals, divergence and wall errors. Energy need not remain constant: the maintained thermal background supplies the growing disturbance. These budgets are for the two-dimensional nondimensional model, not total thermodynamic energy of an apparatus.

The original step 0.01 passed growth-rate checks but missed the unchanged 0.1% budget limit near neutrality: the worst relative residual was about 0.155%. Those [initial failed budget checks](data/initial-budget-checks.json) are retained. The complete finer ensemble and finer-grid controls bring the maximum refined budget residual to **0.01018%**. All 7 coarse runs failing the budget screen remain marked false in [cases.csv](data/cases.csv); passing final refinement does not relabel them as accepted. At step 0.2, a deliberate negative control turns the true positive growth rate into -0.690985: a bounded, decaying numerical trajectory can conceal a genuine physical instability.

An earlier implementation relied on analytical reflection symmetry without enforcing it each stage. Forbidden wall modes amplified roundoff and produced invalid/nonfinite results. That interrupted attempt and original source are preserved in [rejected-boundary-run](data/rejected-boundary-run/README.md), excluded from every physical conclusion. A regression test now injects a forbidden mode and verifies removal. The initial test XML also records a separate near-neutral cancellation tolerance error: the repaired comparison scales roundoff by the cancelling terms, not their nearly zero difference. These histories are not successful fluid runs.

## Does removing inertia lose real oscillations?

![Independent full-inertia and Stokes modal comparison](figures/inertia.png)

This is a separate linear modal calculation using matrix exponentials, cross-checked against an independent adaptive integrator. Removing inertia imposes W=(Ra a/q)Theta and gives Theta'=(Ra a/q−q)Theta. It preserves the linear onset threshold here, but can change the timing and erase oscillations. Matching the threshold alone does not validate a reduction.

For stable stratification Ra=−4 Ra_c and Pr=1, the full eigenvalues are −14.804407 ± 29.608813 i. The thermal disturbance oscillates while its modal envelope decays. The Stokes reduction predicts monotonic decay. Initial velocity and temperature are matched using the algebraic Stokes velocity at time zero, so the difference is not a mismatched initial velocity. Its temperature-curve error over the recorded interval falls from 99.46% at Pr=1 to 0.273% at Pr=1000. This tests selected modes; it is not a proof for arbitrary initial data. [Wang (2004)](https://doi.org/10.1002/cpa.3047) treats the infinite-Prandtl limit rigorously and identifies its initial-layer character.

## What this means for the user's fluid question

**Yes, changing a physical condition can switch a small fluid disturbance from shrinking to growing.** This is a verified example in a thermally driven Boussinesq layer. The energy source is the maintained temperature difference. A footstep could seed a disturbance in a suitable physical setting, but this study has no footstep input and does not measure how much motion would reach a pond or suspended organisms.

The existing unforced Navier–Stokes study, the fitted dynamic-q equation and the pond-vibration model remain distinct. Their equations and results have not been changed. This calculation does not establish that the existing fluid runs undergo the same transition, that a positive cubic-feedback law is supported by their data, or that microbes select a branch. In a horizontally translation-invariant layer the roll phase is a continuum; opposite signed roll amplitudes can be translations of the same pattern. We therefore do not claim exactly two isolated physical states or prove a global pitchfork basin from these runs. Long-time saturation, hysteresis, nonlinear branch continuation, real no-slip/free-surface boundaries and physical calibration remain untested.

## What comes next: test the original flow, then the pond

**PLANNED — NOT RUN — no new fluid or pond runs accompany this plan.**

The heating study establishes a decay-to-growth transition in its specified heated layer. The next question is whether a disturbance can grow in the project's original flow, and which source of energy could support growth in a pond. Neither outcome is assumed in advance.

| Next test | Primary outcome | Go/no-go decision |
| --- | --- | --- |
| 1. Qualify the original starting field | Initial-field and trajectory convergence, divergence, total kinetic-energy loss and its viscous budget. | Candidate grids N=24,48,96 and two successive step halvings are a starting screen, not a promise of resolution. Require full-field discrepancies below 1% on the comparison interval and energy-budget residual below 0.1% of initial energy plus accumulated absolute work/loss. If the start or gradients remain unresolved, refine before extending duration. The existing repaired N=16, t=0.01 checks do not clear this gate. |
| 2. Run paired unforced flows | Primary outcome G(t)=RMS(delta-u(t))/RMS(delta-u(0)); fit log G separately on predeclared intervals 0–1, 1–2 and 2–4 initial turnover times. Save full-field differences, perturbation energy, strain-to-perturbation energy transfer and viscous loss. Scalar q and visible shape are secondary diagnostics. | Accept a finite-time growth/decay classification only after the field, budget, grid, step and amplitude checks pass. Report each seed, the ensemble mean and its 95% t interval separately from numerical error. Require a claimed sign to exceed both the numerical uncertainty margin and the ensemble interval; mixed seeds or a window crossing zero are inconclusive. A viscosity-dependent sign change on these intervals would show finite-time sensitivity of this specified decaying flow, not a global bifurcation or permanent growth. |
| 3. Test what happens after ground motion stops | Separate energy put in during motion from energy remaining afterward. Track total wave energy, surface amplitude and phase, particle displacement, and energy lost to drag. A rising local surface height alone is not evidence of instability: kinetic and surface-potential energy exchange during an ordinary oscillation. | In the current resting, passive, positive-drag pond model, total wave energy must not increase after bed motion ceases, beyond verified numerical error. This is the expected negative control. Continued forcing can amplify waves by supplying energy. Never insert negative damping just to obtain a desired growth curve. A positive post-pulse energy trend must first be treated as a model or numerical discrepancy, not as new pond physics. |
| 4. Add a physical energy source only when specified | A converged perturbation growth rate, its identified energy source and dependence on a single controlled physical condition. For a later transport study, compare passive dye with inert bacterial-sized particles before adding individually justified settling, motility or attachment terms. | Synthetic inputs support conditional model predictions only. A footstep-to-pond or organism-response claim needs measured forcing and independent validation. Growth of fluid motion alone does not establish a biological effect, microbial choice, or the behavior of trillions of organisms. |

### 1. Qualify the original starting field

Use the separately repaired hug-ns-corrected-v1 solver and one documented smooth starting field. First check a finite-Fourier-mode analytic decay control; then qualify the repaired hug starting field. Construct the same continuous field on each grid before projection and compare retained common modes, energy, gradients and spectral tails. Do not compare different grid-dependent starts as if they were the same experiment.

### 2. Run paired unforced flows

For each viscosity, evolve a reference U and an otherwise identical U+delta-u. Proposed positive viscosities are 0.005, 0.01 and 0.02 in the repaired hug model's units; hold geometry and initial U fixed. Record Re=U_ref L/nu with the same measured initial RMS speed U_ref and fixed L=6. Use eight fixed independent perturbation seeds (510–517), divergence-free and zero-mean, with the same physical low-wavenumber spectrum on every grid. Start at RMS size 1e-5 U_ref; halve it for any candidate sign-change comparison.

### 3. Test what happens after ground motion stops

Start with the existing synthetic pond model and its positive drag. Compare a zero-input control, one finite bed-motion pulse, repeated pulses, and a matched repeated-pulse case stopped at a declared time. Keep pulse strength, basin and drag matched. Retain an independent analytic or matrix-exponential reference and halve the time step twice. In a later spatial model, also refine grid and retained modes.

### 4. Add a physical energy source only when specified

If the question concerns an actual pond, measure or explicitly assume its depth and shape, bed/bank vibration, damping, existing current and temperature profile. A background current can transfer stored energy; maintained current or heating supplies continuing energy. Model these as separate, documented conditions with suitable bed and free-surface boundaries. Compare disturbed and undisturbed trajectories within each condition.

### Where the energy could come from

Unforced does not mean motionless or unable to amplify a disturbance temporarily. In periodic incompressible Navier–Stokes with the same positive viscosity in each pair, let delta-u be the difference from reference flow U and E_delta=integral(|delta-u|^2)/2. Then E_delta' = -integral(delta-u_i delta-u_j partial_j U_i) - nu integral(|grad delta-u|^2). The first term can transfer energy from the existing flow into the difference even while each trajectory's total kinetic energy decreases. This identity assumes periodic boundaries and matched forcing (zero here); extra boundary or forcing terms are needed in a real pond.

### What the present pond equations predict

For the existing single-mode pond, surface height is A cos(kx), velocity is U sin(kx), and bed displacement is B cos(kx). Its wave energy per unit transverse width is E=rho L(g A^2+H U^2)/4. The model gives E'=rho L g A B'/2-rho L H d U^2/2. Once B'=0 and d>0, E' is nonpositive. The present equations therefore predict post-pulse energy decay; another run would validate that prediction numerically, not discover a new instability.

### Controls and acceptance criteria

- Original-flow negative controls: zero disturbance stays zero; a perturbation about resting fluid matches analytic viscous decay in an invariant shear mode. Preserve seed pairing across viscosities and perturbation sizes.
- Original-flow numerical margin: target successive-refinement changes in fitted amplitude rate below 0.01 per initial turnover time. For sign classification use at least the larger of 0.01 per turnover time and three times the largest selected grid/step rate change. If this margin overlaps zero, refine or report unresolved; do not choose a favorable fitting window afterward.
- Audit perturbation energy separately from total energy. Normalize its integrated residual by initial perturbation energy plus accumulated absolute transfer and dissipation; proposed limit 0.1%. Keep divergence and initial-field/gradient resolution as separate checks. Budget agreement alone cannot certify profile accuracy.
- Pond controls: zero forcing; drag removed only as a declared conservative control (post-pulse energy then stays constant); positive-drag pulse-off decay; the same pulse schedule with different relative timing. The existing surface-accuracy limit of 0.1% and energy-budget limit of 0.5% remain separate. Refinement must additionally resolve any claimed post-pulse trend.
- Repeat one accepted fine run from a fresh process and preserve all raw outputs, parameters, seeds, versions, failures and source hashes. A reproducible deterministic run is not an independent physical sample. No future pass count or effect size is recorded before execution.

### How long to run

Start the unforced screen with four initial turnover times, L/U_ref, and evaluate the declared windows. Extend only after resolution passes and only if an unresolved timescale or late-window drift motivates it. For the existing underdamped synthetic pond (drag 0.15 per second), five amplitude-envelope e-folding times equal about 66.7 seconds after forcing stops; this also exceeds ten periods of its lowest mode. That is a planning duration, not a completed run or a real-pond timescale. Recalculate the duration for other drag, depth or modes. Benchmark wall-clock cost on one accepted case before budgeting an ensemble; no completion-time guarantee is implied.

### Decision after these tests

First complete the original-field qualification and the passive pond energy control. Proceed to an unforced paired comparison only when its numerical gates pass. Add a driven pond extension only with an explicit, independently justified energy source. Decay, temporary growth, mixed results and unresolved comparisons are all reportable outcomes.

[Machine-readable next-study plan](next-study-plan.json) · [Repaired original solver and its unresolved initial-field checks](https://github.com/NousVolition/My-Sources-Project-ALL/blob/4a6f46b9f082ec41b94400796906effa44508c9d/math/imports/hug-ns/corrected_v1/README.md) · [Existing synthetic pond model](../pond_vibration.py) · [Existing pond assumptions and omitted processes](../pond_vibration_plan.json)

## Reproduce and inspect

Use the pinned dependencies in the parent dynamics package. Run from this directory:

```text
python -m pytest test_convection.py -q --junitxml=data/independent-tests.xml
python run_study.py
python make_report.py
python verify.py
```

Then regenerate the parent package manifest with its `make_report.py` and run its `verify_package.py`. The historical rejected artifacts are preserved, not recreated by a successful rerun. Keep copies before reproducing because current outputs are overwritten.

[Declared plan and refinement rationale](plan.json) · [All case metrics and individual pass/fail](data/cases.csv) · [Summary and environment](data/summary.json) · [Checks](data/checks.json) · [Independent test record](data/independent-tests.xml) · [Delivery audit](data/delivery-audit.json) · [Recorded column names](data/record_columns.json) · [Modal comparison metrics](data/inertia.json).

Every PDE NPZ contains the complete recorded diagnostic time series plus initial, midpoint and final fields (u,w,theta) on the reflected grid. It does **not** contain every internal integration stage. Modal NPZ files contain their time grid, full-state matrix-exponential solution and reduced solution. Sources and the plan are fingerprinted in the summary. Parent manifests cover all delivered files, including clearly rejected historical evidence.
