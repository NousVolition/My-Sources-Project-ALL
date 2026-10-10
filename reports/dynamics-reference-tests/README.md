# Dynamics reference tests

Local computational companion to the textbook screenshots supplied on 9 October 2026. Open **report.html** for the illustrated results. These are tests of stated mathematical models, not laboratory measurements. This package was created locally; it does not create or update a pull request.

An expanded follow-up is now available in **stress/report.html** with methods in **stress/README.md**. It adds 24 long runs, a parameter scan, harder transport resolution checks, and explicit oscillation amplitude/phase acceptance rules. Its results and uncertainty are separate from the original baseline measurements below. The parent archive and manifest include both studies.

## Completed tests and primary outcomes

| Reference / question | Computation | Outcome / benchmark | Controls and replication |
|---|---|---|---|
| Euler, improved Euler (Heun), RK4 | Logistic and RC integration | Maximum error over the complete sampled trajectory; order under step halving | Exact solutions; four steps 0.2, 0.1, 0.05, 0.025; independent vector problem against matrix exponential |
| Logistic growth and slope field | x'=x(1-x), six positive initial values | Approach to x=1; invariant equilibria; reciprocal u=1/x | Exact logistic solution; compare coarse Euler h=0.25, 1.5, 2.8; linear-decay negative control |
| Stable, unstable, half-stable, neutral | -x^3, x^3, x^2, 0 | Trajectory error and direction on either side of zero | Four initial values per field; exact solutions; stop before finite-time blow-up |
| Infer an equation from a sketch | Candidate x'=x(x-1) | Match equilibrium and direction pattern | Four initial values and closed form; candidate is not a unique identification |
| RC relaxation | Q'=(2-Q)/0.75 | Equilibrium 2; 95% charging time 2.246799 | Exact exponential and step refinement |
| Overdamped approximation | m x''+x'+x=0; x(0)=1, x'(0)=0 | Displacement error versus x=exp(-t); energy budget; zero crossings | Four masses 1, 0.1, 0.01, 0.001; Radau versus matrix exponential; integrated drag loss |
| Nonuniqueness | x'=real_cuberoot(x), x(0)=0 | Residual of six departing/waiting branches plus zero | Both signs, waiting times 0,1,2; independent finite-difference derivative; perturbed initial values |
| Saddle-node | x'=r-x-exp(-x) | Roots and stability; error of ±sqrt(2(r-1)) approximation | 51 positive offsets from 10^-6 to 0.1; bracketed roots; analytic classification for r≤1 |
| Cropped autocatalysis example | A+X ⇌ 2X, assumed elementary mass action | Exact logistic agreement; closed-system A+X conservation | Closed A+X versus externally maintained A; explicitly assumed rate constants |
| Lorenz local dynamics | sigma=10, beta=8/3 | Fixed-point residuals, Jacobian eigenvalues, Hopf threshold, symmetry | Nine rho values, both branches; 100 symmetry points; finite-difference Jacobian |
| Lorenz numerics | Short trajectories and tangent matrices | Fourth-order convergence; log-volume identity | Independent DOP853 reference; exact divergence; four time steps |
| Lorenz chaotic sensitivity | QR Lyapunov spectrum | Positive largest growth rate; spectrum sum | Four independent initial-condition seeds at each of two steps; paired initial conditions; 30-unit transient then 180-unit measurement |
| Waterwheel conservation / reduction | Spectral ring mass equation plus angular velocity | Total mass, positivity, higher-mode decay, projection to Lorenz beta=1 | Three grids × three steps; exact mass law; constant-speed transport; no-torque, no-forcing-asymmetry, higher-mode ablations |

Acceptance thresholds, all measured values and every pass/fail are in `data/checks.json`. Independent tests are in `test_independent.py`; their execution record is `data/independent-tests.xml`. A passing negative control means that the intended numerical failure was detected, not that the coarse integrator was accurate.

## Quantitative findings from the pilot

- Logistic observed orders: **1.0008** for Euler, **1.9845** for Heun and **3.9850** for RK4. At step 0.025, maximum errors were **0.003195**, **2.737e-5**, **4.813e-10**, respectively. These numbers apply to x(0)=0.1 and 0≤t≤4.
- Coarse Euler at h=2.8 crossed the logistic equilibrium **25 times in its final 30 samples**. The continuous positive logistic trajectory approaches equilibrium monotonically. The discrete equilibrium multiplier is 1-h=-1.8, so this is a numerical instability.
- At mass 0.1, 0.01 and 0.001, maximum displacement discrepancies from the zero-inertia approximation were **0.07093**, **0.009313** and **0.0009831** on the sampled 0≤t≤12 interval. The full model's numerical errors were about 10^-11; these larger discrepancies are reduction error. The sampled maximum can slightly underestimate the continuous-time maximum, especially during a fast initial layer.
- Seven constructed solutions of the cube-root initial-value problem share x(0)=0. Their largest analytic residual was **2.22e-16**. The numerical solver initialized at exactly zero stayed zero. The equation has infinitely many such branches; seven were checked here.
- The saddle-node root approximation had maximum relative error **0.0002357** (0.02357%) at r-1=10^-6. This is a local approximation, not an exact replacement away from the bifurcation.
- Lorenz's local Hopf threshold was **rho_H=24.7368421053** for sigma=10, beta=8/3. This test concerns local eigenvalues, not the global onset of a chaotic attractor or the criticality of the Hopf bifurcation.
- At the finer Lorenz step 0.0025, the largest finite-time Lyapunov rate averaged **0.89875 per dimensionless time unit**, with an approximate 95% t interval **[0.85747, 0.94003]** across four initial conditions. At step 0.005 it averaged **0.91553**. All eight runs had positive largest rates. The interval does not include every source of numerical or finite-time bias.
- The Lorenz phase-volume log identity had maximum absolute error **6.57e-11** over 0≤t≤0.5. The volume contraction rate is **-13.6666667**; positive trajectory growth and negative total volume growth coexist.
- The finest waterwheel run had mass-budget error **2.27e-13**, projection error **4.43e-9** versus a separate high-accuracy Lorenz beta=1 solve, and observed time order **4.0017**. Resolved-grid differences were **2.98e-13**. These errors use the chosen model units.

## Equations and interpretation

### Scalar stability and fixed-point exercises

At x=0, the derivatives f'(0) vanish for ±x^3 and x^2. Linearization alone is inconclusive. The sign of f determines attraction or repulsion: -x^3 attracts from both sides; x^3 repels; x^2 attracts from the left and repels from the right. For f=0 every starting value remains fixed. These are representative equations consistent with the sketches; the pictures do not determine unique formulas. The x^3 runs end at 0.4, before the earliest blow-up at 0.5 for |x0|=1.

For the separate exercise asking for smooth fields with specified fixed points:

- Every real number fixed: f(x)=0.
- Every integer and no other fixed points: f(x)=sin(pi*x).
- Precisely three fixed points, all stable: impossible for a continuous smooth scalar field on the real line with three isolated two-sided stable equilibria. On the interval between adjacent isolated stable zeros, f would have to be negative near the left zero and positive near the right zero, requiring another zero by continuity. The issue is a sign argument, not an inability to find a formula.
- No fixed points: f(x)=1.
- Exactly 100 fixed points: f(x)=product from k=1 to 100 of (x-k). This is an analytic construction. We do not use floating-point polynomial root-finding to certify that count.

These exercise answers are analytic reasoning, not additional counted numerical tests. For the portrait with equilibria at 0 and 1, f=x(x-1) is one consistent candidate. Above 1 it blows up in finite time; the finite sketch alone says nothing about arbitrarily late times.

For the tested logistic model, u=1/x gives u'=1-u wherever x≠0. Thus u=1+(1/x0-1)e^-t and x=1/u. This supplies the exact solution and a separate change-of-variable check. The isolated cropped instruction about 1/N does not specify its original equation; no missing equation has been inferred.

### Nonuniqueness is different from chaos

For any waiting time tau≥0 and sign s=±1,

    x(t)=0                            for t≤tau
    x(t)=s*[2(t-tau)/3]^(3/2)         for t>tau

satisfies x'=real_cuberoot(x), including at the join. The stationary solution is another possibility. The cube-root field is continuous but not locally Lipschitz at zero. A numerical result of zero cannot select a uniquely mandated physical outcome. In contrast, the Lorenz vector field is smooth and has local uniqueness; its sensitivity compares different initial states.

The perturbation illustration uses x0=10^-3, 10^-6 and 10^-9 and a fixed RK4 step 10^-4. Exact perturbed solutions are also saved. At x0=10^-9 the final error is **8.91e-6**; the increasingly short initial time scale is not uniformly resolved by that fixed step. No high-order convergence claim is made for this non-Lipschitz limit.

### Damping, energy and approximation

For m x''+b x'+k x=0, E=(m v²+k x²)/2 satisfies E'=-b v². The test integrates the dissipated energy as a third state and checks E(t)+integral(bv² dt)=E(0). Here b=k=1. Neglecting inertia gives x'=-x and x=e^-t. The full solution starts at v=0, whereas the reduced model immediately has x'=-1; an initial adjustment layer is therefore expected. A well-posed autonomous scalar ODE cannot sustain a nonconstant periodic orbit. This restriction does not apply to a multi-state or externally driven viscous system.

### Saddle-node and autocatalysis

For f=r-x-exp(-x), f'(x)=exp(-x)-1 and f''(x)=-exp(-x)<0. Its unique maximum occurs at zero and equals r-1. Hence there are no roots for r<1, one nonhyperbolic root at r=1, and two roots for r>1. The negative root is unstable and the positive root is stable. Taylor expansion yields (r-1)-x²/2+O(x³).

The reaction crop does not specify whether A is fixed. We explicitly test two assumed elementary mass-action models, with forward rate kf*A*X and backward rate kr*X² (any combinatorial convention is absorbed into kr). With A fixed, X'=kf*A*X-kr*X². With A+X=T conserved, X'=kf*T*X-(kf+kr)*X². We choose kf=1, kr=0.5, T=1, X0=0.05; the fixed-A comparison uses A=0.8. Equilibria are 2/3 and 1.6, respectively. These illustrative reaction kinetics are not measured bacterial growth rates, a culture procedure or a model of microbes choosing motion.

### Lorenz and phase-space volume

The classic equations tested are:

    x'=sigma*(y-x)
    y'=rho*x-y-x*z
    z'=x*y-beta*z

The equilibria are the origin and C±=(±sqrt(beta*(rho-1)), ±sqrt(beta*(rho-1)), rho-1) for rho>1. The nonzero-branch local Hopf threshold is sigma*(sigma+beta+3)/(sigma-beta-1), requiring a positive denominator. The symmetry sends (x,y,z) to (-x,-y,z).

The Jacobian trace is -(sigma+1+beta). We integrate the variational equation A'=J A, A(0)=I, then compare log|det A| to the exact trace integral over a short interval to avoid determinant ill-conditioning. For long growth rates we repeatedly orthogonalize tangent vectors (QR every 0.1 time unit). Four independent seeds are drawn uniformly in [-0.5,0.5]×[-0.5,0.5]×[20,30]. Each seed is paired across the two integration steps; these are four initial conditions, not eight independent physical samples. Ten-unit block rates are retained but not treated as independent measurements.

Finite-time growth rates and trajectory separation support chaotic sensitivity in this parameter setting; finite computations do not prove an asymptotic theorem. After trajectories separate, pointwise agreement between different step sizes is not a sensible long-time convergence requirement. The short-time trajectory error and long-time growth statistics are reported separately. The pictured orbit uses densely sampled data from the nearby-pair run, t=10 to 40, rather than linearly connecting sparse QR samples.

### Waterwheel: conservation first

With mass per angular coordinate m(theta,t), use the periodic continuum equations

    m_t = Q(theta) - K*m - omega*m_theta
    omega_t = -d*omega + c*a1
    a1 = (1/pi) integral m*sin(theta) dtheta

Q=q0+q1*cos(theta), q0=60, q1=28, K=1, d=c=10. Initial m=45+sin(theta)+27cos(theta)+2cos(3theta), omega=1. The inflow is nonnegative. Total mass M=integral m dtheta obeys M'=2pi*q0-K*M exactly; the computation compares against M=2pi*[60-15exp(-t)]. The positive initial mass and source are resolved in this smooth example.

Write m=m0+a1 sin(theta)+b1 cos(theta)+higher modes. Projection gives a1'=omega*b1-K*a1 and b1'=-omega*a1-K*b1+q1. For scaled time tau=Kt, x=omega/K, y=c*a1/(d*K), z=rho-c*b1/(d*K), sigma=d/K, rho=c*q1/(d*K²), the first-mode system is Lorenz with **beta=1**. Higher harmonics do not feed the sine torque in this idealized model. The unforced third-mode amplitude decays as 2exp(-t). Thus the closure is exact within these equations; it is not proof that all real cups follow this continuum model.

The test assumes constant effective inertia, linear leakage, linear drag, a continuous ring, and prescribed smooth inflow. It omits cup discreteness, overflow, dry friction and changes in inertia as water mass changes. The deliberately smooth source is a validation case, not a reconstruction of a particular physical apparatus. Since only modes 0,1,3 are populated, all three grids resolve the spectrum; this checks agreement of resolved grids, not a spatial convergence order for discontinuous inflow.

## Relationship to the fluid / biology study

These examples test numerical reasoning and reduced dynamical systems. Lorenz phase volume means a volume of possible system states. It is not the spatial volume of water, and its contraction does not contradict incompressible velocity fields satisfying divergence(u)=0. The wheel has an open physical mass budget because water enters and leaks out. Neither this suite nor the cube-root example proves existence, uniqueness, blow-up or a new law for the full three-dimensional Navier–Stokes equations.

The earlier three-question study remains in the neighboring `water-biology-study` package. No new bacteria, freezing or wet-lab experiments were performed here. Water carrying bacteria is miscible with relatively bacteria-free water. Advection, diffusion, density and thermal stratification govern mixing; cells are suspended particles that can additionally settle, attach, aggregate or swim under suitable modeled conditions. Dissolved dye and bacterial-sized particles need not disperse at the same rate. Living cells, dead cells and extracellular ice-nucleating material have distinct possible roles; ice-nucleating activity does not establish pathogenicity. These dynamical-system examples do not measure any of those effects or supply new biology-fluid forces.

## Not completed / limits

The unstable saddle cycles, their continuation, subcritical Hopf criticality, homoclinic threshold near 13.926, attractor coexistence boundaries near 24.06 and escape/return structure shown in the screenshots were **not independently computed**. A local Jacobian calculation cannot establish them. Continuation and basin studies would be separate extensions. No experimental parameters were fitted. No empirical dataset was reanalyzed in this supplementary package. No biological material was handled.

Numerical error is quantified by exact benchmarks, independent solvers and step/grid comparisons. Parameters are chosen examples with no inferred physical uncertainty. The four-seed intervals assume enough between-start independence for a rough t summary; they do not bound finite-time bias or certify a universal exponent. Exact repeatability is tested on the recorded environment; chaotic trajectories can differ across libraries or hardware even when qualitative statistics agree.

## Reproduce

Use Python 3.12 and the pinned dependencies in a fresh virtual environment. From this package directory:

```text
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-lock.txt
.venv\Scripts\python run_tests.py
.venv\Scripts\python -m pytest test_independent.py -q --junitxml=data/independent-tests.xml
.venv\Scripts\python verify_reproduction.py
.venv\Scripts\python make_report.py
.venv\Scripts\python verify_package.py
```

On macOS/Linux use `.venv/bin/python`. The numerical run and report need no network after dependencies are installed. `verify_reproduction.py` repeats the full numerical run and compares CSV bytes, decompressed NPZ arrays and numerical summary/check files. It deliberately excludes timestamps, execution durations, figures and platform metadata. `make_report.py` regenerates the readable report and checksum manifest; `verify_package.py` verifies that manifest. Reproduction overwrites generated files within this package; keep a copy to preserve an earlier run.

`core.py` contains the small integrators and model equations. `run_tests.py` writes raw CSV/NPZ results and figures. `test_independent.py` uses independent or analytic comparisons. `data/environment.json` records versions and source hashes from the completed run. `data/reproduction.json` records the repeated-run comparison. All figures are supplied as PNG and SVG. The archive is a self-contained copy of this package, excluding interpreter caches and virtual environments.

## Sources and provenance

The primary specification is the user-supplied textbook screenshots, including sections on one-dimensional flows, numerical integration, bifurcations and Lorenz equations. Page excerpts were treated as reference material. We did not infer an edition or supply missing cropped exercise text, and the original page images are not republished here.

- E. N. Lorenz (1963), *Deterministic Nonperiodic Flow*, Journal of the Atmospheric Sciences 20, 130–141. [Publisher DOI](https://doi.org/10.1175/1520-0469(1963)020%3C0130:DNF%3E2.0.CO;2). Primary source for the classic deterministic convection model; [accessible original-paper copy](https://samizdat.co/works/do-while/lorenz-1963.pdf).
- L. Illing, R. F. Fordyce, A. M. Saunders and R. Ormond (2012), *Experiments with a Malkus-Lorenz water wheel: Chaos and Synchronization*. [Author manuscript](https://arxiv.org/abs/1202.5508), [journal DOI](https://doi.org/10.1119/1.3680533). The authors describe a real apparatus and test ideal-model assumptions. Their measurements are external empirical evidence; our computations do not reproduce their apparatus or constitute new waterwheel observations.

Source pages were consulted on 9 October 2026. Analytic solutions, model reductions and numerical calculations displayed here were independently implemented for this study. No claim of a new law is made.
