# Next test: do identical waters develop nested collective patterns?

Status: a proposed follow-up, not a completed experiment. The completed experiment is documented in README.md and data/results.json.

## Target

Look for recurring collective configurations, with faster changes within a configuration and slower transitions among configurations, while every water molecule retains the same physical parameters. A collective state is a configuration of many molecules; it is not a privileged molecule or a permanent group identity.

The supplied figures and equations are from [Heteroclinic networks for brain dynamics](https://www.frontiersin.org/journals/network-physiology/articles/10.3389/fnetp.2023.1276401/full). They motivate a hypothesis, rather than provide a water mechanism. In particular, Figure 5 explicitly assigns a pacemaker to one site. The primary water test should contain no assigned molecular pacemaker.

## Water-specific measurement

1. Run longer trajectories under one fixed environment after a defined equilibration interval. Analyze field-free and uniformly driven conditions separately; do not call the imposed on/off schedule spontaneous state switching.
2. Save positions and orientations frequently enough to resolve both hydrogen-bond rearrangement and candidate slower collective changes. Predefine several temporal resolutions and verify that conclusions survive resampling.
3. Build local collective descriptors from spatial neighborhoods: dipole alignment, hydrogen-bond graph structure, density, and their correlations. Discover candidate recurring states on a training trajectory, without assigning three groups or nine states in advance.
4. On separate trajectories, measure state occupancy, dwell-time distributions, transition probabilities, within-state fast dynamics, and cross-scale prediction. Report weak, absent, and ambiguous recurrence as well as positive findings.
5. Compare independent initial states, longer runs, a larger box, another water model, and a polarizable model. Use state-label permutations and spatial/temporal surrogate data that preserve simpler statistics. Account for time correlation when estimating uncertainty.

## Testing a proposed reduced equation

The supplied activity model is

    dA_i/dt = delta * Laplacian(A_i) + sigma*A_i - gamma*A_i^2
              - sum_{j != i} rho_ij*A_i*A_j + eta*abs(xi_i(t)).

To apply it to water, define each nonnegative A_i from an observable collective activity. Estimate coefficients from training molecular trajectories and test predictions on held-out trajectories. A simulation using chosen coefficients is an illustrative activity model until that mapping is validated.

Equal local parameters do not make every network role equivalent. A specially chosen competition matrix or an assigned pacemaker can build in the desired hierarchy. Compare the proposed coupling with zero coupling, symmetry-preserving coupling, and suitable shuffled topology. Compare an assigned-pacemaker case separately with a case containing none.

Specify the noise process before integration. The absolute value of ideal white noise is not an ordinary function. Naively adding the absolute value of a Brownian increment at every step creates a positive drift that grows as the numerical step shrinks. One well-defined candidate control is a rectified, piecewise-constant Gaussian signal with a fixed physical update interval; step refinement must keep the same noise realization and update times. This is a stated regularization choice, not a transcription of an unspecified noise convention in the paper. Include zero noise and several noise amplitudes.

## Evidence needed for a heteroclinic claim

Recurring plots alone are insufficient. Identify candidate invariant saddle states or sets in the reduced dynamics, their local stable and unstable directions, and connecting transition paths. Distinguish an engineered network, a fitted approximation, a noise-driven sequence of metastable configurations, and a deterministic heteroclinic network. Slow and fast time scales alone do not establish hierarchical chunking.

## Decision rule

Accept a proposed reduced model as useful only if it predicts held-out collective measurements more accurately than simpler baselines and its patterns survive the physical and numerical controls. Do not tune it solely until its plot resembles the supplied figure. A failed or unresolved test is an informative result.

## Sequence learning: an additional, unrun experiment

The later excerpt adds two distinct questions: can a fitted network learn to imitate water trajectories, and can an engineered adaptive interaction model learn a sequence? The first is inference about the data. The second changes the dynamical rules; it is not a property already implemented in fixed-charge TIP3P water.

For an adaptive network test, specify the learning rule and allowed trainable parameters before training. Present a sequence, freeze learned parameters, remove all teacher forcing, and test continuation from incomplete or noisy cues. Keep training and test sequences separate. Compare with untrained weights, scrambled training order, no coupling, and teacher forcing left on as a positive tracking control. Measure correct-order transitions, recall duration, failure rate, and sensitivity to cue/noise across independent seeds. Teacher-on imitation alone is not free recall.

For relaxation after a quench, report first-passage brackets or a justified decay fit and its uncertainty. Distinguish a long transient from persistent storage. Repeat multiple quench depths, observation durations and noise levels before attributing a delay to a particular basin or to hierarchical dynamics. The completed water experiment reports only the first half-decay brackets after field removal; its on/off protocol does not establish a heteroclinic bifurcation.
