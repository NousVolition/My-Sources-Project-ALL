# Findings from the saved simulation

**All results below are simulated. No human participants were studied.** Seed 20261008; Python 3.12; 45,852 stochastic/deterministic task rounds plus exhaustive enumeration of 512 B starts. The four pilot rounds per condition are separate from the main analysis.

## Coordination outcomes

| Main condition (4,000 rounds each) | Consensus by 120 s | 95% Monte Carlo interval | Mean capped seconds | Mean switches | Mean modeled refusals |
| --- | ---: | ---: | ---: | ---: | ---: |
| A: fixed leader, perfect compliance | 100.00% | 99.90–100.00% | 9.97 | 3.98 | 0 |
| B: deterministic two-neighbor rule | 31.00% | 29.59–32.45% | 88.95 | 3.07 | 0 |
| C: primary voluntary-choice model | 96.73% | 96.13–97.23% | 34.34 | 18.33 | 2.95 |

The mean capped time includes every timeout as 120 seconds; it is not the mean time to eventual agreement. A's near-10-second average includes already-unanimous starts at zero. A is guaranteed to succeed in this idealized primary implementation. Lower simulated compliance slows it: with .3 compliance per differing directive, 94.17% succeeded and capped time averaged 58.96 s in the separate sensitivity sample.

C exceeded B by 65.73 percentage points in the matched main sample (95% block-bootstrap interval 64.13–67.30 points). This does not isolate the effect of voluntary choice: the conditions also differ in timing, randomness, imitation and opportunities to initiate. The two-neighbor rule is unusually restrictive; it is not a representative sample of decentralized algorithms.

Exact B enumeration yields **152/512 = 29.6875%**, rather than the sampled main estimate of 31%. All 512 starting configurations settle to fixed points within 40 seconds, including 360 disagreements. No lasting cycles occur for this exact nine-seat ring. The Monte Carlo sample and the exhaustive population are distinct denominators, and their difference is expected sampling variation.

## C depends on its behavioral assumptions

In the separate 1,200-start sensitivity sample, baseline C succeeded in 97.33% of rounds. Slowing active opportunities from .25 to .08 per seat per second reduced success to 63.17%; increasing them to .6 yielded 100% within this sample. Simultaneous ten-second C, with everyone active at each tick, reached 42.58%. That changes both timing and opportunity count and is not a pure synchrony test.

With refusal .9 and initiative .1, success fell to 27.25%. Refusal 1 with no initiatives produced **zero switches** and only the seven initially unanimous starts succeeded (7/1,200). Therefore the capacity to refuse is compatible with consensus under some policies, but cannot guarantee consensus under all permitted policies. Raising initiative is not uniformly beneficial: it can also disturb configurations near agreement. These probabilities are chosen examples, not human parameter estimates.

## Apparent leadership: structure dominates in the constructed star

“Influence” here counts actual copied switches credited to their source by the simulator. It is not a measure of authority, intention or human causation. Four separate topology/trait scenarios each include 240 independent starting-state replicates and nine crossed identity/seat mappings.

| Scenario | Held-out prediction improvement: position | Identity | Both |
| --- | ---: | ---: | ---: |
| Ring, homogeneous choices | −0.10% | −0.03% | −0.13% |
| Ring, stipulated identity-0 traits | −0.03% | 0.14% | 0.12% |
| Star, homogeneous choices | 65.53% | −0.002% | 65.52% |
| Star, stipulated identity-0 traits | 64.43% | −0.04% | 64.39% |

Prediction improvement is reduction in held-out mean squared error against the training-mean baseline. Negative values indicate worse prediction. On the symmetric ring there is no intrinsically special seat. Even after inserting identity traits, most variation in copied-switch counts is unexplained by fixed seat or identity effects. On the star, the hub is visible to eight others while each leaf is visible only to the hub, so hub dominance follows directly from access and the local copying mechanism. A different mechanism could give a different result.

In the homogeneous star, the hub averaged 5.41 credited switches per round, while leaves averaged about .066–.080. The same top seat persisted across adjacent identity permutations with a tie-aware score of 99.67%; the same top identity persisted only 0.066%. That near-zero identity persistence partly reflects the design: every permutation moves a new identity into the dominant hub. It must not be compared to 1/9 as if mappings were independent. On the homogeneous ring, identity persistence was 11.24% and seat persistence 14.05%; shared starts can cause the latter difference without a permanently privileged seat.

The identity-0 advantage was deliberately inserted (attention weight 4, refusal .65, initiative .03). Its small ring effect neither establishes nor rules out real personality effects. Position versus identity is not a universal either/or: the answer depends on the network and behavioral mechanism.

## Hierarchy persistence and initiators

A's copied-switch concentration is 1 by construction whenever a follower changes. C's average concentration is .211 over defined rounds; B's is .286. Concentration across a finite round is not equivalent to a recognized hierarchy. Spontaneous initiators and first observed movers are recorded separately, with ties preserved. A first mover need not be the most-copied source, and neither need intend leadership.

In matched C sensitivity rounds occurring after A, the former leader's average share of copied switches was 10.85% under no status memory, 12.51% with attention weight 3, and 12.39% with weight 6. The last two differences are modest, not monotonic, and are not presented as statistically established effects. The weights are inserted assumptions and do not discover persistence in people. Baseline simulations contain no learned social status at all.

## Interpretation and next evidence

The simulation establishes the behavior of specific algorithms. Local voluntary imitation with nonzero refusal can reach agreement in this constructed setting. A strict local rule can freeze. A network hub can look like a leader regardless of which identity occupies it. These are useful mechanisms and testable hypotheses, not empirical conclusions about human freedom.

Optional volunteer testing should preserve the right to decline in every condition, privately distinguish deliberate refusal from waiting, record learning and pressure, and counterbalance both roles and positions across multiple independent groups. The supplied protocol and blank forms specify that separation. The results do not support claims about actual water, survival, welfare, or proving/disproving Freud.
