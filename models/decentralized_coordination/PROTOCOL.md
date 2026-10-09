# Optional volunteer protocol — not yet conducted

This is a draft for an ethically reviewed, low-pressure coordination exercise. The current results are entirely simulated. The purpose is to distinguish first movement, influence, refusal and persistence of authority; it is not to assess personality, diagnose participants or settle a theory of human nature.

## Before recruitment

Use consenting adults who can freely decline without effects on employment, grades, services or relationships. Avoid recruiting people over whom the organizer has evaluative authority. Explain the task, approximate duration, data collected, possible mild frustration/social pressure, the right to skip any action, and the right to stop without penalty. Payment, if any, must not depend on unanimity or obedience. Do not label dissenters publicly or ask a leader to enforce compliance.

For organized research or publication, ask the responsible institution/research ethics body to determine review and consent requirements before collecting data. Do not self-declare an exemption from this document. Consent and voluntariness guidance is available from [HHS OHRP](https://www.hhs.gov/ohrp/regulations-and-policy/guidance/faq/informed-consent/index.html); requirements depend on jurisdiction and institution. Its linked FAQ includes older regulatory citations, so use the institution's current process.

Choose a data custodian, private storage, retention/deletion schedule, withdrawal procedure and contact before consent. Collect anonymous study IDs only; keep consent/contact records separate. Record hands/cards or digital state logs, not faces or voices, when feasible. Obtain separate optional recording consent. Do not upload identifiable participant data to this repository. Tell volunteers the latest point at which their records can be removed, and whether anonymized aggregate results will be shared. A group cannot promise that other members will forget what they saw.

Prepare large high-contrast cards marked **R** and **B** in addition to color, or accessible digital equivalents. Accommodate color vision, motor and visual needs. Provide a clearly accessible private stop/withdrawal control. Stopping or withdrawing is allowed in **all** conditions, including B: a task rule does not remove a person's right to refuse. Log deviations without sanction.

## Setup and operational definitions

Seat nine volunteers in a closed ring, numbering positions 0–8 independently of anonymous identities P0–P8. Arrange screened card holders or a local digital display so B/C participants see **only their two immediate neighbors' cards**. Merely seating people in a circle with every card visible would fail the local-information condition. The facilitator can monitor the full state privately to stop a round, but must not give hints about global color counts or progress. Use the same neutral stop signal in all conditions.

For A, add a visible central command display controlled by the designated leader; followers need only that command. This information difference is part of A and must be reported. No negotiations or discussion in any condition; A's command is permitted. Do not allow additional gestures or signaling beyond changing cards. If this is violated, record it; do not humiliate or exclude a participant retrospectively to improve the result.

Randomize colors independently with fair probability .5 at every new matched block. Use the same colors by **seat**, same seating permutation and same potential A leader for A/B/C within the block. Keep unanimous starts and score them at zero seconds. Randomize which condition occurs first. Between matched blocks, change seating by a saved random permutation. Mark practice separately and never analyze it as a pilot. Use a reproducibly generated schedule, such as `results/pilot-schedule.csv` for a demonstration, or a new recorded seed for a fresh pilot. The saved file contains no human observations.

Unanimity means all nine cards have the same R/B state, verified by the facilitator. End immediately on first unanimity or at 120 seconds. Record timestamps to at least one second; for simultaneous changes preserve ties instead of inventing an order. Mark timeout explicitly; consensus exactly at 120 counts as success. If a participant withdraws, stop the round and mark withdrawal/missing, not nine-person nonconsensus. Do not replace the participant mid-round. Analyze withdrawals separately from timed-out complete rounds.

## Instructions read before every condition

Common script: “The aim of this round is for all nine cards to show one color. Either color counts. You can stop participating or decline any action without penalty. We will stop on agreement or after two minutes. Please do not talk or signal beyond the allowed card actions. This is a task about coordination, not a judgment of you.”

**A — fixed hierarchy.** Randomly designate a leader before revealing the starting colors. “The leader will command the color on their own starting card and keep that command unchanged this round. On each ten-second signal, followers are asked to display the commanded color. The leader holds the commanded color. You still have the right to decline or stop.” The leader may display the command but may not threaten, bargain, criticize or change strategy. Record any declined directive and deviations; do not assume perfect compliance from the simulation. Leader identity is fixed for a round and counterbalanced across later blocks.

**B — deterministic local rule.** “At each ten-second signal, compare your card with your two neighbors. If both neighbors have the opposite color, switch; otherwise hold.” Implement a short decision/reveal procedure at each tick: participants decide from the still-visible old cards, commit their choice face down, then reveal together. Do not let later responders react to already-updated cards. Log actual reveal latencies. This is a practical approximation to the simulator's simultaneous update; deviations and missed ticks are outcomes for feasibility review.

**C — voluntary local coordination.** “You can change your card whenever you choose, keep it as it is, decline to match a neighbor, or start a change. There is no designated leader and no required rule for changing. Please do not talk.” Use no ten-second action signals during the primary C round. Record changes continuously. Do **not** instruct people to follow the simulator's .25/.02/.15 probabilities. Those are hypotheses about a toy mechanism, not the human task. A separately labeled synchronized-C extension may be added only if planned in advance.

Do not use descriptions such as “obedient people,” “water people,” or “free people” during instructions. They encourage demand effects. In debriefing explain that B is only a water-inspired local-interaction analogy and does not reproduce water physics.

## Initial pilot: four rounds per condition

Run four matched blocks, each containing A/B/C in a newly randomized order: 12 total task rounds, at most 24 minutes of task time, plus orientation, practice, short private reports and breaks (allow roughly 45–60 minutes). Take a short break between blocks and offer additional breaks. The goal is to check visibility, comprehension, timing, comfort and record quality; four rounds per condition cannot establish a population effect.

After each round, separately and privately ask whether the participant deliberately initiated a switch, deliberately declined a perceived request to change, strategically waited, did not see a clear cue, or is unsure. Ask whom they attended to and whether they felt free to refuse. Let them skip any question. Do not display answers to peers. Record self-report separately from observed behavior, with an “uncertain” value. Asking about refusal may itself affect later behavior; keep wording and timing identical across conditions and report that limitation.

Pilot success criteria are understandable instructions, local-only visibility where intended, reliable timing and no unaddressed pressure or discomfort. Decide any rule or logging changes before a new confirmatory sample. Preserve the pilot's actual protocol version; never rewrite its record to match a later design.

## Stronger optional human study

Use multiple independent nine-person groups, with group as the sampling unit. Choose sample size through prospective precision or power calculations after defining a practically relevant group-level effect and accounting for repeated measures; 4,000 simulated blocks are not 4,000 human groups. A single group of nine cannot separate personality, existing relationships and network effects reliably.

Preregister outcomes, exclusions, stopping rule, hypothesis direction, missing-data treatment, condition timing and planned sensitivity analyses. Counterbalance all six condition orders across groups and blocks; test or report condition-by-order and previous-condition effects. Matched repeats can be recognized by people, so include block index and learning/carryover checks. If recognition dominates, consider a separately planned between-group matched design instead of treating repeated human rounds as memoryless.

For identity-versus-seat inference, plan at least a full nine-mapping crossing per group where feasible, split into sessions with breaks if needed. Random cyclic shifts of a randomized base mapping put every identity in every position exactly once, but rotate neighbor identities together; include additional independently shuffled mappings if separating stable neighbor relationships is important. Balance leader assignments independently of seat and condition order. The simulation's nine cyclic shifts likewise retain neighbor identity pairs on a ring, a limitation for claims about dyadic relationships.

The ring has no intrinsically privileged seat. To test topology rather than arbitrary seat labels, add a **separate**, consented, preregistered C extension with unequal network access, such as a star implemented through local displays. Rotate every identity through hub and leaf roles while holding starts by position. Do not call a leaf's one-neighbor rule condition B. Star versus ring changes both degree and information; do not claim to isolate abstract centrality.

To test persistence after hierarchy, cross prior A exposure with C and a seat permutation, balancing whether the former leader keeps their seat, moves to a new seat, or (in the separate star extension) moves into/out of the hub. Reset colors. Use counterbalanced comparisons with C preceding A to distinguish earlier status from preexisting identity influence. Collect private attention reports and perceived leadership nominations; account for local access and opportunities to influence. Avoid public rankings.

## Recording and analysis

Use the blank `forms/round-log.csv`, `forms/event-log.csv`, and `forms/private-report.csv`. Store completed forms privately, with a protocol version and seed. Preserve observations even if they conflict with the simulation. Never merge human data into files named simulated results.

Report group-level consensus proportion, capped time through 120 seconds, observed switch counts, tied first movers, flip-backs, global state revisits, intended initiatives, self-reported refusal, perceived freedom, and protocol departures. Distinguish observation from self-report; neither holding nor copying proves intent. For persistent hierarchy, report the same identity's influence before/after permutations separately from the same seat's influence. Include ties, undefined measures and disagreements between coders.

For humans, the simulator's exact donor is unavailable. Predefine an observable lagged response measure (a neighbor changes, followed by adoption of that color within a fixed window, adjusting for other available cues), with ambiguity flags when both neighbors could be sources. Combine this with optional private attention reports; do not label lagged association causal influence. An actual randomized perturbation design would require separate planning and consent. Two blinded coders should independently score a prespecified subset of recordings and report agreement and unresolved ambiguities.

Use paired within-group comparisons with resampling at the independent **group** level; keep all that group's rounds together. Analyze timeouts as right-censored and show restricted mean time and event proportions together. For identity/position, use within-group crossed effects with local initial configuration, round, order, former leader and exposure opportunities considered; evaluate on held-out groups or full blocks. Tiny pilots do not justify complex fitted models or reliable p-values. Predefine multiplicity treatment and label unplanned analyses exploratory.

Debrief without implying that disagreement is failure as a person. Explain that consensus speed is only one outcome and is compatible with different amounts of pressure, intention and satisfaction. Report what was learned and what remains uncertain. No outcome establishes or refutes Freud's broader claims.
