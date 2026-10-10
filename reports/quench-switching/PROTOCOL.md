# Quench and switching-identity diagnostics

Protocol locked before any new evolution. All definitions below are fixed for the first pilot.

## Diagnostics
- Peak vorticity: \(W_\max(t) = \|\omega(\cdot,t)\|_\infty\).
- Core width: half-local-spin chord at the current global-max location (12 transverse directions), as already defined in central-response.
- Energy: box kinetic energy.

## Bands (pre-declared)
1. Absolute: enter when \(W_\max \ge 1.5 \times W_\max(0)\).
2. Relative: enter when \(W_\max \ge 0.90 \times\) running maximum so far.

## Times
- First-arrival: earliest sample at which the diagnostic enters the band.
- Sustained residence: earliest sample at which the diagnostic has remained inside the band for a continuous window of length 0.02 (model time) without exiting.

## Switching identity
Track the discrete grid index (or nearest tracked material marker) that currently holds the global maximum of \(|\omega|\). Count the number of changes of that index over the recorded times. A change is recorded only when the argmax index differs from the previous sample.

## Quench (Phase 2)
- Start from a saved checkpoint of an existing unforced run.
- Copy the complete velocity state.
- Quenched branch: multiply viscosity by a pre-declared factor (pilot: \(\times 2\) and \(\times 0.5\)) at the checkpoint time; no other change.
- Control: continue with original viscosity.
- Optional gradual control: linearly ramp the factor over 0.01 model time.
- Evolve both under identical subsequent forcing (zero) and record the times above.

## Controls required before claims
- Independent solver or half-step on the same initial array.
- Energy-budget residual gate.
- Source fingerprint + bit-for-bit restart check.
- No new initial projection at restart.

## Scope
Finite-time numerical observation on the existing periodic-cube discretization. Does not establish continuum regularity or breakdown.
