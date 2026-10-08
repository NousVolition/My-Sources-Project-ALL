# Matched-stretch numerical controls

Snapshot: 2026-10-08T23:28:22.747514+00:00. **5 of 48 runs complete.** Further runs are active or queued.

The available 64³ runs amplify spin but fail the spatial-resolution screen. Small time-step differences alone do not resolve this issue. The 128³ and 256³ evolution comparisons are pending.

The supplied starting formula is retained: periodic cube side 6, Heun stepping, baseline viscosity 0.001, original mean velocity, and zero external force. Other viscosities are explicit in each run. Its raw strain has unequal curl at opposite box faces, and its initial maximum spin lies outside the central tube. A common smooth periodic limiting start is not established by these runs.

![Spin, integral and resolution](../files/matched-stretch-snapshot.png)

## Saved outcomes

| Run | Status | Time | Peak W at saved times | I at endpoint | Largest cutoff enstrophy % |
| --- | --- | --- | --- | --- | --- |
| baseline-n128-base | running | 0.00 | 59.04298 | 0.00000 | 0.773 |
| baseline-n64-base | complete | 0.40 | 701.78986 | 165.97088 | 67.036 |
| baseline-n64-half | complete | 0.40 | 701.77054 | 166.00092 | 67.035 |
| perturb-n64-a0.001-s101 | complete | 0.40 | 766.57401 | 167.39999 | 66.938 |
| perturb-n64-a0.001-s202 | complete | 0.40 | 768.18715 | 164.80941 | 66.793 |
| perturb-n64-a0.001-s303 | running | 0.04 | 81.44513 | 2.78197 | 3.622 |
| viscosity-n64-nu0.01 | complete | 0.40 | 582.85511 | 130.46923 | 47.426 |

## Time-step comparison

| Comparison | Through | W curve L2 difference % | I curve L2 difference % |
| --- | --- | --- | --- |
| timestep-64 | 0.40 | 0.10865 | 0.01478 |

Errors compare identical saved times; the reference is the half-step run. Full definitions and other observables are in [analysis.json](analysis.json).

![Widths and budgets](../files/matched-stretch-budget-width.png)

## Measurements

W is maximum grid-point vorticity magnitude. I integrates W at every accepted step. Energy is ½∫|u|²; enstrophy is ½∫|curl u|². The spin ratio uses mean vorticity magnitude. Width is the shortest measured transverse half-peak chord; both domain units and cell counts are recorded. Fixed thresholds are 50, 100, 200 and 400. All units are model units.

Finite-amplitude perturbation rates and recurrence diagnostics are in the analysis file. They describe the recorded window and do not establish an asymptotic Lyapunov exponent or recurrent attractor.

![Perturbations and finite-window recurrence diagnostics](../files/matched-stretch-perturbation-recurrence.png)

[Protocol](protocol.json) · [Initial audit](initial-audit.json) · [Implementation checks](implementation-check.json) · [Snapshot verification](../numerical-progress-verification.json) · [Run measurements](runs/)

Source is preserved in numerics.py and run_suite.py. To reproduce the matrix, copy those files, protocol.json, requirements.txt and reproduce.py into a fresh directory, install the requirements, and run `python reproduce.py`. Generated arrays and restart files remain local; the published JSON files are measurement records, not restart checkpoints.
