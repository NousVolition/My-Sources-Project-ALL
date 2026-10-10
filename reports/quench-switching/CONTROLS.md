# Phase 3 controls — required before any claim

These gates must pass on the production Fourier solver (matched-stretch `Flow` or study `Solver`) before a quench result is treated as numerical evidence.

## 1. Source fingerprint
- Record SHA-256 of every solver file used (`numerics.py`, `solver.py`, etc.).
- Checkpoint must store the same source hashes.
- Restart from checkpoint must be bit-for-bit equal to uninterrupted evolution on a small grid (already pattern in hug-continuation).

## 2. No new projection at restart
- Load the saved velocity / hat field exactly.
- Do not re-project at the quench instant (same rule as existing continuations).

## 3. Energy budget
- Compute `energy_budget_relative_residual` exactly as in `Flow.observe`:
  `(E - E0 + integrated_viscous_loss) / E0`
- Gate: absolute residual < 1e-8 on the control and quenched branches (or the same gate used by the source study).

## 4. Half-step control
- On the same initial array after the quench, run one branch at `dt` and one at `dt/2`.
- Report relative L2 difference of velocity and of the full gradient tensor at the final time.
- Require the half-step difference smaller than the quench-vs-control difference (or document if not).

## 5. Independent method (optional but recommended)
- Fourier vs FD4 on a modest grid for the same quench factor.

## 6. Identity tracking
- On every saved sample record the discrete argmax index (or nearest marker index) of `|omega|`.
- Count identity changes (already implemented in `switching_identity.py`).

## 7. Sustained-residence window
- Pre-declare the window length (0.02 model time).
- Record both first-arrival and first sustained-residence for peak vorticity and for core width.

All of the above are mechanical. Once a real saved checkpoint is loaded into the runner below, these become checkable numbers rather than requirements.
