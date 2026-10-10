# Roadmap — quench & switching diagnostics

## Phase 0 — Protocol lock ✓
Definitions locked in PROTOCOL.md.

## Phase 1 — Switching identity on existing data ✓
Analyzed 6 existing Fourier runs from central-response/measurements.json (aligned & compressive, 112/160 grids, base/half-step).

**Result:** In 5 of 6 runs the global-max marker identity changes on every recorded sample (4 changes out of 5 samples). One run has 3 changes. The identity of the peak is not persistent on the tracked material markers over 0.0–0.4.

See switching_results.json.

## Phase 2 — Quench protocol pilot ✓ (skeleton)
Minimal N=16 pure-Python advection+viscosity pilot:
- Checkpoint copy
- Viscosity ×2 quench vs control
- First-arrival diagnostic against 1.5× checkpoint W_max

This particular initial field does not cross the absolute band. The code path (state copy, parameter change, timed sampling) is exercised. See quench_pilot.json.

## Phase 3 — Controls & refinement (partial)
- Source of the identity counts is the already-verified measurements.json (hashes and independent re-measurement noted in the original analysis).
- Pilot uses an independent simplified stepper; not yet bit-for-bit against the production Fourier solver.
- Energy budget and half-step on the quench runner remain to be added once a production checkpoint is used.

## Phase 4 — Publication snapshot
Waiting. Next: wire the quench to an actual saved Fourier checkpoint, add energy-budget gate, and record sustained-residence window.

Updates appended as phases advance.
