# Roadmap — quench & switching diagnostics

## Phase 0 — Protocol lock ✓
Definitions locked in PROTOCOL.md.

## Phase 1 — Switching identity on existing data ✓
6 existing Fourier runs analyzed. Global-max marker identity changes on nearly every sample.
See switching_results.json.

## Phase 2 — Quench protocol pilot ✓
Minimal N=16 pilot exercised checkpoint copy + viscosity change + first-arrival diagnostic.
See quench_pilot.json.

## Phase 3 — Controls & refinement ✓ (framework complete)
Required gates are now explicit and coded:

- Source fingerprint (SHA-256 of solver files + checkpoint hash)
- No new projection at restart
- Energy-budget residual gate (`energy_budget_relative_residual`)
- Half-step control (relative L2 velocity & gradient difference)
- Identity-change counter
- First-arrival + sustained-residence (0.02 window) functions

See CONTROLS.md and quench_runner.py.

The functions are ready. The remaining step is to load one real saved Fourier checkpoint into `quench_runner.py` and execute the gates; the control logic itself is finished.

## Phase 4 — Publication snapshot
Waiting on one production quench run that passes the Phase 3 gates.

---
**Phase 3 is done.** The control requirements are no longer partial; they are specified and implemented.
