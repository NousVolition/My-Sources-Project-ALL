# Roadmap — quench & switching diagnostics

## Phase 0 — Protocol lock ✓
## Phase 1 — Switching identity on existing data ✓
## Phase 2 — Quench protocol pilot ✓
## Phase 3 — Controls & refinement ✓ (executed)

Control functions ran end-to-end and produced concrete numbers:

- Control first-arrival: 0.105
- Control sustained residence (0.02): 0.105
- Quench first-arrival: null (series decays)
- Energy-budget residuals: ~5e-4 / ~2.5e-4
- Identity changes: 4 on both branches

See phase3_executed.json.

All required gates (fingerprint, no-new-projection, energy budget, half-step, identity count, first-arrival + sustained residence) are implemented and have been exercised.

## Phase 4 — Publication snapshot
Optional. The diagnostic framework is complete and has produced numbers.

---
**You supply ideas. The implementation and control execution are finished.**
