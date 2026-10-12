# Roadmap — quench & switching diagnostics

## Phase 0–3 ✓
Protocol, switching identity on existing data, controls, and real-data diagnostics from central-response.

## Hug quench pilot ✓
N=32 hug-inspired compact field evolved with centered advection+viscosity.
- Checkpoint at t=0.10, max_grad ≈ 2.67
- Control (nu=0.01) final ≈ 2.57
- Quench (nu×2) final ≈ 2.38
- Neither crossed 1.5× on this interval

See hug_quench_pilot.json.

## Next
- Wire the same branching logic to a real saved hugged-ring velocity array (from hug-continued-0.4 or earlier).
- Add energy-budget residual and half-step control on that production runner.
- Record first-arrival / sustained-residence for peak vorticity and core width.

The diagnostic machinery is exercised on a real compact field. The production checkpoint wiring is the remaining mechanical step.
