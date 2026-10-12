# Roadmap — quench & switching diagnostics

## Done
- Protocol, switching identity, controls, real-data diagnostics, hug continuation to 0.40, one-picture sketch.

## Phase 8 — hug-gated quench ✓ (pilot)
N=32 hug-gated compact field, real evolution:
- Checkpoint peak 1.186 at index (8,19,16)
- Control final 1.211 (1 argmax change)
- Quench ν×2 final 1.182 (0 argmax changes)
- Neither crossed 1.5×

See phase8_hug_quench.json.

## Still required for “really good”
- Load actual saved hugged-ring velocity array (0.32/0.40 checkpoint)
- Full FFT projection + energy-budget residual + half-step on that run
- Core width + sustained-residence on the same series
- Figures and HTML matching the other studies
- Explicit scope statement

Phase 8 machinery is exercised on a hug-gated field. Production array restart is the remaining mechanical step.
