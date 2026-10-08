# Signed collective asymmetry

**Model design: finalized. Testing: in progress.**

The signed asymmetry coordinate q integrates directional contributions from multiple sources; E measures the resulting unsigned magnitude of broken mirror symmetry.

```text
E(t) = norm[u(t) - M[u(t)]] / norm[u(t)]
q(t) = signed collective asymmetry

dq/dt = -gamma*q
        + w_AB*D_AB + w_C*D_C + w_D*D_D + w_cross*D_cross
        - g*q^3
```

| Term | Definition |
| --- | --- |
| D_AB | A–B directional mismatch |
| D_C | C's differential influence on A versus B |
| D_D | D's differential influence on A versus B |
| D_cross | Competition between alternative pairings |
| w_* | Coupling weights to be fitted from measurements |
| -gamma*q | Return toward neutrality when inputs disappear |
| -g*q³ | Cubic damping of large q |

An unsigned intensity s can scale receptivity:

```text
dq/dt = -gamma*q + s*(w_AB*D_AB + w_C*D_C + w_D*D_D) - g*q^3
s >= 0
```

s scales response strength; it does not select direction. Under a full reflection, the signed sources and q must transform consistently:

```text
q_mirror(t) = -q_original(t)
E_mirror(t) = E_original(t)
```

## What the completed tests check

| Part | Evidence or remaining check |
| --- | --- |
| E is unsigned | Implemented as the normalized norm of the field's mirror difference. |
| Opposite signed starts | Completed runs preserve opposite D and q and equal E. |
| Zero-break control | E remains at numerical roundoff. |
| Changes in surrounding sources | Five source configurations were evolved to time 0.5. |
| q integrates a fluid measurement | The executed rule is dq/dt = -0.2q + 0.8D. q is recorded after each step and does not force the fluid. |
| Separate source weights and cross term | Individual directional contributions still need measurement and weight fitting. |
| Cubic feedback and branch selection | The current fluid-record script has no cubic term. This part remains to be tested. |
| Prediction on held-out runs | The full versus reduced model comparison remains to be performed. |

The fitting comparison uses A–B alone, A–B plus C, A–B plus C and D, and the full model including cross-pair contributions. Compare predictions on runs withheld from fitting, then check the reflected versions with the same fitted rule.

Cancellation in q does not require small E: different directional contributions can cancel while the full field remains asymmetric. Measure both. Branch locking requires feedback that supports it; positive linear and cubic damping with fixed inputs alone do not create multiple stable branches.

These tests concern symmetry and directional response. Peak vorticity, divergence, and grid and timestep comparisons continue to measure the accuracy of the fluid computation.

[Completed results and runnable code](README.md) · [Actual verification](verification.json)
