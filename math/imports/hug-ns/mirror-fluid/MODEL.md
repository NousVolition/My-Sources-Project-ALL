# Signed collective asymmetry

**Model status: being examined and tested.**

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

## C and its mirror partner D in the existing tests

In `fluid_sources.py`, define the source D as the reflected field `M[C]`. The code already runs these combinations of the symmetric A–B field and the surrounding sources:

| Code name | Initial surrounding field | Mirror mismatch E at time 0.5 |
| --- | --- | --- |
| ab | None | 4.30e-15 |
| ab_cd | 0.5 C + 0.5 D | 4.19e-15 |
| ab_uneven | 0.7 C + 0.3 D | 0.0131763 |
| ab_c | C | 0.0329282 |
| ab_c_mirror | D | 0.0329282 |

Equal contributions preserve mirror symmetry. The unequal and single-source starts have measurable asymmetry. Reflected single-source starts give opposite signed measurements and equal E. In these tests D is constructed by reflecting C; its shape is not independently varied.

The JSON column named `D` is a signed measurement, distinct from the surrounding source named D above. Source weights in this table set initial velocity contributions; they are not fitted coupling coefficients in the q equation.

## What the completed tests check

| Part | Evidence or remaining check |
| --- | --- |
| E is unsigned | Implemented as the normalized norm of the field's mirror difference. |
| Opposite signed starts | Completed runs preserve opposite D and q and equal E. |
| Zero-break control | E remains at numerical roundoff. |
| Changes in surrounding sources | Five source configurations were evolved to time 0.5. |
| q integrates a fluid measurement | The executed rule is dq/dt = -0.2q + 0.8D. q is recorded after each step and does not force the fluid. |
| Separate source weights and cross term | [Dynamic test completed](dynamic-q/README.md): four explicit interaction features fitted from 10 training starts, with independent AB bias and independently shaped D. This is a six-template implementation; its weights are not unique causal contributions. |
| Cubic feedback and branch selection | Cubic damping was fitted and compared with the linear model on four unseen starts. This does not establish branch selection: nonnegative linear and cubic damping cannot create two isolated stable branches under constant inputs. |
| Prediction on held-out runs | Four complete starts withheld from fitting, with reflected and half-step controls, plus 49- and 65-grid repeats. Both measured-input predictions and forecasts with inputs frozen at time 0.1 are reported. |

The new source-addition test measures an interaction contrast relative to the fixed A–B background. It does not uniquely decompose the evolving fluid or establish competition between pairings. Its 33³ results pass the time-step comparison; spatial convergence is not established.

The completed fitting comparison includes AB alone, AB plus C, AB plus C and D, all source interactions, and all interactions plus cubic damping. The fluid coordinate q_fluid is measured independently; the historical passive q_record is retained as a separate control. [Definitions, prediction errors and numerical limits](dynamic-q/README.md).

Cancellation in q does not require small E: different directional contributions can cancel while the full field remains asymmetric. Measure both. Branch locking requires feedback that supports it; positive linear and cubic damping with fixed inputs alone do not create multiple stable branches.

These tests concern symmetry and directional response. Peak vorticity, divergence, and grid and timestep comparisons continue to measure the accuracy of the fluid computation.

[Completed results and runnable code](README.md) · [Actual verification](verification.json)


