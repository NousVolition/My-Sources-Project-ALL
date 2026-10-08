# Done & next

[← Project home](README.md)

## Already completed

| Work | Where it lives | What is established |
| --- | --- | --- |
| Smooth initial field | [Derivation](math/notes/smooth-initial-field.md) · [code](math/initial_field.py) | Smoothness through the axis, incompressibility, decay, and exact energy formula |
| Initial-field checks | [Tests](math/tests/test_initial_field.py) · [energy check](math/verify_energy.py) | Seven original checks plus independent quadrature already completed |
| Stokes geometry | [Legacy field](math/stokes.py) · [checks](math/tests/test_stokes.py) | Reuses the submitted checks; the old Gaussian remains a comparison with an axis cusp |
| Cube & projection | [Implementation](math/navier.py) · [checks](math/tests/test_projection.py) | Compatible centered operators, residual stopping, and optional Fourier solve |
| Small & large box runs | [Results and interpretation](math/EXPERIMENTS.md) | Lengths 6 and 18, a fixed-spacing comparison, and a half-time-step check already run |
| Mirror and boundary checks | [Norm and Entropy](math/notes/norm-entropy-mirror.md) | The 65 generated pairs match after alignment; K-means depends on representation. The endpoint-matched fluid copy still has a corner and nonzero divergence |
| Soft-envelope prototypes | [Hug, imprint, and pressure](math/notes/soft-envelope.md) | Timed mirrored arms and a separate chosen pressure-opening rule pass nine checks; material-response labels and pressure for conversation data remain undefined |
| Hug joined to the ring | [Construction and results](math/notes/hug-boundary.md) | The actual closed hug now shapes the smooth ring's stream function. Smooth, divergence-free initial data with zero face mismatch; nine new checks and four short unforced box runs completed |
| Hug evolution refinement | [Saved grid and time comparisons](math/notes/hug-boundary.md#finer-grid-evolution-check) | Six saved runs through 256³ at time 0.08, with 128³ and 256³ half-step controls. Finest grid differences: velocity 0.23%, gradient 1.75%. The 256³ time-step differences are 0.0140% and 0.0582%. All 146 repository tests passed |
| Longer hug evolution | [Continuation to time 0.16](math/notes/hug-boundary.md#longer-run-to-model-time-016) | Three exact checkpoint restarts completed; two restart checks passed. 256³ time-step differences: 0.0240% velocity, 0.0819% gradient. External force 0 |
| 384³ grid comparison | [Same-time refinement](math/notes/hug-boundary.md#finer-grid-at-time-016) | At time 0.16: 0.0619% velocity difference and 0.3509% gradient difference vs 256³, using dt 0.001. Four checkpoints saved; six targeted comparison/restart checks passed |
| 384³ time-step control | [Half-step comparison](math/notes/hug-boundary.md#half-time-step-on-the-384-grid) | At time 0.16: 0.0242% velocity difference and 0.0840% gradient difference; identical starting arrays, zero force, four saved checkpoints, eight targeted checks passed |
| Hug continuation to 0.24 | [New interval](math/notes/hug-boundary.md#continuation-to-model-time-024) | Three exact restarts completed; grid and time-step comparisons retained; eight targeted checks passed |
| Hug shape measurement | [Spread of motion](math/notes/hug-boundary.md#how-the-motion-changes-shape) | Nineteen observations from saved arrays; all three runs narrow in both directions; three targeted checks passed |
| Hug continuation to 0.32 | [Evolution and shape](math/notes/hug-boundary.md#continuation-to-model-time-032) | Three exact restarts; six new checkpoints; prior shape observations reused; fourteen targeted checks passed |
| Hug continuation to 0.4 | [Evolution and shape](math/notes/hug-boundary.md#continuation-to-model-time-04) | Three exact restarts; six new checkpoints; prior shape observations reused; new data validated using unchanged code |
| Saved-flow grid-pattern check | [Measurements and limits](math/notes/hug-boundary.md#grid-pattern-check-through-model-time-04) | 31 observations from 30 distinct saved fields; checkpoint hashes and energies agree; 15 focused checks passed; no added evolution |
| Viscous energy accounting | [Measured loss and estimates](math/notes/hug-boundary.md#viscous-energy-check-through-model-time-04) | 31 measurements reused; 15 focused checks passed; remaining gap shrinks with half-size time step; no added evolution |
| One-step energy trace | [Contributions and limits](math/notes/hug-boundary.md#one-step-energy-trace-through-model-time-04) | 31 disposable one-step probes; 14 focused checks passed; exact step accounting closes to roundoff; saved trajectory remains at 0.40 |
| Local time-method comparison | [Measurements and limits](math/notes/hug-boundary.md#local-time-method-comparison-at-model-time-04) | 18 disposable controls from time 0.40; Euler and Heun use identical inputs and spatial operators; 17 focused checks passed |
| Longer time-method comparison | [Measurements and limits](math/notes/hug-boundary.md#longer-time-method-comparison-to-model-time-041) | Four branches from identical 384³ state at 0.40 to 0.41; eight reusable checkpoints; eight new focused checks passed |
| Separate trigonometric illustration | [Derivation, runs, and checks](math/notes/hidden-flow-evolution.md) | Assistant-chosen fields, distinct from the ring construction. One has an analytic smooth reduction; a 3D variant has refined finite-time results. The conversation-to-fluid mapping was not established |
| Language reports | [Report guide](language/README.md) | Existing tables and charts preserved; available chart links repaired |
| Conversation exploration | [Wordless map and observations](language/exploration/README.md) | Sequence, overlapping saved categories, reply length, and unmarked replies exposed without exporting reply text |
| Conversation study | [Full study and sources](language/conversation-study/README.md) | Three preserved snapshots, 53 selected observations, eight word paths, six context comparisons, and nine integrity regression checks |
| Optional K-means comparison | [Method and saved results](language/clustering/README.md) | Count and length-adjusted representations compared; seed sensitivity and limitations recorded |
| Two working papers | [Paper guide](papers/README.md) | Existing documents preserved, with missing source material identified |

**Verification:** all **146 repository tests passed** in the preceding publication check; the two checkpoint-continuation checks also passed. The 384³ addition passed six targeted checks, including the new 3:2 grid-transfer check. This includes the field, cube, projection, array-backend, comparison, conversation-study, and earlier prototype checks. The long numerical runs were reused when their source fingerprint matched. The conversation study separately retains its source hashes, 66 exact fragments, and 6,811 validated local HTML links. These checks establish implementation and evidence integrity; they do not establish every interpretation or a global Navier–Stokes result. The standalone LaTeX source is preserved; its earlier editor compilation remained unverified because the compiler could not find its platform directories.

The 384³ time-step control passed eight targeted checks, including rejection of a changed initial field before evolution and exact reuse of a completed control.

## Next conversation step

Review selected context comparisons with a second reader and record disagreements. For new pastes, match overlap first and extend existing evidence IDs when appropriate. [The study method](language/conversation-study/METHOD.md) keeps the new snapshots separate from the earlier 68-reply map.

Continue collecting and observing. Keep each conversation identifiable, preserve reply order, and record changes to category definitions. [The collection guide](language/exploration/COLLECTING.md) explains how to add observations without overwriting earlier work. A fixed hypothesis or fixed number of categories is not required to explore.

## Next math step

The [separate trigonometric illustration](math/notes/hidden-flow-evolution.md) retains seven saved runs and ten checks, now included in the passing 146-test suite. Its fields were chosen by the assistant and were not derived from the ring or the contributor's A/B/C meaning. The current ring study continues below; those earlier long runs were not repeated.

The user resumed the ring work by explicitly requesting the existing hug for the boundary. [That connection is now implemented](math/notes/hug-boundary.md): a C∞ gate from the closed hug is applied to the stream function, preserving incompressibility and leaving an unchanged inner core. Opposite physical faces match exactly, and external force remains zero. Nine new checks and seven existing integration checks passed. Four new runs used the existing cube solver, with no rerun of the unrelated trigonometric illustration.

The refinement now includes [six saved runs](math/notes/hug-boundary.md#finer-grid-evolution-check): 32³, 64³, 128³, and 256³ grids at time step 0.001, plus 128³ and 256³ controls at time step 0.0005, all reaching time 0.08. Completed runs were reused. Successive final velocity differences are 3.18%, 0.86%, and 0.23%; full-gradient differences are 17.03%, 6.15%, and 1.75%. On the 256³ grid, halving the time step changes velocity by **0.0140%** and gradient by **0.0582%**. The two time-step runs start from identical arrays. External force remains zero. These are finite-time numerical comparisons, without a rigorous error bound.

**Longer interval completed:** Three existing runs now continue to time **0.16**, reusing their saved 0.08 states. On 256³ with dt 0.0005, the largest sampled gradient changes from **10.500958** to **9.600445**. At the new final time, the grid comparison gives **0.3287%** velocity difference and **1.8051%** gradient difference; the 256³ time-step comparison gives **0.0240%** and **0.0819%**. External force remains zero. [Read the continued study](math/notes/hug-boundary.md#longer-run-to-model-time-016).

The [384³ comparison](math/notes/hug-boundary.md#finer-grid-at-time-016) is complete at time **0.16**, with time step **0.001** and zero external force. Relative to the saved 256³ run, the final velocity difference is **0.0619%** and the full-gradient difference is **0.3509%**. The earlier 128³/256³ values were 0.3287% and 1.8051%. The refinement ratios differ (2 and 1.5), so this alone does not establish a convergence order.

The [384³ half-step control](math/notes/hug-boundary.md#half-time-step-on-the-384-grid) is complete at time **0.16**. Starting from exactly identical arrays, time steps **0.001** and **0.0005** give **0.0242%** velocity difference and **0.0840%** full-gradient difference. The hug construction, equation, viscosity, and zero external force are retained. The previous full-step run was reused.

The [continuation to model time 0.24](math/notes/hug-boundary.md#continuation-to-model-time-024) is complete. The 256³/384³ grid differences are 0.0749% in velocity and 0.3639% in its gradient. The 384³ time-step differences are 0.0328% and 0.0977%. All three exact 0.16 checkpoints were reused; six new checkpoints preserve progress at 0.20 and 0.24. Eight targeted checks passed. The equation, field, viscosity, and zero external force are retained. These finite-time measurements supply no rigorous error bound or all-time smoothness proof.

**Shape measurement completed:** [all three saved flows narrow in both directions](math/notes/hug-boundary.md#how-the-motion-changes-shape). The finest control changes by −4.81% sideways and −6.75% vertically. No solver steps were added. A widening phase has not been observed at these saved times. That continuation is now [completed through 0.32](math/notes/hug-boundary.md#continuation-to-model-time-032).

**Continuation to 0.32 completed:** The largest sampled velocity gradient decreased at every recorded new time in all three runs. Grid differences are 0.0829% in velocity and 0.3697% in its gradient. Time-step differences are 0.0410% and 0.1099%. Both energy-weighted widths decreased at the new saved times in all three runs. [All measurements](math/notes/hug-boundary.md#continuation-to-model-time-032).

**Continuation to 0.4 completed:** The largest sampled velocity gradient decreased at every recorded new time in all three runs. Grid differences are 0.0875% in velocity and 0.3710% in its gradient. Time-step differences are 0.0491% and 0.1245%. Both energy-weighted widths decreased at the new saved times in all three runs. [All measurements](math/notes/hug-boundary.md#continuation-to-model-time-04).

**Grid-pattern check completed:** The largest measured energy fraction in the seven exact checkerboard modes was 4.061e-22. The broader short-wave band held at most 1.266e-08 of the energy. All three final derivative comparisons are recorded [here](math/notes/hug-boundary.md#grid-pattern-check-through-model-time-04). The equation and saved velocities were unchanged.

**Viscous energy check completed:** Viscosity accounts for most of the observed energy decrease. On the 384³ half-step run, the two time-integration estimates leave a gap of 0.234% and 0.279% of the observed energy loss. Halving the time step reduces both estimates of the gap. Sampling sensitivity is recorded; the balance is not exact. [All estimates and limits](math/notes/hug-boundary.md#viscous-energy-check-through-model-time-04).

**One-step energy trace completed:** The positive finite-step contribution, after pressure correction removes part of it, is the largest non-viscous contribution in every sampled step. Transport adds or removes a smaller amount. Every individual step balances within 5.527e-14 energy units; the interval reconstruction remains approximate. [Contributions and limits](math/notes/hug-boundary.md#one-step-energy-trace-through-model-time-04).

**Local time-method comparison completed:** Heun gives smaller energy and velocity changes when the time step is halved in all three local controls. These are short controls from time 0.40, without replacing the saved Euler trajectory. [Measurements and limits](math/notes/hug-boundary.md#local-time-method-comparison-at-model-time-04). Further tests await the contributor's go-ahead.

**Longer time-method comparison completed:** The averaging method remains less sensitive to halving the time step over the longer interval. Four branches share the same starting arrays and now reach 0.41. [Measurements and limits](math/notes/hug-boundary.md#longer-time-method-comparison-to-model-time-041). Further tests await the contributor's go-ahead.

## Material still missing

- The original pasted inputs for the two language scripts and the generator for the third report.
- The original connection-map image; its existing table is readable.
- The supporting source documents for the California paper and the notebook underlying the boundary paper.
- The video associated with the subtitle file.

The conversation reconstruction is preserved as a reconstruction; it is not substituted for a missing original input.

## Work history

The [consolidation record](archive/CONSOLIDATION.md) identifies which Copilot branches were reused, which duplicate implementations were retired, and which sessions produced no file changes. It is the place to check before restarting an old task.

**Original grid-comparison series complete through 0.40; separate Euler/Heun controls complete through 0.41.** Further evolution is not scheduled.
