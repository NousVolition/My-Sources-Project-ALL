# Fresh independent confirmation: declared before new research simulations

Status: prepared; the existing frozen full-history refinement remains active. **No fresh pilot, train, validation or test research trajectory has been generated.** Small analytic and synthetic unit-test fixtures are explicitly separate from research data. This document and protocol.json define the next phase; neither contains a finding of benefit.

## Question and observations

Does the past reorganization of neighboring passive markers improve prediction of future strong deformation at their current material locations, above their present arrangement and the full present local velocity gradient?

Predictor A uses the parent 31 current features: three velocity entries, all nine local velocity-gradient entries, three strain eigenvalues, strain norm, rotation norm, speed, four neighbor-distance summaries, six shape-tensor entries, strain alignment and two neighbor-angle summaries. B uses exactly A plus the parent23 historical features: distance changes/path/rate/curvature, angle changes/path, neighbor turnover, alignment/shape changes and fitted finite-neighborhood past deformation. Exact names/formulas remain in [features.py](../features.py). Neither has particle ID, future velocity, future positions or a future-derived neighborhood.

The default past window is .2, ending at each anchor .3/.4/.5/.6/.7. A .025 gap separates inputs from the forward .2 deformation interval. The continuous target is forward logarithmic maximal singular stretch divided by horizon, derived from tangent-map composition over that strictly future interval. Event labels use the top10% **training** target quantile, which is then frozen. This predicts finite-time deformation at sampled material trajectories; it does not locate the exact Eulerian maximum or the first onset time.

## Independent runs and sealed testing

The fixed split is32 training initial conditions (51000–51031),16 validation (52000–52015) and64 untouched test conditions (53000–53063). Four separate numerical pilot fields (50000–50003) are excluded from every fit and confirmatory interval. All markers and all times from an initial condition stay in its assigned split. Initial markers use seed+80000 as in the parent. There is no optional expansion of test seeds after results are observed.

Train-only StandardScaler and ridge/logistic fits use the same fixed grids of regularization choices for A/B. Validation selects parameters and an event decision threshold. Training and validation are not merged for a later fit. All primary, control and observation-sensitivity models, scalers and thresholds are hashed and published before **any** test trajectory is generated. The code refuses fitting if a test archive, claim or checkpoint exists; test generation requires a verified model-lock publication. Test outcomes cannot select a model, resolution, threshold, endpoint or further seed batch.

## Primary endpoint and practical size

At each anchor, compute average precision for ranking markers that will have a high-deformation event. Average the five anchor scores within a run, then average equally across64 runs. This compares locations at the same time, avoiding an apparent gain obtained only by ranking different times or runs. A no-positive anchor has AP0 for both models and is counted. Pooled AP, precision and recall at a fixed10% marker budget, calibration and continuous metrics are secondary.

The primary effect is **spatial AP(B) minus spatial AP(A)** for the linear logistic models. A useful improvement is predeclared as **.005 absolute AP (0.5 percentage points)**. A paired whole-run percentile bootstrap uses10,000 draws, seed84017, two-sided95% intervals. Useful benefit is supported only if the lower bound exceeds .005 and all required numerical audits pass. An improvement of at least .005 is excluded for this fixed comparison if the upper bound is below .005 and the audits pass. Otherwise evidence remains unresolved, or indicates a smaller positive effect. Any numerical failure prevents either usefulness or exclusion claims.

The secondary continuous effect is mean per-run relative RMSE reduction, with useful size .01 (1%). It cannot substitute for the event-location conclusion. MAE, R², Brier score, log loss, calibration bins and ECE10 are reported. Secondary/control intervals are exploratory, with no claim of multiplicity-adjusted discoveries.

## Numerical prerequisite and preserved fluid model

Finish, independently verify and publish the existing 56/80/112/conditional160/224 frozen experiment first. Its result may remain numerically unresolved at its hard bound; that limitation is preserved. The new phase then has its own accuracy pilot rather than claiming that the old trilinear predictions have converged.

The Eulerian solver, same smooth random low-mode initial-condition family, projection/dealiasing, joint Heun method, periodic side6, primary viscosity.02, zero mean and zero force are unchanged. The new observation method is periodic four-point tensor Lagrange interpolation of velocity with the exact piecewise Jacobian of that interpolant. This is a documented tracer-accuracy change, motivated by the observed current-gradient drift. It does not alter the Eulerian dynamics, add a force, introduce material memory or substitute a toy system. The interpolant is not exactly divergence-free and its derivatives can jump at cell faces.

Four pilot fields run at56/80/112 with two predefined half-step controls at80 and112. If necessary, four160 fields and two160 half-step controls follow. The smallest candidate80/112 must pass both adjacent-grid comparisons, contracting target differences and its half-step checks. Pilot comparisons require target relative RMS≤.25%, event Jaccard≥.98, gradient relative RMS≤.5%, velocity relative RMS≤.1%, standardized present/history drift≤.5%, old frozen prediction drift≤1%, accumulated probe-volume error≤.5%, and temporal target error≤10% of spatial error and≤.025% relative. Original train feature scales and original training event cutoff are used only for these pilot measurement screens.

After fresh train/validation fitting, the frozen fresh predictors must also pass pilot probability drift≤1%, spatial AP contrast drift≤.000625 macro and≤.00125 per pilot run, and RMSE contrast drift≤.25% of validation A RMSE. A numerical failure at80 may advance to112 and refit using new train/validation trajectories **before any test data exist**. If112 fails, stop with a published numerical limitation. No statistical benefit can justify a grid change.

After model freezing, test seeds53000/53001 also receive a predefined finer-grid and half-step control, with the same gates. A failed test audit withholds the usefulness/exclusion claim; it does not trigger refitting or retrospective resolution selection. Primary histories run to1.2, save every.025, default step.0025 and512 markers. The absolute maximum is194 new trajectories, including every pilot, conditional train/validation replacement and control. No longer intervals or added seed batches are authorized automatically.

## Controls and prechosen sensitivities

All history controls are refitted independently with the same fitting budget: history rows shuffled within run/anchor using seeds91/92/93, recorded histories reversed, farthest-current-marker donor histories, and mirrored/mirrored-reversed recorded histories. The mirrored histories keep current fluid inputs fixed and are diagnostic representations, not physical reflection interventions. A retrained reversal can preserve information, so its success cannot establish chronology's irrelevance.

There are shuffled training-label controls, consistent complete marker relabeling tests with no identity input, a31+23 current quadratic-feature control (seed513), past-deformation-only/no-past/no-exchange ablations and a separate impoverished baseline. Secondary tree A/B/matched-current models have identical depth, leaves, iterations and tuning budget. This distinguishes improvement above full current gradient from improvement above an impoverished observation or extra fitting flexibility.

Ten declared observation sensitivities reuse raw trajectories: lookback.1/.3, horizon.1/.4, markers256, neighbors4/16, every-other past sample, and event quantiles.8/.95. Their A/B fits and cutoffs are frozen before tests. Viscosity.01/.04 variants reuse four fixed test initial conditions with frozen primary models; they are dependent variants, not eight additional independent test groups.

## Interpretation and reproducibility

Publish protocol, code/tests, raw archives with configurations and source hashes, fitted models, lock/publication records, per-run predictions/metrics, uncertainty, numerical failures and plots. A preselected spatial example uses test seed53000 at anchor.5. Preserve failed checkpoints and do not duplicate active numerical jobs. Report accurate execution status and environment. The original sharp-field periodic-join/off-central-maximum limitations remain separate from this smooth-flow family.

Any supported benefit is out of sample within this particular distribution, observation method, model class and time window. History can improve predictions from an incomplete local state; deterministic evolution does not permit fundamental extra information beyond a complete present state. No correlation establishes causation, material memory, leadership, singularities or a Navier–Stokes breakthrough. A small/null comparison cannot rule out all possible histories, predictors or fluids.
