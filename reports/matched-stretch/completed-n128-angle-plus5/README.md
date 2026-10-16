# Completed grid-128 strain-axis tilt +5 degrees — October 16, 2026

[Matched-stretch reports](../README.md) · [Previous completed parameter batch](../completed-n128-parameter-batch-20261015/README.md)

The original **angle_degrees-n128-v5.0** run reached **0.40** with **41 saved fields**, all independently remeasured. The frozen supplied initial-field family, per-run viscosity, mean, domain, Heun method and zero force were preserved. No new simulations or horizon extensions were introduced.

This authorized one-factor case tilts the supplied raw strain axis by +5 degrees and retains viscosity 0.001, domain side 6, grid 128 and zero force. It preserves its original Heun step of 0.01/290 = 3.4482758620689657e-05. The exact prescribed initial field was checked using the frozen initial-field routine with angle_degrees=5.0. Its initial mean is [0.37293302996088135, 0.5203188902456848, -1.1643092014970942]; the zero-angle baseline has a different prescribed initial mean. Each trajectory conserves its own initial mean. Peak vorticity, energy and enstrophy curves are normalized by their own initial values; I(t) remains absolute. Absolute quantities remain in the raw diagnostics. This is a parameter-control experiment, not a disturbed/reference separation test.

| Run | Preserved dt | W(0.40)/W(0) | E(0.40)/E(0) | Z(0.40)/Z(0) | Failed spatial screens |
|---|---:|---:|---:|---:|---:|
| [angle_degrees-n128-v5.0](diagnostics/angle_degrees-n128-v5.0.json) | 3.448276e-05 | 3.611467 | 0.7506145 | 34.153 | 41/41 |

![Parameter curves](parameter-curves.png)

## Existing controls and numerical limits

The matching completed +5-degree grid-64 control and completed zero-angle grid-128 baseline (both viscosity 0.001) are included at identical saved times. Completed baseline grid-64 and half-step records are also supplied as context. Each retains its own recorded viscosity, initial field and time step. The 64/128 endpoint comparison below does not establish spatial convergence; baseline half-step checks do not establish time convergence of this tilted trajectory.

![Matching grid controls at the endpoint](grid-endpoints.png)

**Raw strain is nonsmooth across periodic joins, and the initial maxima lie outside the central tube.** The published spatial screens require high-band energy below 0.1%, high-band enstrophy below 1%, and minimum half-peak chord at least six cells. These are screening criteria, not proof of resolution. The table distinguishes spatial limitations from budget consistency.

| Run | Final global width, cells | Final central width, cells | Final high-band Z | Maximum absolute E budget residual | Maximum absolute Z budget residual |
|---|---:|---:|---:|---:|---:|
| axis tilt +5 degrees | 1.6165 | 1.431 | 67.98017% | 0.000106 | 0.000485 |

![Budgets and resolution checks](budgets-resolution.png)

Individual budgets use recorded Heun-stage integrals, checked against saved endpoint energies and enstrophies. Independent physical-space checks cover norms, means, divergence, enstrophy production and viscous losses. All declared numerical hard stops pass; neither small residuals nor a viscosity change proves a physical stability transition. Unresolved transfer to small scales remains a limitation.

![Final shell spectra](spectra.png)

The initial-return diagnostic measures distance from each trajectory's own initial velocity field. It does not establish recurrence, periodicity or a search over all state pairs. This unperturbed strain-axis sweep does not supply disturbance-growth or Lyapunov measurements.

![Initial-return diagnostics](initial-return.png)

## Verification and reproducibility

The [verifier](verify_matrix.py) runs as `python -B verify_matrix.py STUDY_DIRECTORY JOB_ID OUTPUT_JSON`. It checks all 41 saved fields per new run, the exact prescribed initial field, unchanged source/protocol hashes, mean, divergence, physical norms, stretching and losses, saved budget endpoints, spectra, widths and declared hard stops. Width and spectral remeasurement uses the frozen routine and is a consistency check, not a second discretization. Original large Fourier fields are required locally and are not embedded; their SHA-256 values are retained below.

- [angle_degrees-n128-v5.0 independent verification](verification/angle_degrees-n128-v5.0.json)

[Raw diagnostics](diagnostics/) include the newly completed strain-axis run and existing completed controls. Existing controls are reused without rerunning or claiming a new independent full-field verification. [Analysis](analysis.json) records summary values and the current completion snapshot. [Frozen protocol](protocol.json), [numerical source](reproduction/numerics.py), [matrix controller](reproduction/run_suite.py) and [publication builder](build_publication.py) retain provenance.

At this report's snapshot, the original matrix has **47/48 complete**, **1 active**, **0 queued**, and **0 recorded failures**. Active jobs: baseline-n256-half. This is a completed batch report; the original matrix's final report awaits its remaining jobs. Historical intake reviews remain unchanged.
