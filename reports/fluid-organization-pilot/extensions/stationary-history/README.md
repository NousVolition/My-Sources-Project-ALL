# Stationary fluid history experiment — protocol registration

**Status: Models frozen; independent test integrations have not started. There is no held-out prediction result in this commit.**

This study tests whether causal histories of partial observations improve prediction of the response to localized impulses in a continuously forced, three dimensional incompressible fluid. It extends the [completed transient-flow experiment](../stress-history/README.md).

The independent design has 24 training, 8 validation and 32 held-out baseline fields. Each baseline supplies four interventions, with two impulse strengths, and all eight observations remain together in the uncertainty calculation. The model family, input comparisons, stationarity gates and numerical requirements are specified in [protocol.json](protocol.json).

Main settings are a periodic cube of side 6, viscosity 0.025, constant input power 0.2, a 64-cubed grid and time steps 0.02 for the baseline and 0.01 for disturbed copies. Four predeclared test runs have matched half-step and 80-cubed controls. This targets a higher, but still modest, Reynolds number; the achieved stationary-ensemble Reynolds number is a measurement to be reported later.

History must outperform richer and global current-state models, wrong-run histories, and training-only averages by strength and by time/strength. Independent runs, not individual snapshots, are the sampling units. The full instantaneous fluid field determines its future under this forcing law, so any partial-observation history benefit would not imply an additional physical memory force.

## Completed and pending work

Solver and prediction-design tests have been executed locally. Long screening integrations and independent calibration are separate from the prospective prediction test. A completed 64-cubed screening field supplies the response-resolution pilot. Their full observations, numerical receipts and provenance will accompany the final result.

The independent prediction result, ensemble stationarity assessment, final resolution summary and report remain pending in this registration. The completed result will be added to [Streams and Rocks](https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization), with runnable code and recordings in this source repository. Laboratory validation and substantially higher Reynolds numbers remain future experiments.

## Execution

Use Python 3.12 and install requirements.txt. Run the two test files with pytest. Calibration is run with `python pipeline.py calibration --workers 8`; it refuses to freeze a model until the required screening, stationarity and numerical gates pass. Confirmation verifies the frozen source and model hashes before integrating new test fields. This freeze release includes the exact calibration and screening observations required to verify the model hashes. Run `python restore_recordings.py` before checking or using the frozen models. The full result audit additionally requires the later completed test observations and result manifest.

The protocol can be developed during calibration. The final frozen source and model commit must precede all held-out integrations. Git history records intermediate changes; only the final pre-test freeze is the confirmatory specification.
