# Fluid Organization Under Stress: completed pilot

## Completed recovery and network extension — 10 October 2026

[Read the new experiments, results and code links](extensions/recovery-network/README.md). This separate extension adds 36 three-arm fluid recovery runs, a specified three-state network with frozen learning controls, independent test runs, and 14 passing automated tests. It does not alter the original pilot results. The intervention-specific fluid history benefit remains unresolved.

Does past marker organization help predict a fluid's response to a local disturbance?

This is a physical-fluid pilot, separate from the social SIMS studies. Completed work comprises **93 paired fluid configurations**, **7 particle-tracking configurations**, and **22 passing automated tests**, with separate phase-locking, linear-system, pendulum, relaxation, weakly nonlinear, parametric-forcing and bifurcation controls.

## Main result

Across eight held-out initial conditions, adding past marker geometry reduced mean per-run RMSE for overall future neighborhood deformation from **0.07082 to 0.02040 (71.19%)**. The paired whole-run bootstrap interval for the RMSE difference was **[-0.05664, -0.04431]**. Shuffling histories removed the benefit, and the finer-grid check preserved it.

For the **additional deformation caused by the disturbance**, history did not demonstrate an improvement: RMSE changed from **0.006877 to 0.006897**, with an interval for the difference spanning zero. The first result therefore does not establish the main disturbance-response hypothesis. History can help infer state missing from an instantaneous local observation; it is not evidence for a new physical memory force.

## Scope and numerical limits

- Three-dimensional incompressible Navier–Stokes in a periodic box, smooth matched initial conditions, identical baseline/disturbed copies and a localized divergence-free impulse. Markers are passive, not interacting molecules.
- Assigned water and air viscosities, forcing and speed controls, independent random starts, half-time-step and finer-grid comparisons. A Stokes ablation is explicitly a mathematical control.
- Fog and idealized smoke are dilute inertial particles carried by air, with gravity present/removed. Fixed droplets, no evaporation, condensation, combustion or two-way coupling.
- Perturbation energy can grow by transfer from background flow while total unforced kinetic energy falls. Marker separation alone does not establish amplification.
- The original fast-water coarse-grid marker comparison **failed**: RMS change 0.03482 against a 0.005 screen. For one seed, the 50-to-62-grid continuation passed at 0.003797. This does not validate every stress realization or prove continuum convergence.
- Ice, phase changes, wall-boundary comparisons and broader stressful regimes remain proposed. The finite-dimensional examples are established mathematical controls, not evidence of a fluid bifurcation or a solution to the Clay problem.

## Report and reproduction

The presentation is part of [Streams and Rocks](https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization). This directory is the code and complete-recording source linked from that study.

- [Full report](https://github.com/NousVolition/Nous-Volition/blob/main/studies/fluid-organization/report.html) — download this folder and open locally; GitHub displays HTML source.
- [Run instructions and file guide](README.txt)
- [Parameter regimes](results/parameter_regimes.csv), [all fluid configurations](results/all_runs.csv), [prediction results](results/prediction.json), [numerical controls](results/numerical_controls.csv)
- [Saved-data verification](results/verification.json), [22-test record](results/tests.xml), [original snapshot SHA-256 manifest](manifest.json)
- [Original protocol](protocol.json), [subsequent additions](protocol_addendum.json), [upstream provenance](UPSTREAM.txt)

All raw paired marker trajectories, initial/final fluid states, diagnostics and prediction outputs are included in the checksum-verified [recorded package](recorded-package.json). Run `python restore_recordings.py` once after cloning to restore the complete original folder, including all numerical arrays and the offline report. No network or extra dependency is needed for restoration. The original delivered files retain their exact bytes; the original manifest describes those files. This README, Git attributes/ignores, the archive-parts index, restoration script and `publication.json` are publication additions.

### Reproduce

Use Python 3.12 and run from this directory:

```text
python restore_recordings.py
python -m pip install -r requirements.txt
python -m pytest test_pilot.py -q --junitxml=results/tests.xml
python run_experiments.py --stage all
python particles.py
python prediction.py
python phase_locking.py
python dynamics.py
python pendulum.py
python relaxation.py
python weak_nonlinear.py
python parametric.py
python bifurcations.py
python extend_stress.py
python analyze.py
python verify.py
```

The fluid runner checks source, configuration and data hashes before reusing saved runs. A changed solver needs a fresh data directory; see [README.txt](README.txt). Git attributes preserve original bytes across platforms.

The upstream Fourier operators come from `NousVolition/My-Sources-Project-ALL` at commit `1b12558739736d472b0eb75845cf721371fcbb71`. The previously published history study has a different target and numerical limitations; its effect size is not directly comparable to this pilot.

## Added heteroclinic-network reference

[Identified source, nine-saddle diagram, oscillator equation and proposed tests](network_reference_addendum.html). This reference addendum is not an additional completed simulation.
