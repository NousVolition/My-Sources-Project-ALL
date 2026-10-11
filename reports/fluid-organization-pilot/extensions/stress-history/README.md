# Stronger fluid disturbances: frozen protocol, study in progress

Published before held-out simulations. The numerical screening and 17 analytic/physical/leakage tests are complete; calibration runs are underway. No held-out prediction result is claimed here.

The [locked prospective protocol](protocol.json) specifies24 training,8 validation and32 independent test fields per viscosity, paired across two viscosities. Each field is copied into unchanged and localized-impulse arms at0.25 and0.75 baseline velocity RMS. Models are frozen before generating the test fields. Strong current inputs include full velocity and velocity-gradient tensors at the impulse center, current marker geometry, global/spectral summaries and initial perturbation-energy transfer rates. The study compares causal histories with this stronger present observation, training-only priors and matched shuffled-history controls, with field-level uncertainty and separate numerical refinements.

The lower-viscosity cases are transient modest-Reynolds-number 3D flows, not demonstrated fully developed high-Reynolds turbulence. Markers are one-way tracers. The retained trilinear marker screen failed; cubic tracking passed finer-grid and half-step revalidation, with identical fluid-energy arrays. See [screen events](screen_events.json), [screen verification](results/screen_verification.json) and [source provenance](UPSTREAM.txt).

Python3.12 or later, from this folder:

```text
python -m pip install -r requirements.txt
python -m pytest test_study.py -q
python pipeline.py
python report.py
python manifest.py
python audit.py
```

Completed results, recordings, plots and a report will be added here and linked from [Streams and Rocks](https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization). The earlier120-field study has a different forcing-change target and remains a separate result.
