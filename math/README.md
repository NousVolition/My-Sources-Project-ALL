# Math & fluid experiments

[← Project home](../README.md) · [Done & next](../STATUS.md)

## Read

1. **[The smooth initial field](notes/smooth-initial-field.md)** — the corrected construction and exact energy.
2. **[The cube experiments](EXPERIMENTS.md)** — saved results, pressure correction, and what the tests establish.
3. **[Review of the original framework](notes/math-review.md)** — mathematical claims still requiring work.
4. **[Norm and Entropy: two mirrored halves](notes/norm-entropy-mirror.md)** — paired samples, K-means before and after alignment, and the separate boundary check.
5. **[A soft envelope: hug, imprint, and pressure](notes/soft-envelope.md)** — small motion and memory prototypes, with explicit opening rules.
6. **[The hug closes the ring's boundary](notes/hug-boundary.md)** — the actual closed hug shapes the maintained ring's stream function, giving a smooth, divergence-free start with matching cube faces. [Finer-grid evolution](notes/hug-boundary.md#finer-grid-evolution-check) follows the same field through a short unforced run; [the longer interval](notes/hug-boundary.md#longer-run-to-model-time-016) continues three saved states to time 0.16. [The 384³ comparison](notes/hug-boundary.md#finer-grid-at-time-016) adds spatial refinement at that same final time.
7. **[Separate trigonometric illustration](notes/hidden-flow-evolution.md)** — assistant-chosen fields with an analytic reduction and refined numerical runs; these do not evolve the project's ring field or establish the proposed conversation-to-fluid mapping.

![Original and corrected ring profiles](figures/smooth-initial-profile.png)

The right panel zooms into the center and shows `value / center value - 1` for each curve. The figure uses equal parameter numbers, which give different ring widths; the cube experiments separately match the local peak widths.

The editable derivation is [smooth-initial-field.tex](notes/smooth-initial-field.tex). The [original report](https://github.com/NousVolition/Nous-Volition/blob/main/stokes_stream_report.pdf) remains in the companion repository.

## Use the code

| File | Role |
| --- | --- |
| [initial_field.py](initial_field.py) | The single maintained implementation of the smooth field and its energy |
| [stokes.py](stokes.py) | The original Gaussian, kept for reproducing historical comparisons |
| [navier.py](navier.py) | The original cube update plus the corrected discrete projection |
| [box_experiment.py](box_experiment.py) | Configurable comparisons for the original, smooth, and hug profiles and all box sizes |
| [verify_energy.py](verify_energy.py) | The existing independent quadrature check |
| [reflection_probe.py](reflection_probe.py) | Face matching, slope reversal, and divergence for a flipped neighboring copy |
| [mirror_family.py](mirror_family.py) | Norm/Entropy pairing and the illustrative K-means comparison |
| [hug_envelope.py](hug_envelope.py) · [pressure_envelope.py](pressure_envelope.py) | Shared mirrored-arm geometry, timed release, and a separate pressure-triggered opening prototype |
| [hugged_ring.py](hugged_ring.py) | The existing closed hug applied to the existing smooth ring through its stream function; used by `box_experiment.py --profile hug` |
| [hug_refinement.py](hug_refinement.py) | Grid and time-step comparisons of the same hugged ring, using the existing cube stencils; caches completed runs |
| [extend_hug_refinement.py](extend_hug_refinement.py) | Extends the grid study; `--half-step` checks 256³ at time 0.08, and `--finer-grid 384` compares the saved 256³ field with 384³ at time 0.16, preserving checkpoints |
| [continue_hug_refinement.py](continue_hug_refinement.py) | Continues the saved 128³ and 256³ states to time 0.16, preserving the original step sizes and zero force; checks restart identity |
| [hidden_flow_experiment.py](hidden_flow_experiment.py) | Unforced periodic Fourier evolution for A/B/C and the three-dimensional continuation |
| [tests/](tests/) | Consolidated checks from the repositories and this review |

From the repository root:

```sh
python -m pip install -r requirements.txt
python -m pytest -q
python math/box_experiment.py --profile smooth --start projected --box 18 --points 16
```

Read the [saved CSV](results/box-experiment-results.csv) for the completed full suite. To intentionally reproduce it, use `python math/box_experiment.py --suite --csv scratch/box-results.csv`.

The field construction and the fluid evolution have separate scopes: the analytic initial field is established; the cube is an experimental discretization with known boundary and centered-grid limitations. [Next step →](../STATUS.md#next-math-step)
