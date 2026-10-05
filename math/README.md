# Math & fluid experiments

[← Project home](../README.md) · [Done & next](../STATUS.md)

## Read

1. **[The smooth initial field](notes/smooth-initial-field.md)** — the corrected construction and exact energy.
2. **[The cube experiments](EXPERIMENTS.md)** — saved results, pressure correction, and what the tests establish.
3. **[Review of the original framework](notes/math-review.md)** — mathematical claims still requiring work.

![Original and corrected ring profiles](figures/smooth-initial-profile.png)

The editable derivation is [smooth-initial-field.tex](notes/smooth-initial-field.tex). The [original report](https://github.com/NousVolition/Nous-Volition/blob/main/stokes_stream_report.pdf) remains in the companion repository.

## Use the code

| File | Role |
| --- | --- |
| [initial_field.py](initial_field.py) | The single maintained implementation of the smooth field and its energy |
| [stokes.py](stokes.py) | The original Gaussian, kept for reproducing historical comparisons |
| [navier.py](navier.py) | The original cube update plus the corrected discrete projection |
| [box_experiment.py](box_experiment.py) | One configurable runner for both profiles and all box sizes |
| [verify_energy.py](verify_energy.py) | The existing independent quadrature check |
| [tests/](tests/) | Consolidated checks from the repositories and this review |

From the repository root:

```sh
python -m pip install -r requirements.txt
python -m pytest -q
python math/box_experiment.py --profile smooth --start projected --box 18 --points 16
```

Read the [saved CSV](results/box-experiment-results.csv) for the completed full suite. To intentionally reproduce it, use `python math/box_experiment.py --suite --csv scratch/box-results.csv`.

The field construction and the fluid evolution have separate scopes: the analytic initial field is established; the cube is an experimental discretization with known boundary and centered-grid limitations. [Next step →](../STATUS.md#next-math-step)
