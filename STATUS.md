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
| Language reports | [Report guide](language/README.md) | Existing tables and charts preserved; available chart links repaired |
| Two working papers | [Paper guide](papers/README.md) | Existing documents preserved, with missing source material identified |

**Verification:** the consolidated suite has 80 passing tests in Python 3.12.14 with NumPy 2.3.5 and pytest 9.1.1. The longer box runs reuse their saved results. The standalone LaTeX source is preserved; PDF compilation remains unverified because the editor compiler could not find its platform directories.

## Next math step

Choose the boundary treatment for the evolution study. The current field is defined on all of space; the cube wraps periodically, which creates a seam when the field is still appreciable at its faces.

Then use the **existing** box runner to refine space and time at a fixed physical domain and boundary treatment. Record an error measure against a reference solution or a finer run. The present short runs do not establish convergence, long-time stability, or a singularity theorem.

## Material still missing

- The original pasted inputs for the two language scripts and the generator for the third report.
- The original connection-map image; its existing table is readable.
- The supporting source documents for the California paper and the notebook underlying the boundary paper.
- The video associated with the subtitle file.

The conversation reconstruction is preserved as a reconstruction; it is not substituted for a missing original input.

## Work history

The [consolidation record](archive/CONSOLIDATION.md) identifies which Copilot branches were reused, which duplicate implementations were retired, and which sessions produced no file changes. It is the place to check before restarting an old task.
