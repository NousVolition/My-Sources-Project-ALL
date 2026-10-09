# One lean: reproduced ball and fluid series

**The supplied `one_lean.py` reproduces all four series in `one-table.json`.** Both ball series match exactly. The fluid values agree within 3.64e-16.

| Example | Starting value | Final value | Final time |
| --- | ---: | ---: | ---: |
| Ball, one low point | x = 0.2 | x = 0.00863322 | 8 |
| Ball, two low points | x = 0.2 | x = 1.07347523 | 8 |
| Even fluid | E near roundoff | E near roundoff | 0.40 |
| Leaning fluid | E = 0.01999900 | E = 0.02143760 | 0.40 |

The ball's signed position x and the fluid's unsigned mirror error E are recorded as different measurements. The supplied script does not fit a numerical mapping between them.

## Executed rules

The ball starts at x = 0.2 with zero velocity. Both examples use drag 0.45 and dt = 0.02. The update changes velocity first, then position. The one-low-point force is -x; the two-low-point force is x - x^3.

The fluid uses the same 33-grid mirrored start, perturbation template, viscosity 0.01, Heun step and prescribed per-step breathing push as the verified six-case experiment. That push is not multiplied by dt, so this is driven evolution with timestep-dependent driving.

The fluid's midpoint E is 0.01768094 at time 0.20, before its later rise. [The full step-by-step analysis and actual fluid animations](../timeline/README.md) show why: the absolute mirror difference and the overall flow strength change at different rates.

## Files and reproduction

- [Supplied script](one_lean.py), [solver](clay_hug.py), [supplied reference](reference-one-table.json), [reproduced series](one-table.json), [verification](verification.json).
- Run `python reproduce.py` from this folder after installing `requirements.txt`. New results go into `rerun/`.

This verifies the included time-8 ball examples and time-0.40 fluid examples. The separate longer-time, strength/place, and removal-sweep JSON files in the uploaded bundle are not reproduced by this script. The supplied `plot_steps.py` renders their reported numbers but does not evolve those simulations.

[Breathing overview](../README.md)
