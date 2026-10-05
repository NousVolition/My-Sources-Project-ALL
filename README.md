# Nous Volition

**One home for the math, studies, and work already completed.**

Start with a section below. [What is done and what comes next →](STATUS.md)

| Explore | Start here |
| --- | --- |
| **Math & fluid experiments** | [Smooth initial field, derivation, code, and saved results](math/README.md) |
| **Language studies** | [Existing reports, charts, and source notes](language/README.md) |
| **Working papers** | [The two papers and their source status](papers/README.md) |
| **Earlier work** | [Original drafts, Copilot branches, and run history](archive/README.md) |

## Current focus

The smooth initial field and its energy calculation are finished. The cube experiments and corrected projection are implemented and checked. Their saved results are ready to read.

The next math step is to choose the boundary treatment, then study accuracy as the grid and time step are refined. [Continue from here →](STATUS.md#next-math-step)

## Run the existing work

Use Python 3.12 from this repository's root:

```sh
python -m pip install -r requirements.txt
python -m pytest -q
python math/box_experiment.py --profile smooth --start projected --points 8 --steps 1
```

The [full comparison results](math/results/box-experiment-results.csv) are already saved; the short command above is just a way to try the runner.

## How the two repositories fit

This repository maintains the code, corrected note, and project status. [Nous-Volition](https://github.com/NousVolition/Nous-Volition) keeps the original stream-function report and points here for current work.

Check [STATUS.md](STATUS.md) before beginning a task. Update the existing implementation when extending it. Prior branch work is traced in the [consolidation record](archive/CONSOLIDATION.md).
