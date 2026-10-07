# Nous Volition

**One home for the math, studies, and work already completed.**

Start with a section below. [What is done and what comes next →](STATUS.md)

| Explore | Start here |
| --- | --- |
| **Math & fluid experiments** | [Smooth initial field, derivation, code, and saved results](math/README.md) |
| **Conversation patterns** | [Explore the conversation without the reply text](language/exploration/README.md) |
| **Conversation study** | [Sources, confusion, context comparisons, and evidence checks](language/conversation-study/README.md) |
| **Working papers** | [The two papers and their source status](papers/README.md) |
| **Earlier work** | [Original drafts, Copilot branches, and run history](archive/README.md) |

## Current focus

The [conversation study](language/conversation-study/README.md) is now a separate, source-linked collection: 53 observations, eight word paths, and six context comparisons across three supplied snapshots. Its reading pages work directly on GitHub. The original text, annotations, interactive edition, and checks are together in the folder. It remains distinct from the earlier 68-reply table below.

Explore how a conversation changes over time: repeated categories, overlaps, unmarked replies, and reply length. The [new map](language/exploration/README.md) preserves the existing annotations and sequence. It is an exploratory view that can grow as more conversations are collected.

The smooth initial field and its energy calculation are finished. The existing hug now shapes the ring's starting field, giving smooth matching boundaries without external force. [Read the hug construction, grid comparisons, and time-step checks →](math/notes/hug-boundary.md)

The saved study now includes [a 384³ comparison at time 0.16](math/notes/hug-boundary.md#finer-grid-at-time-016), alongside the earlier 32³, 64³, 128³, and 256³ grids. [Three saved runs now continue to model time 0.16](math/notes/hug-boundary.md#longer-run-to-model-time-016), including the 256³ half-time-step control. Code, checks, numerical records, and figures are together in [math/](math/README.md). [What is finished and what comes next →](STATUS.md#next-math-step)

**Latest:** [The 384³ half-step comparison](math/notes/hug-boundary.md#half-time-step-on-the-384-grid) is complete. Starting from exactly identical arrays, time steps **0.001** and **0.0005** give **0.0242%** velocity difference and **0.0840%** full-gradient difference. The hug construction, equation, viscosity, and zero external force are retained. The previous full-step run was reused.

## Run the existing work

Use Python 3.12 from this repository's root:

```sh
python -m pip install -r requirements.txt
python -m pytest -q
python language/conversation_map.py --output scratch/conversation-map
python math/box_experiment.py --profile smooth --start projected --points 8 --steps 1
```

The [full comparison results](math/results/box-experiment-results.csv) are already saved; the short command above is just a way to try the runner.

## How the two repositories fit

This repository maintains the code, corrected note, and project status. [Nous-Volition](https://github.com/NousVolition/Nous-Volition) keeps the original stream-function report and points here for current work.

Check [STATUS.md](STATUS.md) before beginning a task. Update the existing implementation when extending it. Prior branch work is traced in the [consolidation record](archive/CONSOLIDATION.md).
