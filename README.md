# Nous Volition

**Exploring patterns in fluid motion and human–AI conversations.**

This personal research project turns visual ideas and observations into diagrams, code, and recorded comparisons. The pages below explain what was tried, show the evidence, and identify questions that remain open.

## Start here

You can read the studies on GitHub without installing anything.

| Explore | The question | Start reading |
| --- | --- | --- |
| **Fluid motion** | How does a smooth, ring-shaped flow change over time, and how much do the computer's calculation settings affect the answer? | [The math in plain language](math/README.md) |
| **Conversation patterns** | What repeats or changes across replies, and how does earlier wording carry into later answers? | [A concrete conversation example](language/conversation-study/CONTEXT.md#c01) · [See categories without the reply text](language/exploration/README.md) |

## Current math result

The **“hug”** is the name for a soft enclosing shape used to prepare the flow's starting pattern. The simulation then follows that flow using the Navier–Stokes equations with zero external force.

The latest calculation continued three saved runs from **model time 0.32 to 0.4**. The finer-grid comparison gives **0.0875% in velocity** and **0.3710% in its spatial gradient**. Halving the time step on 384³ gives **0.0491%** and **0.1245%**, respectively.

**Latest check:** traced energy through 31 individual calculation steps. The positive finite-step contribution is partly removed by pressure correction; their remainder is larger than transport work in every sampled step. Individual energy accounting closes to roundoff. [Measurements and limits](math/notes/hug-boundary.md#one-step-energy-trace-through-model-time-04).

These are measurements over a limited simulated interval. The Clay problem still requires a mathematical proof.

**[See the comparison and chart →](math/notes/hug-boundary.md#continuation-to-model-time-04)**

## Explore the conversations

The [conversation study](language/conversation-study/README.md) links observations about confusion, repeated wording, and changed interpretations back to the supplied text. One example follows the same comparison request receiving different interpretations after different preceding exchanges.

The [earlier category map](language/exploration/README.md) shows reply order, overlapping categories, and reply length across 68 replies. The two source collections are documented separately.

## What is finished, and what comes next?

- **Completed:** the smooth starting-field construction and energy calculation; saved grid and time-step comparisons; conversation maps with evidence links and checks.
- **Latest math check completed:** [one-step energy tracing](math/notes/hug-boundary.md#one-step-energy-trace-through-model-time-04); earlier energy, grid, evolution and shape results are preserved.
- **Next conversation work:** review interpretations with another reader and add new material while recording overlaps.

[Full project status and completed work →](STATUS.md)

## Code, papers, and history

- [Math code and instructions](math/README.md#use-the-code)
- [Conversation study and its checks](language/conversation-study/README.md#run-the-checks)
- [Working papers and their source status](papers/README.md)
- [Earlier work and consolidation record](archive/README.md)

This repository holds the current code, reports, and status. [Nous-Volition](https://github.com/NousVolition/Nous-Volition) preserves the original stream-function report. Check [project status](STATUS.md) before starting work so completed experiments can be reused.

<details>
<summary>Run a small example locally</summary>

Use Python 3.12 from this repository's root:

```sh
python -m pip install -r requirements.txt
python -m pytest -q
python language/conversation_map.py --output scratch/conversation-map
python math/box_experiment.py --profile smooth --start projected --points 8 --steps 1
```

The measurements from completed experiments are already saved with their reports.

</details>
