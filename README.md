# Nous Volition

**Sources, conversation studies, and mathematical experiments.**

Start with a section below. Each study connects its questions to the available source material, working code, checks, and limits.

## Explore the work

| Area | Start here |
| --- | --- |
| **Conversation study** | [Categories, confusion, wording, and what changes after a correction](language/conversation-study/README.md) |
| **Original conversation material** | [Three supplied snapshots, source links, and overlap notes](language/conversation-study/SOURCES.md) |
| **Earlier conversation map** | [68 replies viewed through overlapping categories](https://github.com/NousVolition/My-Sources-Project-ALL/tree/copilot/put-stokes-start-on-cube/language/exploration) |
| **Math & fluid experiments** | [Smooth initial field, derivation, code, tests, and saved results](https://github.com/NousVolition/My-Sources-Project-ALL/tree/copilot/put-stokes-start-on-cube/math) |
| **Working papers** | [The two papers and their source status](https://github.com/NousVolition/My-Sources-Project-ALL/tree/copilot/put-stokes-start-on-cube/papers) |
| **Project progress** | [What is finished and what comes next](https://github.com/NousVolition/My-Sources-Project-ALL/blob/copilot/put-stokes-start-on-cube/STATUS.md) |

## Latest math: the hug-shaped ring

**[Read the latest 384³ grid comparison](https://github.com/NousVolition/My-Sources-Project-ALL/blob/copilot/put-stokes-start-on-cube/math/notes/hug-boundary.md#finer-grid-at-time-016)** · [Code and tests](https://github.com/NousVolition/My-Sources-Project-ALL/blob/copilot/put-stokes-start-on-cube/math/README.md)

The closed hug shapes the maintained smooth ring's starting field. A new **384³** run is compared with the saved **256³** field at **model time 0.16**, using the same time step **0.001** and **zero external force**.

- 256³ vs 384³: **0.0619% velocity difference**, **0.3509% gradient difference**.
- Earlier 128³ vs 256³: **0.3287%** and **1.8051%**. These grid changes have different refinement ratios (2 and 1.5).
- Six targeted comparison/restart checks passed, including a new check of the 3:2 grid alignment. Four checkpoints preserve the new work; previous completed runs were reused.

The 256³ time-step control and the earlier results remain available. A 384³ half-step check has not been run. These are finite-time numerical comparisons, without a rigorous error bound or an all-time smoothness or breakdown proof.

## Latest study: following a conversation

The new [conversation study](language/conversation-study/README.md) follows **request → answer → correction → next answer**, with 53 selected observations, eight word paths, and six comparisons of prior context. Its category view makes patterns visible before returning to the wording.

One concrete example: the same complete request appears twice, but the answers use different frames after different preceding discussion. [Read the comparison and inspect the original lines →](language/conversation-study/CONTEXT.md#c01)

The GitHub pages are ready to read. The full source snapshots, machine-readable annotations, interactive HTML, and executable checks are together in the same folder. The work is exploratory and prepared with AI assistance; interpretations and source claims are kept distinct.

## Check the study

After downloading or cloning the repository, use Python 3.12 or later from its root. These checks use only the standard library:

```sh
python language/conversation-study/validate.py
python -m unittest discover -s language/conversation-study/tests -v
python language/conversation-study/build.py --check
```

For the interactive view, open `language/conversation-study/index.html` from the downloaded folder in a browser. The reading pages above work directly on GitHub.

## Project organization

The conversation study is available here on the main branch. The links to math, earlier language work, papers, and progress open the organized working version while [the broader reorganization](https://github.com/NousVolition/My-Sources-Project-ALL/pull/1) remains under review. They can all be read without using the pull-request interface.

[Nous-Volition](https://github.com/NousVolition/Nous-Volition) holds the original stream-function report. This repository brings the sources, implementations, and studies together. The [consolidation record](https://github.com/NousVolition/My-Sources-Project-ALL/blob/copilot/put-stokes-start-on-cube/archive/CONSOLIDATION.md) traces earlier work so completed pieces can be reused.


