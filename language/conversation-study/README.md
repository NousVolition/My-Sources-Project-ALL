# Conversation study

**What changes between a question, an answer, and the next attempt to be understood?**

An exploratory study of an evolving human–AI conversation. The starting interest is simple: look at the categories and paths through a conversation, then return to the wording to see what fits. There is no required outcome or fixed classification target.

**Three supplied snapshots · 53 selected observations · 13 overlapping categories · 8 word paths · 6 context comparisons**

## Start here

| To explore | Open |
| --- | --- |
| A quick example | [The same request receives different interpretations](CONTEXT.md#c01) |
| Confusion and repair | [Prior context → request → answer → clarification](CONTEXT.md) |
| Categories without the reply text | [Category map and definitions](CATEGORIES.md) |
| Wording that changes direction | [Eight paths with exact excerpts](PATHS.md) |
| Every selected observation | [Evidence and reading limits](OBSERVATIONS.md) |
| The original material | [Three unchanged text snapshots and overlap record](SOURCES.md) |
| How to inspect or extend it | [Method, limitations, and checks](METHOD.md) |

These pages render directly on GitHub. For search, filters, expandable evidence, and numbered source pages, download this repository, unzip it, and open `language/conversation-study/index.html` in a browser. Keep its `sources/` folder beside it. No server, account, or installation is needed. GitHub's file view displays HTML source rather than running the interactive report.

## One concrete result

The **complete comparison request is identical** at [source-03 line 1898](sources/source-03.txt#L1898) and [line 2210](sources/source-03.txt#L2210). The two sentences quoted in the answers are also unchanged.

| What came before | What the answer does |
| --- | --- |
| Variables, distance, thresholds | Assigns the sentences roles in a cause-and-effect engine. |
| Layer, skin, counting the middle | Explicitly applies the recent skin metaphor, making one sentence the outer layer and the other the activity inside. |

The user then rejects that interpretation. This is an inspectable example of earlier wording carrying into an answer. It does **not** isolate history as the sole cause of the difference. [Read the full comparison and its limits →](CONTEXT.md#c01)

## What this project contributes

- An evidence trail from each observation to its original source lines.
- A way to distinguish uncertainty, ambiguous references, changed interpretation, and partial repair.
- Separate tracking of presentation improvements and whether the intended meaning was understood.
- Preserved overlap and reconstruction status, so repeated text is not counted as independent evidence.
- Executable checks for source integrity, references, quoted fragments, and the repeated request.

The contributor supplied the conversations and guided the questions. Organization, annotations, and tooling were prepared with AI assistance. Readings remain provisional and inspectable. The conversation records statements, questions, quotations, and experiments; it does not authenticate personal facts or the assistant's explanations of its internal behavior.

## Run the checks

From the repository root, using Python 3.12 or later; these checks need only the standard library:

```sh
python language/conversation-study/validate.py
python -m unittest discover -s language/conversation-study/tests -v
python language/conversation-study/build.py --check
```

To regenerate the GitHub reading pages from the saved annotations:

```sh
python language/conversation-study/build.py
```

The interactive HTML is the preserved reviewed snapshot. The builder regenerates the GitHub Markdown and category graphic; it does not re-annotate the conversation or rebuild the interactive HTML.

## Where this fits

This folder belongs to the language work in **My-Sources-Project-ALL**. It is a separate source collection from the older 68-reply conversation map; their rows and categories have not been merged.

[Project home](../../README.md) · [Organized language work](https://github.com/NousVolition/My-Sources-Project-ALL/tree/copilot/put-stokes-start-on-cube/language) · [Earlier 68-reply map](https://github.com/NousVolition/My-Sources-Project-ALL/tree/copilot/put-stokes-start-on-cube/language/exploration)
