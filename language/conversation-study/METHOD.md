# Method and continuing the study

[← Study home](README.md)

## Unit of observation

Read **prior context → request → answer → clarification → next answer**. Preserve the source order and distinguish speaker turns from quotations, displayed code, interface text, and AI reconstructions. The available material is partial; later speaker boundaries have not all been certified.

The `E` records are 53 selected episodes, sometimes spanning multiple turns. The `P` paths and `C` comparisons reuse those passages. Adding another view does not add another independent conversation event. Categories can overlap. They are reading aids, not verified class labels.

## How confusion is tracked

| Track separately | Evidence to record |
| --- | --- |
| Question still forming | The participant explicitly says they are finding the question. |
| Reference unclear | A word such as “it” could point to more than one earlier passage. |
| Request interpreted differently | Compare the requested operation with what the answer actually supplies. |
| Earlier framing carried forward | Locate the prior wording and its later reuse. |
| Repair | Record the correction, next answer, and what changed. |
| Remaining uncertainty | Keep “partial,” “unclear,” or “unconfirmed” available when the text does not establish resolution. |

A request for shorter columns can be satisfied while the intended relationship is still misunderstood. A changed question is not necessarily a misunderstood question. Agreement, length, and confident language each need context. Absence of acceptance alone is not proof of failure.

## Data and views

| File | Role |
| --- | --- |
| `sources/source-01.txt` through `source-03.txt` | Unchanged supplied source snapshots; byte hashes and line ranges in the manifest. |
| `manifest.json` | Provenance, overlap, coverage, and limits. |
| `evidence.json` | Category definitions and selected observations. |
| `word-paths.json` | Short exact excerpts and interpreted transitions. |
| `confusion-tracks.json` | Context, request, answer, next clarification, and limits. |
| `record-index.json` | Earlier exchange index and candidate alignment of embedded reconstructions. |
| `index.html`, `sources/*.html` | Saved interactive reading edition and escaped source displays. |
| `build.py` | Rebuilds Markdown reading pages and the category SVG from saved data. |
| `validate.py`, `tests/` | Integrity and evidence checks. |

The text at source-03 lines 1–565 repeats source-02 lines 777–1341 exactly. Both source files remain unchanged. The overlap is explicitly recorded and does not become 565 new utterances. AI-authored reconstructions remain reconstructions even when their wording matches an earlier paste.

## What is checked

The validator checks source hashes and line totals, overlap, unique IDs, source spans, category and episode references, exact quoted fragments, the identical request, and the two unchanged quoted premises. It also checks local HTML links and generated source displays. Regression tests deliberately alter a quote, a source, an ID, an overlap span, and the repeated request to verify that damage is detected. `build.py --check` detects stale GitHub views.

Passing checks means the evidence trail is intact. It does not establish the correctness of an interpretation, a scientific analogy, a statement in the conversation, or a causal explanation.

## Limits of this edition

- One evolving conversation, supplied as three partial snapshots, with missing history and embedded repeats.
- Selected manual annotations prepared with AI assistance; no independent inter-rater study has been completed.
- No full-conversation frequency estimates, speaker totals, diagnostic scores, or classification benchmark.
- Four categories were introduced later; older episodes have not been reassessed for them. A blank cell means no tag was assigned, not that the phenomenon is absent.
- An identical repeated request provides an observed comparison, not a controlled replay. Topic, wording, and other intervening context changed together.
- Claims about medicine, hidden system behavior, motives, and technical mechanisms remain claims in the source. The saved report contains limited citation checks, not a full fact-check of every passage.

## Add another paste without duplicating work

1. Save a new source snapshot with a new ID; preserve its bytes and record its hash. Do not overwrite an earlier snapshot.
2. Locate exact overlaps and changed versions before adding observations. Preserve uncertain matches as uncertain.
3. Link new evidence to an existing `E`, `P`, or `C` record when it belongs there. Add an ID for a genuinely new observation.
4. Keep short definitions for new categories and record which episodes were assessed. Do not turn unchecked cases into zeros.
5. Record wording, speaker in context, the requested operation, the answer's interpretation, the next clarification, and the visible outcome.
6. Review the interpretation, rebuild the reading pages, and run the integrity checks. If changing the interactive edition, update that snapshot deliberately as well.

The next useful review is to check a small set of context comparisons with a second reader, recording disagreements and ambiguous cases. Additional conversations can broaden the collection after each retains its own provenance and category version. There is no requirement that they produce the same patterns.
