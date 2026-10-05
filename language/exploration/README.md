# Conversation without the reply text

[← Language studies](../README.md)

## Purpose

Look at the shape of a conversation: what repeats, where categories overlap, what changes, and where the existing categories leave gaps. This personal exploratory project starts from saved observations. Categories can evolve as more material is collected.

![Eight saved categories and reply length across 68 replies](conversation-map.svg)

Each column is one reply, in its original order. A filled cell means a saved count above zero. The upper bars show reply length. Counts are available in the [wordless data export](map-data.json); the chart uses presence so multiple matches do not dominate its appearance.

## Two kinds of labels

| Layer | What the saved table contains | How this view uses it |
| --- | --- | --- |
| Descriptive notes | 66 distinct `primary_mode` labels across 68 replies | Preserved as individual observations; variety is useful and is not itself a data error |
| Tracked categories | Eight numeric columns such as validation, apology, and institutional | Shown together; a reply can have several marks or none |

The 66 descriptions and the eight count columns serve different purposes. Neither is established ground truth for a classifier. The rare descriptions are still useful notes about individual replies.

## Things worth looking at

These are observations about the **saved annotations**, not new interpretations of the original conversation.

- **Overlap:** four of the five apology-marked replies also have validation marks: replies 24, 49, 50, and 58. Keeping multiple categories on a reply preserves this relationship.
- **A repeated stretch:** validation is marked on four consecutive replies, 9–12.
- **Several marks together:** reply 58 has six of the eight categories. It is a useful place to inspect how the categories fit together.
- **A change in length:** replies 57–60 contain 280, 364, 288, and 281 words. Replies 61–68 contain 46, 8, 77, 75, 122, 7, 7, and 82. The map shows the change; it does not tell us why it happened.
- **Coverage gaps:** 30 of the 68 replies have no recorded matches in these eight categories. That means the current measurements leave them unmarked, not that nothing happened in them.

These are starting points for looking again. One conversation cannot establish how common a pattern is across conversations, people, or models. More data can show which patterns recur when sources and annotation methods remain traceable.

## Reproduce and extend

From the repository root, after installing `requirements.txt`:

```sh
python language/conversation_map.py --output scratch/conversation-map
python -m pytest -q
```

The [generator](../conversation_map.py) reads the existing [third-version table](../data/ai_language_turns_3.csv). It exports only reply IDs, word counts, descriptive labels, and category counts. The source table itself still contains text; this map is a view without reply text, not an anonymization of the repository.

The input's SHA-256 is recorded in `map-data.json`. The third-version annotation generator and original source input are missing, so the counts and word totals have **not** been independently re-extracted. Category names are retained as supplied; their underlying matching rules need recovery or fresh documentation before applying them consistently to new material.

See [Collecting more observations](COLLECTING.md) for the next step. The [optional clustering study](../clustering/README.md) offers another lens on these same counts. It stays separate from the sequence map because K-means does not use reply order.

This extension was developed with AI assistance. The code checks reproduction and representation; interpreting whether a category fits remains a human review task.
