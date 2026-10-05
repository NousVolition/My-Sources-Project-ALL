# Collecting more observations

[← Conversation map](README.md)

Keep looking, keep your notes, and let categories develop. The useful discipline is preserving what each observation refers to.

## Keep three connected layers

1. **Source:** the original conversation, with a conversation ID, date, source, and reply order. Keep a private source private; decide separately what belongs on GitHub.
2. **Observations:** your description of each reply and any category marks. More than one mark is allowed. An unmarked reply remains visible.
3. **Views:** charts and optional computed groups derived from those observations. A new view does not overwrite the source or notes.

## Small record to keep with each reply

| Field | Purpose |
| --- | --- |
| `conversation_id` + `turn` | Stable identity, even when another conversation also has a reply 1 |
| `source_id` | Reference to the original material |
| `primary_mode` | Your existing free-form descriptive note; it need not be one of a fixed set |
| Category counts or marks | Reusable observations that can overlap |
| `annotation_version` | Which definitions or matching rules were used |
| Review note | Uncertainty, a disagreement, or a reason for changing a category |

This is a collection guide, not an implemented multi-conversation importer. The current runner accepts one saved table at a time. Keep new conversations separate until their identifiers and annotation versions can be joined deliberately.

## When a category changes

Write down a short definition, an example that fits, and an example that does not. Record the date/version and keep the original label. If a broad pattern becomes clearer later, add a shared tag alongside the original observation. Do not silently merge different meanings just because their names look similar.

For counts, record the matching rule and how words are counted. Record zero only when the category was checked and no match was found; keep “not checked” distinct. The existing eight numeric columns contain zeros but do not establish whether the annotation coverage was complete.

## What more data can answer

- Does a pattern recur in another conversation?
- Does an apparent overlap persist when the category definition stays the same?
- Which replies remain difficult to describe with the current categories?
- Does a pattern depend on reply length, topic, or where it occurs in the conversation?

More rows create more opportunities to compare. They do not automatically make a category accurate. Preserve the ability to go back from a mark to its source and reconsider it.
