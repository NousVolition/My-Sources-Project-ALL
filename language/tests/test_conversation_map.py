"""The wordless export must preserve chronology and exclude source prose."""
import csv
import json

from conversation_map import FEATURES, map_data


def test_wordless_map_keeps_overlap_and_unmarked_turns_without_reply_text(tmp_path):
    rows = []
    for turn in (3, 1, 2):
        rows.append(dict(turn=turn, primary_mode=f"observation {turn}", ai_words=20,
                         user_words="PRIVATE USER SENTINEL", opening="PRIVATE REPLY SENTINEL",
                         **{key: int(turn == 2 and key in ("apology", "validation")) for key in FEATURES}))
    source = tmp_path / "input.csv"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = map_data(source)
    assert [row["turn"] for row in result["turns"]] == [1, 2, 3]
    assert result["summary"]["no_recorded_matches"] == [1, 3]
    assert result["summary"]["co_occurrence"] == [
        dict(first="validation", second="apology", turns=[2])]
    assert "SENTINEL" not in json.dumps(result)
    assert result["turns"][1]["label"] == "observation 2"
