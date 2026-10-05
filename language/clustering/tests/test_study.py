"""Checks for the study's concrete risks: invalid counts and label leakage."""
import copy
import csv

import numpy as np
import pytest

from study import FEATURES, feature_matrix, fit_partition, load_rows, partition_change


def example_rows():
    return [dict(turn=str(i + 1), ai_words=str(20 * (i + 1)), primary_mode=f"label {i}",
                 **{key: str((i + j) % 4) for j, key in enumerate(FEATURES)}) for i in range(6)]


def write_input(tmp_path, rows):
    path = tmp_path / "table.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_descriptive_labels_and_prose_cannot_change_features_or_partition():
    rows = example_rows()
    modified = copy.deepcopy(rows)
    for row in modified:
        row.update(primary_mode="completely different", opening="arbitrary prose", turn="999")
    _, original, _ = feature_matrix(rows, "counts")
    _, changed, _ = feature_matrix(modified, "counts")
    np.testing.assert_array_equal(original, changed)
    model_a, _, _ = fit_partition(original, 2, 0)
    model_b, _, _ = fit_partition(changed, 2, 0)
    assert partition_change(model_a.labels_, model_b.labels_)["adjusted_rand_index"] == 1.0


@pytest.mark.parametrize("field,value", [("ai_words", "0"), ("apology", "-1"),
                                         ("validation", "nan"), ("clinical_pr", "1.5")])
def test_invalid_observations_fail_instead_of_becoming_data(tmp_path, field, value):
    rows = example_rows()
    rows[0][field] = value
    with pytest.raises(ValueError):
        load_rows(write_input(tmp_path, rows))


def test_duplicate_ids_cannot_silently_merge_review_records(tmp_path):
    rows = example_rows()
    rows[1]["turn"] = rows[0]["turn"]
    with pytest.raises(ValueError, match="unique"):
        load_rows(write_input(tmp_path, rows))


def test_missing_feature_is_not_silently_treated_as_zero(tmp_path):
    rows = example_rows()
    for row in rows:
        del row["apology"]
    with pytest.raises(ValueError, match="Missing columns"):
        load_rows(write_input(tmp_path, rows))


def test_equal_phrase_rates_have_equal_representations_despite_length():
    rows = example_rows()[:2]
    for i, row in enumerate(rows):
        row["ai_words"] = str(100 * (i + 1))
        for key in FEATURES:
            row[key] = str(i + 1)
    rates, scaled, _ = feature_matrix(rows, "per_100_words")
    np.testing.assert_array_equal(rates[0], rates[1])
    assert np.isfinite(scaled).all()


def test_cluster_renumbering_is_not_reported_as_membership_change():
    result = partition_change([0, 0, 1, 1], [8, 8, 3, 3])
    assert result["adjusted_rand_index"] == 1.0
    assert result["changed_pairs"] == 0
    changed = partition_change([0, 0, 1, 1], [0, 1, 0, 1])
    assert changed["changed_pairs"] == 4


def test_too_few_distinct_vectors_cannot_produce_fake_clusters():
    with pytest.raises(ValueError, match="distinct vectors"):
        fit_partition(np.zeros((8, 3)), 3, 0)
