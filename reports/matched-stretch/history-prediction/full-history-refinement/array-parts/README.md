# Exact final prediction arrays

GitHub rejected the single46 MB archive through its blob API. These three binary parts concatenate to the original archive exactly; [delivery.json](delivery.json) records each part and the complete checksum. They contain the same saved outputs and features, with no numerical recomputation.

From the repository root, run:

```text
python -B reports/matched-stretch/history-prediction/full-history-refinement/restore_prediction_arrays.py
```

The tool checks every part, recreates both final archive names and verifies the original SHA-256 6bd1090182f84b28e5725403956aec33bfb57f6e6d9ea396c24683ca462c8b68. It refuses to overwrite differing existing files. Then run the saved-data verification/reproduction commands in the report.
