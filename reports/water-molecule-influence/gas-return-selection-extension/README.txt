SELECTING OLDER PASSAGES WITH PAST-ONLY VELOCITY

Open report.html. This is an analysis of existing recorded trajectories, not
a new molecular simulation. Method definitions were written in data/protocol.json
before new results were computed. Start 1 selects shared weights; starts 2 and 3
evaluate them. Those starts were inspected in prior work, so this is exploratory
separate-start validation, not untouched confirmatory evidence.

Run with Python, NumPy and Matplotlib:
  python analyze_selection.py
  python build_report.py

Input paths and hashes: data/results.json.
Frozen weights and all training scores: data/locked_training.json.
Per-anchor outputs: data/predictions_r*_gamma*.npz.
Each output stores each method's forecast and squared error at horizons
0.02, 0.05 and 0.10 ps. The three analogue selectors also store their old times,
initial match measurements and absolute paired-path distances and angles.
Velocity uses a backward difference. Future data cannot influence the choice
of old passage or its translated-continuation forecast.

The primary error is mean squared two-dimensional position error, not an
independent-event success probability. See the report for baseline failures,
longer horizons and finer-recording sensitivity.
