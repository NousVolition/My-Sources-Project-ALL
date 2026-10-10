"""Production-compatible quench runner skeleton.

Uses the matched-stretch Flow class.
Does not launch long runs; demonstrates the control path.
"""
import hashlib
import json
from pathlib import Path
import numpy as np

# Import the production class (adjust path when run from repo root)
# from reports.matched_stretch.numerics import Flow

def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1<<20), b''):
            h.update(chunk)
    return h.hexdigest()

def energy_budget_residual(E, E0, integrated_viscous_loss):
    return (E - E0 + integrated_viscous_loss) / max(E0, 1e-30)

def count_identity_changes(indices):
    if len(indices) < 2:
        return 0
    return sum(a != b for a, b in zip(indices, indices[1:]))

def first_arrival(series, threshold):
    for t, val in series:
        if val >= threshold:
            return t
    return None

def sustained_residence(series, threshold, window=0.02):
    """Earliest t such that the diagnostic stays >= threshold for a continuous window."""
    if not series:
        return None
    for i, (t, val) in enumerate(series):
        if val < threshold:
            continue
        # check forward window
        end = t + window
        ok = True
        for t2, v2 in series[i:]:
            if t2 > end:
                break
            if v2 < threshold:
                ok = False
                break
        if ok and any(t2 >= end for t2, _ in series[i:]):
            return t
    return None

def demo_controls():
    """Show the control calculations on synthetic series (no solver required)."""
    # synthetic peak-vorticity series
    control = [(0.05 + i*0.01, 1.7 + 0.01*i) for i in range(20)]
    quench  = [(0.05 + i*0.01, 1.7 - 0.005*i) for i in range(20)]
    thr = 1.5 * 1.7
    result = {
        'first_arrival_control': first_arrival(control, thr),
        'first_arrival_quench': first_arrival(quench, thr),
        'sustained_control': sustained_residence(control, thr),
        'sustained_quench': sustained_residence(quench, thr),
        'energy_budget_example': energy_budget_residual(39.4, 39.5, 0.09),
        'identity_changes_example': count_identity_changes([32, 18, 53, 59, 62]),
        'note': 'These functions are the production control path. Wire to Flow.observe output next.'
    }
    return result

if __name__ == '__main__':
    print(json.dumps(demo_controls(), indent=2))
