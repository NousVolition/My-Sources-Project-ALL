"""Skeleton: count how often the global-max identity changes among tracked markers.

Reads existing measurements.json style data; does not evolve.
"""
import json
from pathlib import Path

def count_identity_changes(marker_max_indices):
    """Count transitions in the sequence of argmax marker indices."""
    if len(marker_max_indices) < 2:
        return 0
    return sum(a != b for a, b in zip(marker_max_indices, marker_max_indices[1:]))

def analyze_existing(measurements_path):
    data = json.loads(Path(measurements_path).read_text())
    results = {}
    for run, rows in data.get('local_terms', {}).items():
        indices = [r.get('marker_max_index') for r in rows if 'marker_max_index' in r]
        results[run] = {
            'samples': len(indices),
            'identity_changes': count_identity_changes(indices),
            'unique_identities': len(set(indices))
        }
    return results

if __name__ == '__main__':
    # Example: point at a local copy of central-response/measurements.json
    print('Skeleton ready. Supply path to measurements.json to run.')
