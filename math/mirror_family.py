"""Norm / Entropy: named, paired halves of the existing reflection experiment.

The names are supplied by the contributor. They are not physical entropy,
normality, personality assessments, or labels supplied to K-means.
This controlled example uses a center-line slice of the actual smooth field.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from reflection_probe import BOX, FIELD, neighboring_copy


def paired_halves(points=65):
    if points < 3:
        raise ValueError('At least three paired sample locations are required.')
    x = np.linspace(-BOX/2, BOX/2, points)
    norm = np.array([[q, FIELD.velocity(q, 0., 0.)[0]] for q in x])
    entropy = np.array([[q+BOX, neighboring_copy(q+BOX, 0., 0.)[0]] for q in x])
    offset = FIELD.velocity(-BOX/2, 0., 0.)[0] + FIELD.velocity(BOX/2, 0., 0.)[0]
    aligned = entropy.copy()
    aligned[:, 0] -= BOX
    aligned[:, 1] = offset - aligned[:, 1]
    return norm, entropy, aligned


def pair_distances(x):
    return np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)


def fit_groups(features, points, seed):
    with threadpool_limits(limits=1):
        model = KMeans(n_clusters=2, n_init=20, random_state=seed,
                       algorithm='lloyd', max_iter=300).fit(features)
    labels = model.labels_
    counts = [dict(group=int(k), Norm=int(np.count_nonzero(labels[:points] == k)),
                   Entropy=int(np.count_nonzero(labels[points:] == k))) for k in range(2)]
    return dict(seed=seed, groups=counts,
                corresponding_pairs_in_same_group=int(np.count_nonzero(labels[:points] == labels[points:])),
                total_pairs=points,
                inertia=float(model.inertia_))


def run(points=65):
    norm, entropy, aligned = paired_halves(points)
    raw = np.vstack([norm, entropy])
    canonical = np.vstack([norm, aligned])
    # Hold scales fixed: fit once to the placed data and reuse after alignment.
    scaler = StandardScaler().fit(raw)
    result = dict(
        question='Can two differently oriented halves retain the same sampled shape?',
        names={'Norm': 'original slice', 'Entropy': 'neighboring endpoint-matched flipped slice'},
        naming_limit='Names chosen for this experiment, not computed categories or measured entropy.',
        construction='Paired, generated samples. The mirror relationship is constructed, not discovered in independent observations.',
        data='Existing SmoothRingField, y=z=0, x in [-3,3]; the right copy is placed at x+6.',
        features=['position_x', 'velocity_x'], units='Both coordinates are nondimensional.',
        scope='One center-line slice. The 3-D flip has position-dependent offsets away from this slice; global rigid-shape equivalence is not established.',
        points_per_half=points,
        max_pair_error_after_alignment=float(np.max(np.abs(norm-aligned))),
        max_within_half_distance_change=float(np.max(np.abs(pair_distances(norm)-pair_distances(entropy)))),
        kmeans=dict(k=2, reason='Two groups requested as an illustrative fixed setting; not an inferred number of natural families.',
                    features_exclude=['half_name', 'pair_id'],
                    scaling='StandardScaler fitted once to placed data and reused after mirror alignment.',
                    seeds=list(range(5)), n_init=20,
                    placed=[fit_groups(scaler.transform(raw), points, seed) for seed in range(5)],
                    aligned=[fit_groups(scaler.transform(canonical), points, seed) for seed in range(5)]),
        rows=[dict(pair=i, Norm=norm[i].tolist(), Entropy=entropy[i].tolist(), aligned_Entropy=aligned[i].tolist()) for i in range(points)],
        interpretation='One constructed paired family can be displayed as two named halves. K=2 still partitions the values after alignment; that partition does not create two different personalities.',
        boundary_status='Naming or aligning data does not remove the physical join corner or the 3-D divergence reported by reflection_probe.py.',
    )
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    result = run()
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({key: result[key] for key in ('points_per_half', 'max_pair_error_after_alignment', 'max_within_half_distance_change', 'kmeans')}, indent=2))
