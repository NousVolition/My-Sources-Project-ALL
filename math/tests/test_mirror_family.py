"""Check the controlled mirror relationship without assuming cluster identities."""
import numpy as np

from mirror_family import paired_halves, pair_distances, fit_groups


def test_actual_slice_matches_after_inverse_reflection():
    norm, entropy, aligned = paired_halves(65)
    np.testing.assert_allclose(norm, aligned, rtol=0, atol=2e-14)
    # Preserves geometry of this 2-D slice under a translation and reflection.
    np.testing.assert_allclose(pair_distances(norm), pair_distances(entropy), rtol=0, atol=2e-14)


def test_changing_one_copy_breaks_the_claimed_pairing():
    norm, entropy, aligned = paired_halves(65)
    aligned[12, 1] += .1
    assert np.max(np.abs(norm-aligned)) > .09
    entropy[12, 1] += .1
    assert np.max(np.abs(pair_distances(norm)-pair_distances(entropy))) > .01


def test_kmeans_assigns_coincident_pairs_together():
    norm, _, aligned = paired_halves(33)
    rows = np.vstack([norm, aligned])
    for seed in (0, 1):
        result = fit_groups(rows, len(norm), seed)
        assert result['corresponding_pairs_in_same_group'] == len(norm)
        assert all(group['Norm'] == group['Entropy'] for group in result['groups'])
