# Immutable pilot trajectories

These 64 numerical simulation archives and 64 metadata records are the exact files used for the published study. They contain generated fluid-marker positions, velocities, local gradients and tangent matrices, grouped into 32 independent seeds. They contain no measurements of people or physical locations.

File hashes and independent train/validation/test groups are listed in [verification.json](../results/verification.json). They are kept separate from the fresh-run `data/` cache: source hash checks intentionally reject byte-different source snapshots, including different line endings. Reproduce with the study commands to create a new local cache.
