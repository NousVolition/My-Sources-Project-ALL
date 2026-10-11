# Solver source versions

`dense-jac` preserves the exact source used for the first completed final-matrix
runs. Later runs store the same Jacobian sparsely and group non-conflicting
spatial probes according to a conservative dependency graph. Every boundary
position is allowed to influence every other boundary position for self-contact.
Triangles supply local dependencies. Grouping changes neither forces nor
continuous equations. Exact rest/port/work blocks are retained. Validation
compares the entire colored matrix with the uncolored version.

Each row resolves its recorded hashes to current or archived source bytes.
Completed dense rows are reused only if all three force/model source hashes
match the current files. To reproduce an older solver version, replace the named
source files in a fresh package copy with the archived versions, retaining the
historical baseline data and creating a data directory before validation.

During the storage change, three unfinished calculations were interrupted:
amplitude32-one, long-s1-two and steady-mesh72. Their partial in-memory states
were not saved and are not numerical results or counted runs. Those requested
configurations are rerun to completion with the sparse version. All completed
dense recordings and their hashes remain intact. No result is claimed for the
interrupted calculations.
