# Consolidation record

[← Archive](README.md)

Reviewed 5 October 2026 across both repositories and all 11 branches. The audit covered 181 distinct repository/path combinations and 221 blob versions. This is the record of how saved work was reused.

## Copilot branches

| Saved branch | Commit at review | Resolution |
| --- | --- | --- |
| copilot/add-divergence-and-project-methods | [2fd83ea](https://github.com/NousVolition/My-Sources-Project-ALL/tree/2fd83ea3f24199417b71954af191b39958651c28) | Distinct tests retained in test_navier.py. Projection corrected for dx squared and convergence; the old return-None assertion now checks the returned residual. |
| copilot/completed-runs-in-project | [af8d93d](https://github.com/NousVolition/My-Sources-Project-ALL/tree/af8d93dbff392eebc9b8dc4031892fac58bd2438) | Completed archive of 18 historical runs remains at its existing commit. Link to it; do not copy its 104 files into active code. |
| copilot/create-dummy-dataset | [aaa52ad](https://github.com/NousVolition/My-Sources-Project-ALL/tree/aaa52ad542822eadb8295327818211580cfc09c3) | Same tree as original main. Retrieved sessions contain a housing-regression example and installation guidance; zero saved file changes and no artifacts. |
| copilot/fix-pull-request-requests | [aaa52ad](https://github.com/NousVolition/My-Sources-Project-ALL/tree/aaa52ad542822eadb8295327818211580cfc09c3) | Same tree as original main. Retrieved session is an explanation of pull requests; zero saved file changes and no artifacts. |
| copilot/initial-axisymmetric-field | [ae551b8](https://github.com/NousVolition/My-Sources-Project-ALL/tree/ae551b87f6c141c6de6b3717eb9c460337f8e8c4) | Original Gaussian consolidated in math/stokes.py; six distinct alpha=0.7 checks retained in test_legacy_parameters.py. |
| copilot/put-stokes-start-on-cube | [3e699d4](https://github.com/NousVolition/My-Sources-Project-ALL/tree/3e699d445d38893d8b3fb908fc3caad9e7e89c22) | Existing PR #1 is the integration point. One sampler and diagnostic runner reused. Tests retain the raw start metrics; the old nonzero residual fixture is replaced with the corrected projection tolerance. |
| copilot/stokes-start-box-length-18 | [d430acf](https://github.com/NousVolition/My-Sources-Project-ALL/tree/d430acf20cb4e6eaef99d266b51968cde815a0d0) | Length-18 runs are options in box_experiment.py. Axis, non-unit spacing, wrap, and evolution checks retained in test_box_integration.py; repeated uniform check is in test_navier.py. |
| copilot/test-meridional-velocity-functions | [9ff7a76](https://github.com/NousVolition/My-Sources-Project-ALL/tree/9ff7a766a17b5be5862ecdc8d97ed61134410531) | Seven geometry checks retained in test_stokes.py; duplicate streamlines implementation retired. |
| copilot/test-velocity-functions | [4a8fdf4](https://github.com/NousVolition/My-Sources-Project-ALL/tree/4a8fdf447e6de270f2d076f00db9bc322c27ae3f) | Seven submitted checks retained in math/tests/test_stokes.py; duplicate stokes implementation retired. |
| main | [aaa52ad](https://github.com/NousVolition/My-Sources-Project-ALL/tree/aaa52ad542822eadb8295327818211580cfc09c3) | Original cube update, regime/survival logic, reports, papers, media, and earlier source retained under the new layout. |

## One implementation per role

- `math/initial_field.py` is the maintained smooth-field implementation. The earlier local companion copy has been consolidated here.
- `math/stokes.py` preserves the original radial Gaussian for comparisons. It is a different profile with an axis cusp.
- `math/navier.py` retains the original explicit update and supplies a compatible projection. The three older pressure variants remain visible in history.
- `math/box_experiment.py` handles both profiles, both box lengths, grid sizes, and time-step options.
- Tests live alongside their subject; one root pytest configuration collects them all.

## Corrections recorded once

For centered divergence and centered gradient, the composed Poisson operator uses neighbors two cells away with the factor 4 dx². The maintained damped Jacobi and Fourier backends solve that same equation. The original nearest-neighbor variants did not. Projection can change gradients; the earlier PR description claiming unchanged gradients was incorrect.

The original Gaussian is divergence-free for positive radius but its Cartesian axial component has a cusp at the axis. The smooth replacement and energy proof were already completed during this review and are reused here.

## What the GitHub session records establish

PR #1 had an approval on 4 October for its earlier version; that is not an approval of later corrections. Retrieved logs show one Copilot failure due to an unavailable model, not a failing math assertion. Two branch families produced guidance with no committed files. Full private Copilot reply text was not exposed by the retrieved logs.

## Preserved material

Papers, original images, CSV results, text examples, and saved raw logs retain their original bytes. HTML changes repair links and label the missing connection-map chart. Generated bytecode and redundant nested test configuration are omitted. No earlier GitHub branch was deleted.

The old report remains in [Nous-Volition](https://github.com/NousVolition/Nous-Volition). Current code and project status have one maintained home here.
