# Current repository review and connections

Reviewed October 11, 2026 (UTC). [Project home](README.md) · [Reports](reports/README.md) · [Status](STATUS.md#current-numerical-status) · [Historical review](REPOSITORY-REVIEW.md) · [Streams and Rocks](https://github.com/NousVolition/Nous-Volition)

## Storage follow-up

The subsequent [verified release migration](DATA-RELEASES.md) moves the numerical arrays and archive parts off main while preserving original bytes and complete offline snapshots. The inventory figures and disposition discussion below describe the earlier pinned audit trees. Use the release guide for current download and restore instructions.

## What changed in this review

- Corrected the current review's obsolete description of PR #3 as a draft. It was merged on October 11; its material-to-fluid feedback study is published and the merge checks passed.
- Replaced old links that presented the historical review as current, corrected two damaged punctuation strings, and removed the outdated 29/48 count from the offline report landing page. That page now directs readers to the maintained study index.
- Updated the reading routes for coupled material, fresh-flow memory, and the completed full-history refinement. The original dated reviews, failed controls and source records remain available.

The audit covers the complete published trees at My-Sources-Project-ALL [`90ce758`](https://github.com/NousVolition/My-Sources-Project-ALL/tree/90ce7583ddd28eb1ed08b939345443ffe82b7348) and Nous-Volition [`1811f92`](https://github.com/NousVolition/Nous-Volition/tree/1811f9273e572d0fc3685a4dcd9288393eda3ca0). Later publications can change counts. [Machine-readable inventory and all identical-file groups](reports/repository-audit/2026-10-11.json).

## What can be reduced, and what should be kept

| Published tree | Files | File bytes, decimal | Main storage use |
| --- | ---: | ---: | --- |
| My-Sources-Project-ALL | 6,330 | 2.024 GB | Compressed numerical arrays: 1.496 GB; matched-stretch and its separate follow-ups together: 1.045 GB |
| Nous-Volition | 531 | 138.55 MB | Social-organization, recovery and dynamics study packages |

These are sums of file contents at the two commits, not compressed clone size, transfer size, Git history size or measured browsing time. The matched-stretch directory total includes separate prediction studies; it is not the size of the original 48-job matrix alone.

There are **349 nonempty identical-blob groups, with 584 additional file copies totaling 18.99 MB** across both trees. This is under 0.9% of their combined file bytes. Git already shares identical blob objects within a repository, so deleting these paths would not provide the same saving in Git storage. This audit detects byte-identical files, not all repeated numerical content inside differently packaged archives.

| Example | Interpretation and disposition |
| --- | --- |
| `clay-hug/data/contact-100.npz`, `forming-clay.npz`, `mesh-24.npz` | One identical baseline appears in three named comparisons. The runner uses that same parameter setting in contact, forming and mesh controls. Keep the comparison identities; a future package format could reference one canonical array. |
| Nine identical balanced-start files in `hug-dynamics/data/sweep-*.npz` for a given load | The runner explicitly tests that exactly balanced starts stay balanced across stiffness/damping settings. Matching outcomes are a null control, not evidence that distinct runs were accidentally copied. |
| Reversal `recorded-data` and `recovery` branch-state copies | A possible delivery consolidation, but recovery provenance and references must be preserved. Do not delete a restart record based only on its matching hash. |
| Fluid-organization figures and tables shared with Nous-Volition | Deliberate reader-facing mirrors. Runnable primary fluid code/data live here; the companion repository preserves its own reading path. |
| Numbered archive parts, failed attempts, dated reviews and original solver source | Keep. Archive parts are complementary pieces, and historical records explain corrections and evidence limits. |

The useful storage change is **study-specific bulk-data delivery**, following the existing [fresh-flow raw-recordings release](https://github.com/NousVolition/My-Sources-Project-ALL/releases/tag/fluid-memory-erasure-2026-10-10): retain compact reports, code, checksums and exact restore instructions on main, with verified bulk recordings in versioned releases. Existing arrays should move only after download/restore verification, reference repair and reproduction checks. Removing them in a new commit would not shrink old Git history. No bulk deletion or history rewrite was made by this audit. Plain `.log` files total only 0.23 MB in the main tree; deleting those would have negligible storage benefit.

## Connections worth seeing

| Connection | What the existing evidence supports | Limit or next unresolved question |
| --- | --- | --- |
| **Material organization → fluid response** | The merged [coupled material test](reports/water-biology-study/coupled-feedback/README.md) transports concentration, makes viscosity depend on it, and feeds viscous stress back into velocity. All 30 pairs / 60 trajectories reached t=1. Feedback-off and uniform controls give zero organization contrast. | The viscosity law is assumed. A frozen heterogeneous-viscosity control is still needed to isolate continued transport from the initial viscosity pattern; biological calibration is absent. |
| **History helps prediction → history changes response** | [Fresh-flow memory tests](reports/fluid-memory-erasure/README.md) report a 7.36% improvement for ordinary forward deformation, but a 0.77% worse error for probe-induced change. Scrambling attenuates the prepared response contrast. | Those are different outcomes. Extra attenuation beyond a strength-matched intact displacement is not established; the positive erasure result does not establish a predictive-response benefit or an extra fluid-memory force. |
| **Better numerical resolution → stronger evidence** | The [fixed-state refinement](reports/matched-stretch/history-prediction/spatial-convergence/README.md) passes its primary short-interval screens. The later [full-history review](reports/matched-stretch/history-prediction/full-history-refinement/FINAL-REVIEW.md) passes 64/96 screens at grid 224; all 32 remaining failures concern individual predictor-output drift. | These are different scopes. Stable short-interval targets do not establish spatial convergence of the entire fitted predictor or the original sharp-field matrix. History-error gains remain small and uncertain in the reused-seed full-history experiment. |
| **Clay → retained shape → a measured HUG** | The [wet-kaolin Burgers model](https://github.com/NousVolition/Nous-Volition/blob/main/studies/dynamics-fractals-sims/CLAY_SIMS_METHODS.md) separates immediate deformation, delayed recovery and retained strain. Published fitted coefficients provide a material starting point. The [calibration notebook](https://github.com/NousVolition/Nous-Volition/blob/main/studies/dynamics-fractals-sims/HUG_CALIBRATION_METHODS.md) provides a route to specimen-specific measurements. | Network coupling is chosen; the notebook's completed examples are synthetic. Physical recordings are pending. The separate [2-D molding mesh](reports/clay-hug/README.md) is dimensionless and is not that measured kaolin model. |
| **Clay memory → sustained rhythm** | The [active HUG extension](https://github.com/NousVolition/Nous-Volition/blob/main/studies/dynamics-fractals-sims/ACTIVE_HUG_METHODS.md) adds inertia, restoring stiffness and powered feedback to a clay internal law, with explicit cycle and energy checks. | Sustained active motion relies on that added power source and selected mechanics. It is not established as a spontaneous rhythm of passive clay or a physically calibrated specimen. |
| **Names, positions, mirrors and response** | [SIMS](https://github.com/NousVolition/Nous-Volition/blob/main/studies/social-organization/README.md), [molecular placement](reports/molecular/isotope-mass/README.md), and [fluid marker controls](reports/matched-stretch/history-prediction/spatial-convergence/README.md) distinguish relabeling from physically changing a position, geometry or field. | This is a shared experimental method. Their variables and coefficients are not interchangeable, and no quantitative social-to-molecular-to-fluid coupling is established. |

## Current work and priorities

1. **Finish the existing matched-stretch matrix.** At this audit, [34/48 completed results are verified and published through 0.40](reports/matched-stretch/README.md); 14 remain. Preserve the supplied field, Heun method, domain, mean, per-run viscosity, zero force, failures and checkpoints. The [separate vortex matrix](reports/study/final-matrix/README.md) is complete at 48/48 and must not be repeated.
2. **Keep numerical limits prominent.** The original raw strain is nonsmooth across periodic joins, initial maxima are outside the central tube, and the [128-grid timestep screen fails](reports/matched-stretch/completed-timestep128/README.md). The final vortex matrix retains 5.19–24.21% finest-grid global-W curve differences despite tiny timestep differences. Conservation and reproducibility do not establish spatial accuracy. Finite-time perturbation rates are not asymptotic Lyapunov exponents.
3. **Follow existing controls before interpreting a new memory mechanism.** The separately managed [equal-work/location/independent-solver extension](reports/fluid-memory-extension/README.md) is an in-progress publication, with preserved failed pilots and partial corrected qualification. Its pending production is not a completed result. This audit adds no simulations and does not take over that work.
4. **Keep solver versions distinct.** The archived hug-ns has three documented expected failures. [corrected_v1](math/imports/hug-ns/corrected_v1/README.md) repairs them, but its four fluid checks only reach t=0.01 at N=16; a separately specified long common-field refinement remains unresolved. The matched-stretch solver has no boundary gate, uses grids excluding the implicated N/3 endpoint, and subdivides output intervals with a ceiling. The vortex matrix has another starting field and method. None of these distinctions is an accuracy proof.
5. **Use physical measurements for physical claims.** Preserve missing-equation/source qualifications, and distinguish synthetic calibration, published fitted coefficients and measurements of a particular specimen. New controls or bulk migration need a concrete plan; this review does not expand any matrix.

## What was checked

Both full Git tree inventories were inspected without truncation. All **178 published Markdown files** were read from bytes matched to their pinned Git blob hashes. Inline local Markdown paths were checked against those trees. Four unavailable links are inside the already documented, hash-frozen historical `obstacle-only/MODEL-NOTES.md`; their source is preserved and the current [contact report](reports/clay-contact-stress/README.md) remains the reading route. An equation resembling Markdown was excluded as a parser false positive. No other missing inline local Markdown file targets were found. This is not a complete audit of external URLs, all HTML pages or every heading fragment.

At the audited heads, [main Project checks](https://github.com/NousVolition/My-Sources-Project-ALL/actions/runs/38171455384), [CodeQL](https://github.com/NousVolition/My-Sources-Project-ALL/actions/runs/38171455555), and [Nous-Volition checks](https://github.com/NousVolition/Nous-Volition/actions/runs/38171926476) passed. Both open-issue listings, which also include open PRs, were empty. [PR #3](https://github.com/NousVolition/My-Sources-Project-ALL/pull/3) is merged, and its [post-merge Project checks](https://github.com/NousVolition/My-Sources-Project-ALL/actions/runs/38106222184) passed. A green workflow covers its configured checks, not every scientific claim or saved experiment. This review does not independently rerun the numerical ensembles.

The supplied q remains a passive signed record. Fitted-q and cubic-feedback tests are separate models, E is unsigned, and mirror averaging is an explicit intervention. The mathematical model remains under examination.
