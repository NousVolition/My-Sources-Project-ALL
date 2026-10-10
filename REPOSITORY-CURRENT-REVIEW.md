# Current repository review and priorities

October 10, 2026. [Project home](README.md) · [Reports](reports/README.md) · [Historical review](REPOSITORY-REVIEW.md) · [Nous-Volition](https://github.com/NousVolition/Nous-Volition)

This follow-up checks the specific concerns raised in an external review against published files and current study records. It does not independently reproduce every study or establish physical accuracy. The review bases are My-Sources-Project-ALL `6d584e12fe49aa0c53f6201a2dd01ee7c0d9f1a5` and Nous-Volition `8e8a84e783a15b699f5e80542584502f82c24b8f`. Later commits may add results.

## Findings and dispositions

| Concern | Verified scope and next action |
| --- | --- |
| Repository size | The published tree contains about 1.168 GB of file contents in My-Sources-Project-ALL and 118 MB in Nous-Volition. These are sums of current Git blob sizes, not clone sizes, compressed transfers, or timing measurements. Large arrays, generated reports, figures and historical records are present. Avoid routinely committing duplicate bulk outputs; use versioned external data deliveries with hashes and compact summaries when an existing study supports that delivery format. |
| Preserved solver defects | The archived hug-ns source has three deliberate strict expected-failure tests. [corrected_v1](math/imports/hug-ns/corrected_v1/README.md) repairs the gate, cutoff and scheduler, but its four fluid checks only reach t=0.01 on N=16. A longer common-field grid/time study for that version remains unresolved. The existing matched-stretch and vortex matrices have separate frozen sources; assigning all historical defects to all their trajectories is incorrect. |
| Spatial convergence | This concern is supported. The [final vortex matrix](reports/study/final-matrix/README.md) has 5.19–24.21% global-W curve differences from grid 112 to 160, although finest-grid timestep differences are below 0.0014%. The [matched-stretch 128 timestep check](reports/matched-stretch/completed-timestep128/README.md) also fails its stated curve screen. Conservation, reproducibility and small timestep sensitivity do not establish spatial accuracy. |
| Unfinished controls and weak prediction gains | The [published matched-stretch index](reports/matched-stretch/README.md) has 29/48 verified published results; the local queue had 30/48 complete at this review, with 18 remaining. The [vortex matrix](reports/study/final-matrix/README.md) is complete and must not be restarted. The [prediction refinement](reports/matched-stretch/history-prediction/resolution-followup/README.md) and [mirror follow-up](reports/matched-stretch/history-prediction/mirror-followup/README.md) do not establish a robust incremental benefit. |
| Connections between domains | Similar testing methods do not make SIMS, molecular and fluid variables interchangeable. The new [clay/SIMS experiment](https://github.com/NousVolition/Nous-Volition/blob/8e8a84e783a15b699f5e80542584502f82c24b8f/studies/dynamics-fractals-sims/CLAY_SIMS_METHODS.md) does use published wet-kaolin fitted coefficients in a numerical network, with an explicitly chosen reciprocal coupling. That coupling is not calibrated by the clay measurements and does not demonstrate a physical social-to-fluid mechanism. The water-biology coupled-material additions remain in [draft PR #3](https://github.com/NousVolition/My-Sources-Project-ALL/pull/3). |
| No open issues | Both repositories had no open issues in the API listing at this review. My-Sources-Project-ALL had draft PR #3 open. An empty issue list says nothing about scientific completeness; unresolved work is documented here and in each study. |
| Missing source inputs | Exact source equations and inputs remain necessary for reconstructions that lack them. Later separately specified heteroclinic examples do not reconstruct a missing exercise. Preserve source provenance and identify any newly assumed equations. |
| Volume of work | Finish the existing authorized controls and interpret their accuracy before expanding experiments or strengthening claims. No new simulations or matrices are added by this review. |

### Which solver is which?

The original gate corner belongs to the archived gated start. The matched-stretch solver introduces no boundary gate: it preserves a supplied raw strain with nonsmooth periodic joins. Its grids 64, 128 and 256 are not divisible by three, so its inclusive integer cutoff excludes the N/3 endpoint implicated in the archived alias reproducer. Its output intervals use ceiling-based equal subdivisions, rather than the archived rounding-up step-size defect. Its initial-speed timestep choice still requires guards and empirical timestep checks; this distinction is not a stability or accuracy proof. The vortex matrix instead uses a separate smooth periodic start, strict retained cutoff and SSP RK3. Keep those solver and starting-field scopes explicit.

## Cleanup made in this follow-up

- Corrected the offline report landing page, which still displayed 5/48 matched-stretch and 28/48 vortex results. It now shows the current published counts and links the completed vortex report.
- Kept the older matched-stretch preview explicitly historical rather than presenting its image as a current result.
- Added navigation to the newer HUG feedback and prediction checks, with their interpretation limits.
- Repaired confirmed punctuation-encoding damage in two navigation guides. The original dated repository review and its historical numbers remain available.

## Work still to resolve

1. Finish and publish the remaining existing matched-stretch controls through 0.40. Preserve checkpoints, failures, starting field, mean, method, viscosity and zero force. Do not change active numerical source or duplicate simulations.
2. Give physical interpretations only within demonstrated spatial/time accuracy. The nonsmooth raw strain and initial maxima outside the central tube remain material limitations. Finite-time perturbation rates are not asymptotic Lyapunov exponents.
3. Specify a common-initial-field refinement design before any separately authorized long-run corrected_v1 study. Repair regression tests alone are insufficient.
4. Extend automated package verification to newer report families where the root test-discovery paths do not include them. A green repository workflow does not mean every saved experiment was tested by that workflow.
5. Plan bulk-data delivery per study before moving existing files. Preserve stable reading links, exact source versions, immutable hash inventories and a usable reproduction route. Removing files in a new commit does not remove their bytes from Git history. Any history rewrite or bulk storage migration requires a separate reviewed plan.

The supplied q remains a passive signed record; fitted-q and cubic-feedback models are separate completed tests. E is unsigned, q is signed, and mirror averaging is an explicit intervention. The mathematical model remains under examination.
