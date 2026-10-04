# Workflow run logs

Downloaded logs for every **completed** GitHub Actions run in this repository
(as of 2026-10-04). Each folder is named by its run ID and contains the raw log
files exactly as produced by GitHub (`<job>/system.txt`, numbered step logs, and
a combined `0_<job>.txt`).

| Run ID | Workflow | Branch | Conclusion | Started (UTC) |
|---|---|---|---|---|
| 37176739266 | Python application | main | failure | 2026-10-04 04:20 |
| 37177590823 | Python application | main | failure | 2026-10-04 04:38 |
| 37178003892 | Python application | main | failure | 2026-10-04 04:47 |
| 37178092719 | Python application | main | failure | 2026-10-04 04:49 |
| 37178460865 | Python application | main | failure | 2026-10-04 04:56 |
| 37185031535 | Python application | main | failure | 2026-10-04 07:11 |
| 37185160840 | Python application | main | failure | 2026-10-04 07:14 |
| 37185160887 | CodeQL Advanced | main | success | 2026-10-04 07:14 |
| 37186324746 | CodeQL Advanced | main | success | 2026-10-04 07:36 |
| 37186324810 | Python application | main | success | 2026-10-04 07:36 |
| 37187754657 | Python application | main | success | 2026-10-04 08:04 |
| 37187754680 | CodeQL Advanced | main | success | 2026-10-04 08:04 |
| 37189983831 | Copilot cloud agent | copilot/initial-axisymmetric-field | success | 2026-10-04 08:46 |
| 37189994854 | Copilot cloud agent | copilot/test-velocity-functions | success | 2026-10-04 08:46 |
| 37190089353 | Copilot cloud agent | copilot/test-meridional-velocity-functions | success | 2026-10-04 08:48 |
| 37191352111 | Copilot cloud agent | copilot/add-divergence-and-project-methods | success | 2026-10-04 09:11 |
| 37191399986 | Copilot cloud agent | copilot/put-stokes-start-on-cube | success | 2026-10-04 09:12 |
| 37191548811 | Copilot cloud agent | copilot/test-velocity-functions | success | 2026-10-04 09:15 |

Note: the two `action_required` runs on `copilot/put-stokes-start-on-cube`
(37191571518 and 37191571496) have no downloadable logs (GitHub returns 404
for them), so they are not included.
