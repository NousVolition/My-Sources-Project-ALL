QUANTUM WATER: COMPLETED RESEARCH PILOT AND STRESS TESTS

Start with report.html. It is self-contained and opens offline in a browser.
It distinguishes passed implementation checks from failed research screens.
The simulations establish a pilot vibrational prediction and a measured-property
flow connection. They do not establish a validated nuclear-quantum explanation
of bulk viscosity or an end-to-end quantum-to-fluid prediction.

CONTENTS
  report.html                    Illustrated scientific report and limitations
  *_comparison.csv, *_results.csv Numerical tables
  *.png, *.svg                   Exportable scientific figures
  vibration-final/               Actual xTB logs, Hessians and mode calculations
  molecular-data/                All molecular trajectories, configs and summaries
  fluid-data/, fluid-stress/, fluid-refinement/
                                Continuum results, Fourier fields and controls
  protocol*.json, *protocol.json Advance criteria and designs
  amendment.json                Numerical repair and isotope-output failure
  followup-amendments.json       Additional checks after initial findings
  provenance.json, sources.json Versions and primary literature references
  verification.json             Execution-integrity checks; not hypothesis proof
  checksums.json                 SHA-256 of the delivered files

QUICK REPRODUCTION OF ANALYSIS (Windows PowerShell, inside this directory)
Use Python 3.12. The delivered data do not require xTB to regenerate charts.

  py -3.12 -m venv .venv
  .\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
  .\.venv\Scripts\python.exe verify_results.py --data .

Verification reads the existing results and checks their recorded hashes,
recomputes Hessian frequencies and material properties, checks saved trajectory
shapes and compares saved channel fields with the analytic solution.

To regenerate the analysis and report (this changes derived file hashes):
  .\.venv\Scripts\python.exe analyze_molecular.py --data molecular-data
  .\.venv\Scripts\python.exe build_report.py --data .

To validate a regenerated copy and record its new hashes:
  .\.venv\Scripts\python.exe verify_results.py --data . --write-manifest
Keep the original archive if you need the original file hashes for comparison.

FULL REPRODUCTION
Download the xTB executable with the supplied installer. This writes only to a
new tools subdirectory, verifies the archive hash and preserves the license files.

  powershell -File install_xtb.ps1
  .\.venv\Scripts\python.exe run_research.py --xtb tools\xtb-6.7.1\bin\xtb.exe --out rerun

The default molecular platform is OpenCL. The original run used an RTX 3080.
For a different supported OpenMM platform, add --platform CPU or --platform
Reference. CPU support depends on your OpenMM build; Reference is much slower.
Expect hours if the GPU is shared, and substantially longer on CPU.
There are 24 molecular cases, 12 primary vibration executions, 84 three-dimensional
vortex runs, 36 channel cases, 10 channel convergence runs, an exact vortex check
and four harmonic-oscillator sampling controls. HDO reuses the primary Hessian.
No cloud compute, paid service, publishing or upstream repository write is used.

Every simulation output must use a new directory. run_research.py refuses to
overwrite an existing output root. The isolated rerun folder contains new results;
the executable source remains here. Run-level seeds, parameters and serialized
OpenMM systems are included. Stochastic GPU trajectories can differ across
hardware or library builds; numerical/statistical agreement, not identical bits,
is the meaningful reproduction target.

INDIVIDUAL COMMANDS
  python run_vibrations.py --xtb <absolute-path-to-xtb.exe> --out new-vibrations
  python check_physics.py --data . --out new-physics-checks.json
  python run_molecular.py --validate --out new-force-check.json
  python run_molecular_ensemble.py --out new-molecules --stage main
  python run_molecular_ensemble.py --out new-molecules --stage controls
  python run_molecular_ensemble.py --out new-molecules --stage temperature
  python analyze_molecular.py --data new-molecules
  python run_fluid.py --out new-fluids
  python extend_fluid_checks.py --data new-fluids
  python channel_error_budget.py --data <complete-result-root>
  python stress_more.py --data <complete-result-root> --fluid
  python refine_fast_flow.py --data <complete-result-root>

The report builder expects the directory structure produced by run_research.py.
Protocols are retained as originally written plus explicit amendments, rather
than retrospectively rewriting the advance criteria to match results.

NUMERICAL INTERPRETATION
- A classical oscillator on the same Hessian predicts the harmonic isotope
  frequency shift too. It is not unique proof that quantum nuclei are essential.
- q-TIP4P/F is unchanged between isotopes and nuclear treatments. Classical
  equilibrium positions cannot depend on masses at fixed potential and state;
  short finite trajectories can show apparent structural differences.
- Nuclear quantum sampling is PIMD equilibrium; RPMD dynamics is approximate.
- Three independent main seeds are the statistical replicas, not beads or frames.
- The liquid pilot is too short to validate bulk diffusion or viscosity.
- Experimental neutron distances are not mean bond lengths. The reported radial
  peak comparison is approximate because no full scattering forward model is used.
- IAPWS inputs are measurement-based correlations. They are not predicted by xTB
  or this molecular simulation. Historical viscosity comparisons are not holdouts.
- Continuum energy conservation is necessary but does not prove grid convergence.
- At 100 mm/s, the original 32-to-48 grid comparison fails. A bounded follow-up
  checks 64 and 96 points per axis, plus a half time step at 64. The report retains
  both the original failure and the follow-up outcome. Agreement between grids
  does not by itself bound error against the unknown exact solution.
- A small isotope-difference field needs its own resolution check; accurate
  total flow does not guarantee the same relative accuracy for a smaller effect.

PROVENANCE AND DISTRIBUTION
The unchanged Navier-Stokes solver comes from the named user's project at the
commit recorded in provenance.json. The other experiment code was authored for
this study. Third-party libraries and xTB retain their respective licenses.
The package includes source URLs and short transcribed numerical reference data;
it does not redistribute complete journal papers or the xTB binary.

This repository copy is stored under reports/quantum-water-pilot. The original
local research archive is preserved separately. README.md provides a GitHub
reading guide; all raw results and original scientific conclusions are retained.
These results do not establish validated quantum transport.
