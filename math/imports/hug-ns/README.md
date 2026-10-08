# Hug flow: code and results

[Project home](../../../README.md) · [All result files](results/) · [Code checks](REVIEW.md) · [Result checks](RESULTS-REVIEW.md)

This experiment prepares a flow with an enclosing gate, then evolves it using an approximation of the unforced, incompressible Navier–Stokes equation in a periodic cube.

**Current status:** the pressure diagnostic is corrected. Twenty-two result files are available. The newest six are reviewed in [the October 8 control update](RESULTS-REVIEW.md#october-8-control-update). The focused checks record **32 passes and three known failures**, listed in [the code review](REVIEW.md#findings-that-need-correction-before-relying-on-larger-runs).

**Report:** [Original two-page PDF](hug-runs.pdf) · [Checks and wording corrections](RESULTS-REVIEW.md#pdf-report). The PDF is preserved as received; the linked review states which conclusions the supplied files support.

## Completed mirror tests

[Four reproduced fluid tests and charts](mirror-fluid/README.md) · [20-run signed calibration](mirror-calibration/README.md). These use the documented mirror-pair starting field, with their own code and measured results. The earlier imported results below retain their original scope.

## Earlier controls

| Received file | Checked finding |
| --- | --- |
| [timestep-half.json](results/timestep-half.json) | Halving dt changes the saved endpoint maximum vorticity by 0.001915% at time 0.16. The producing configuration is only partly documented. |
| [budget-48.json](results/budget-48.json) | Snapshot energy, enstrophy and counts match the supplied 48³ arrays. Integrated dissipation is not included. |
| [ladder-48.json](results/ladder-48.json) | All fixed and relative threshold counts match the supplied arrays. |
| [finer-longer.json](results/finer-longer.json) | 64³ summary reaches maximum vorticity 701.797 at time 0.24; dt and viscosity are omitted. |
| [fd-check.json](results/fd-check.json) | Reports a 32³ finite-difference-advection example. Same-grid method agreement is not established. |
| [reverse-status.json](results/reverse-status.json) | Reports a peak followed by decline and renewed growth through 0.40; the precise starting field is unspecified. |

These are received records, not results from the separately running aligned/reverse/exodus experiment. [Definitions, verification and limits](RESULTS-REVIEW.md#october-8-control-update).

## Earlier results

| File | What is recorded |
| --- | --- |
| [matched-finer.json](results/matched-finer.json) | N=64; eight samples through time 0.12. Maximum vorticity rises from 59.28 to 200.62. |
| [matched-nostop.json](results/matched-nostop.json) | 16 samples through time 0.35. Maximum vorticity starts at 60.20, reaches a saved high of 611.14, and ends at 513.81. |
| [matched-continue.json](results/matched-continue.json) | Nine samples through time 0.1999. The file reports stopping after the maximum passed 400. |
| [matched-strain.json](results/matched-strain.json) | Eight samples through time 0.08. Maximum vorticity rises from 60.20 to 102.46. |

The starting maximum matches across `matched-strain`, `matched-continue` and `matched-nostop`. The N=64 table starts at 59.28 rather than 60.20. Later sample times differ. Their settings and full velocity fields were not included, so the records do not establish a single reproduced trajectory or a comparison between resolutions. [Detailed findings](RESULTS-REVIEW.md#matched-no-stop-attachment).

## What the boundary does

The gate is a function in `solver.py`. Its value is one in the interior and zero near the cube faces. It shapes the starting velocity through the stream function. It is used only at initialization; subsequent steps evolve the velocity and pressure.

The cube has side length 6. Opposite faces connect, so flow leaving one face re-enters through the opposite face. The gate remains in the code. Its written formula has corners in its transition region, and the solver does not measure a gap between two arms.

## Equation and measurements

```text
∂t u + (u · ∇)u = −∇p + ν Δu
∇ · u = 0
external force = 0
```

Here `u` is velocity, `p` is pressure and `ν` is viscosity. The code uses Fourier derivatives, pressure projection and two-stage Heun time stepping. This imported solver has different starting fields from the [earlier ring study](../../notes/hug-boundary.md).

| Term in these files | Meaning |
| --- | --- |
| Peak speed | Largest recorded velocity magnitude on the grid |
| Vorticity, sometimes called “spin” | Curl of velocity: a measure of local rotation |
| Maximum vorticity / `biggest` | Largest reported vorticity magnitude |
| `average` | Matches vorticity RMS at the two independently rebuilt starting states; later attachments do not always specify the normalization |
| `peak_split` | Carrying, pressure and viscosity contributions to the change in local kinetic energy, all evaluated at the fastest grid point |

The pressure correction fixes the reporting function. It does not alter the evolution equation. The older [pressure table](results/pressure-at-peak.json) is preserved as received and needs recalculation before its magnitudes or total can be used.

## Run the code

From this folder:

```sh
python -m pip install -r requirements.txt
python run.py hug --n 32 --time 0.4
python run.py sharper --n 64 --time 2
python run.py tube --n 48 --time 2
python run.py tubes --n 48 --time 2
python run.py compare --n 48 --time 2 --out compare.json
```

`compare` runs the pulled-tube case and records both global term magnitudes and the signed split at the fastest point. The three known numerical issues affect interpretation of longer runs; see [the code review](REVIEW.md).

## What is established

- The code passes the listed projection, short-step, energy and corrected-pressure checks.
- The result files pass checks for finite entries, increasing sample times and recorded arithmetic. The 48³ and 64³ starting fields were also rebuilt.
- The later trajectories have not been independently replayed. Finite saved values establish finite numerical observations; they do not prove smoothness of the continuous equation.

For the mathematical argument, read [the vorticity bound](BOUND.md), [the proposed estimate](step-12.md), and [the short-time calculation](step-12-bound.md), and [the climbing-flow rate comparison](climb-bound.md).

Original upload hashes and reviewed-file hashes are recorded in [provenance.json](provenance.json). [The patch](review-changes.patch) shows every change to the uploaded code and notes, including this wording revision.

