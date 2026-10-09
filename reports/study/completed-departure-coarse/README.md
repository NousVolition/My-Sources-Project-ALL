# Departure case: first six completed runs

All ordinary-step Fourier and FD4 runs at **48, 64 and 80 grids** reached **t=0.40**. Both 80-grid half-step controls are now [complete and verified](../completed-departure-half80/README.md). Finer grids and their half-step controls remain in the existing authorized matrix. These six results bring verified published vortex records to **38/48**.

![Curves, spectra, widths and budgets](curves.png)

## Completed measurements

| Run | Final W | Final I | High-band enstrophy | Global width (cells) | First spectral warning | First width warning |
| --- | --- | --- | --- | --- | --- | --- |
| exodus-fourier-n48-base | 88.15757 | 22.56190 | 17.385% | 2.014 | 0.18 | 0.18 |
| exodus-fd4-n48-base | 79.14726 | 22.08779 | 17.214% | 1.769 | 0.2 | 0.18 |
| exodus-fourier-n64-base | 102.64848 | 24.27390 | 14.257% | 2.692 | 0.22 | 0.22 |
| exodus-fd4-n64-base | 97.55390 | 23.93819 | 13.842% | 2.821 | 0.22 | 0.22 |
| exodus-fourier-n80-base | 107.31294 | 24.33402 | 13.652% | 2.891 | 0.22 | 0.24 |
| exodus-fd4-n80-base | 103.95238 | 24.13931 | 13.077% | 2.965 | 0.24 | 0.24 |

## Grid and method comparisons

| Comparison | Maximum W-curve difference | Endpoint velocity L2 difference |
| --- | --- | --- |
| exodus-fourier-n48-base to exodus-fd4-n48-base | 11.2629% | 7.4148% |
| exodus-fourier-n48-base to exodus-fourier-n64-base | 19.9854% | 9.3035% |
| exodus-fd4-n48-base to exodus-fd4-n64-base | 23.5508% | 8.7089% |
| exodus-fourier-n64-base to exodus-fd4-n64-base | 5.2223% | 5.6173% |
| exodus-fourier-n64-base to exodus-fourier-n80-base | 4.3466% | 7.5169% |
| exodus-fd4-n64-base to exodus-fd4-n80-base | 6.1552% | 6.1526% |
| exodus-fourier-n80-base to exodus-fd4-n80-base | 4.7104% | 4.2191% |

Each W percentage is the maximum absolute curve difference divided by the second listed run's maximum saved W. Field differences use the common retained Fourier modes; omitted fine-only modes are not included in that error.

The 64-to-80 W-curve differences are **4.3466% for Fourier** and **6.1552% for FD4**, smaller than the corresponding 48-to-64 differences. Nevertheless, all six runs fail spectral and width screens by the endpoint. The endpoint global widths are below three cells, and high-band enstrophy remains above 13%. Spatial convergence is not established.

## Verification and scope

All **30 retained fields** at t=0, 0.1, 0.2, 0.3 and 0.4 were checked across the six runs. Independent physical-space sums reproduce W, energy and both enstrophies. Source/settings provenance, finite arrays, zero mean, saved accumulators and budget gates pass. Frozen diagnostics reproduce endpoint widths, strain and spectra. Recorded stage integrals were retained; no trajectory was rerun.

The prescribed departure starting field, domain side 6, per-run viscosity 0.001, zero force and SSP RK3 are unchanged. The Fourier and FD4 methods share FFT pressure/filter infrastructure and the displayed Fourier-curl W diagnostic. The separate matched-stretch start retains its nonsmooth periodic joins and outside-tube initial maxima.

[Measurements](measurements.json) · [Summary](summary.json) · [Verification](verification.json) · [Comparisons](comparisons.json) · [Protocol](../protocol.json) · [Study index](../README.md)

## Reproduction

Use the existing [runner](../run_study.py), [solver](../solver.py) and [initial construction](../initial_design.py). An example command from the study folder is:

    python run_study.py --case exodus --n 80 --method fourier

Use the corresponding grid and method for each listed control. The runner protects active work with a process lock. Raw fields stay local and are identified by SHA-256. Repeat the saved-field audit with NumPy and SciPy:

    python verify.py /path/to/adversarial-vortex-study

Rebuild the charts from the bundled JSON with NumPy and Matplotlib:

    python plot.py
