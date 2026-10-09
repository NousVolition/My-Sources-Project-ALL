# Departure case: both 160-grid ordinary-step runs complete

The **Fourier and FD4 ordinary-step runs reached model time 0.40**. The vortex study now has **46/48 completed results verified and published**. Both final 160-grid half-step controls are running.

## Global and central peaks respond differently to refinement

| Method | Global W, 80 to 112 | Global W, 112 to 160 | Central W, 112 to 160 |
| --- | ---: | ---: | ---: |
| fourier | 19.1009% | 24.2055% | 1.3768% |
| fd4 | 16.6786% | 22.5728% | 0.5090% |

Global-peak curve differences increase under the latest refinement for both methods. Spatial convergence remains unestablished. Central-region peak differences are smaller, but each percentage uses that diagnostic's own maximum on the finer grid as its normalization. These values describe different regions and different scales.

![Central and global peak comparison](central-global.png)

At t=0.40, the global peaks are at **(0, -3, 0.1875)** for Fourier and **(0, -3, 0.15)** for FD4, on a periodic face. Their widths are approximately **1.48 and 1.56 cells**, compared with central core widths of **9.33 and 9.57 cells**. The wider central core does not resolve the narrower whole-domain maximum.

## Endpoint measurements and budgets

| Measurement at t = 0.40 | Fourier 160 | FD4 160 |
| --- | ---: | ---: |
| Global W | 168.15608 | 159.51047 |
| Central-region W | 78.46496 | 76.766362 |
| I | 24.855359 | 24.361568 |
| High-band enstrophy (%) | 8.0987532 | 7.552198 |
| Global peak width (cells) | 1.481137 | 1.562455 |
| Central core width (cells) | 9.3266664 | 9.5747383 |
| Relative energy-budget residual | -1.1453806e-07 | -1.17808e-07 |
| Relative native-enstrophy-budget residual | -4.603057e-05 | -3.6036759e-05 |

![Curves, spectra, widths and budgets](curves.png)

Both runs fail the 1% high-band-enstrophy and six-cell global-width screens. Their maximum absolute relative energy-budget residuals are **1.14538e-07** and **1.17808e-07**, respectively. Budget consistency does not establish spatial resolution.

At grid 160 the maximum Fourier/FD4 W-curve difference is **8.1501%**, normalized by the FD4 saved peak. The corresponding central-region W difference is **2.1232%**. These are comparisons of the implemented discretizations; they do not establish continuum accuracy.

The FD4 projection enforces its native discrete divergence. At the endpoint its recorded maximum native divergence is **1.29415e-15**, while the Fourier-derivative divergence diagnostic is **3.84113**. W and the spectral tail here use the Fourier curl for both methods; the FD4 enstrophy budget uses its native curl. The two enstrophy definitions are retained separately in the records.

The existing reporter's 112-to-160 common-mode velocity comparisons give endpoint L2 differences of **3.0671%** and **2.7907%**. Those comparisons restrict both fields to coarser retained Fourier modes and exclude fine-only modes; they are not full fine-grid errors. This batch independently recalculated **21 scalar curve comparisons**, while retaining the reporter's field-comparison values.

## Verification and preserved model

All **five newly completed FD4 fields**, at t=0, 0.1, 0.2, 0.3 and 0.4, passed independent physical-space summation checks of W, energy and native/Fourier enstrophy. Source/settings hashes, finite arrays, mean, CFL and energy-budget gates pass. Frozen endpoint diagnostics reproduce widths, strain, central-region W and spectra. Previous audits of the five Fourier 160 fields and ten 112-grid reference fields were reused after confirming unchanged result and source hashes. No trajectory was rerun; I and budget integrals retain recorded RK-stage accumulations.

The prescribed departure start, domain side 6, viscosity 0.001, zero mean, zero force, SSP RK3 and ordinary step 0.0004 remain unchanged. This starting field differs from the matched-stretch field, whose raw strain is not smooth across periodic joins and whose initial maxima lie outside its central tube. A peak on a periodic face in this separate study does not establish the same initial-field defect. The mathematical model remains under examination.

[Measurements](measurements.json) · [Summary](summary.json) · [Verification and field hashes](verification.json) · [Comparisons](comparisons.json) · [Earlier 112-grid audit](../completed-departure112-base/verification.json) · [112-grid timestep controls](../completed-departure-112-half/README.md) · [Study index](../README.md)

## Reproduction

Use the preserved [runner](../run_study.py), [solver](../solver.py), [initial construction](../initial_design.py) and [protocol](../protocol.json). Raw fields remain local; hashes are recorded. With NumPy and SciPy, run `python verify.py /path/to/adversarial-vortex-study` to audit the two completed 160-grid ordinary runs. With NumPy and Matplotlib, run `python plot.py` to rebuild both figures from the bundled JSON.
