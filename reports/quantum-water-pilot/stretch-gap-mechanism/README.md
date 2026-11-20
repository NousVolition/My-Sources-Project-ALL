# Separate diagnostic addition: why the stretch bands are close

This is a new analysis of the saved GFN2-xTB force matrices. It leaves the published band-position result and stretch-gap diagnostic unchanged. No parameters were tuned, no alternative method was run, and no new xTB or Lean calculation was performed. No pass mark is assigned to the gap.

The calculation makes two effects nearly cancel. Movement through the shared oxygen tends to separate the two stretching frequencies. The calculated force coupling between the two O–H bonds brings the frequencies closer together. After including the interaction with bending, the saved H2O calculation gives 3643.31 and 3651.44 cm^-1: a gap of 8.13 cm^-1.

This identifies the source of the small separation within the saved model. It does not identify a particular GFN2-xTB parameter or electronic-energy term as the underlying cause of the discrepancy with the reference.

## What the method does

The published protocol uses unscaled harmonic normal modes of an isolated water molecule, after optimizing its geometry with GFN2-xTB. The program estimates the force response to small atomic displacements. These responses form a force matrix, called the Hessian. Atomic masses enter when its vibration frequencies are calculated. The conventional xTB procedure is documented in the [xTB Hessian guide](https://xtb-docs.readthedocs.io/en/latest/hessian.html).

The study uses its own isotope mass weighting because its xTB build reverted to default masses in the printed frequencies. This addition reads the saved Cartesian Hessians and applies the isotope masses in the original protocol. Both stretch labels were checked from whether the two bonds lengthen together or in opposite directions.

The force coupling is present in the saved matrix; the model does not treat the bonds as two independent springs. For the primary H2O case, the average diagonal bond force constant is 0.4782145792 Eh/bohr^2 and the cross-bond constant is +0.0080904382 Eh/bohr^2. The latter makes the together stretch stiffer than the opposite stretch. Movement of the shared oxygen has the opposing effect on their frequencies. Bending is also included.

## A breakdown of the saved result

The table expresses an algebraic decomposition of the frequency difference in bond-length and bond-angle coordinates. Positive entries widen the gap; negative entries narrow it. These are contributions to one saved result, not outputs from separate modified simulations. This coordinate-based decomposition is not a unique attribution to electronic-energy terms.

| Contribution to the gap, cm^-1 | H2O | D2O |
|---|---:|---:|
| Atomic masses and shared-oxygen geometry | +64.1213 | +88.2146 |
| Cross-bond force coupling | -61.8233 | -45.0781 |
| Coupling to bending, with numerical remainder below 0.00001 | +5.8319 | +11.5310 |
| Total | **8.1299** | **54.6676** |

For H2O, the cross-bond term cancels about 96.4% of the mass-and-geometry stretch term before the bending contribution is included. For D2O, the cancellation is about 51.1%. The saved primary H2O and D2O Cartesian Hessians are identical. Changing the isotope masses changes this balance and explains the different predicted gaps within the same force model.

A second way to see the H2O cancellation: the together-stretch force constant is 3.44% larger than the opposite-stretch force constant, while its kinetic factor is 3.45% smaller. Their products are close. The complete calculation also includes stretch–bend coupling.

## Checks made in this addition

Ten saved GFN2 calculations were reanalyzed: five for each isotope, covering three initial geometries and three Hessian displacement sizes. These are the existing calculations, not ten newly run electronic-structure calculations.

| Isotope | Smallest saved gap | Largest saved gap | Unit |
|---|---:|---:|---|
| H2O | 8.1273568 | 8.1630995 | cm^-1 |
| D2O | 54.6659927 | 54.7029808 | cm^-1 |

The small gap persists through the saved refinements. The new calculation reproduced the saved frequencies within 0.000000075 cm^-1. Reconstructing the vibrations in bond-length and bond-angle coordinates reproduced them within 0.00000045 cm^-1. The coordinate derivatives were also checked numerically, and the reconstructed Cartesian force matrix matched the saved one to a relative norm error below 0.000000001.

Thirty-five source files were checked against the original study's SHA-256 manifest before use and checked for preservation afterward. The existing gap results and CSV were also checked against the separate published gap manifest. The numerical tolerances used here check data reproduction and algebra; they are not acceptance criteria for agreement with physical reference gaps.

## What this adds, and what remains open

The proximate cause is the balance of force coupling, geometry, masses, and bending in the saved GFN2-xTB calculation. The displacement refinements do not explain away the small gap. The reference comparison also cannot be explained solely by comparing harmonic predictions with measured fundamentals: the published like-for-like H2O harmonic reference gap is 3943 - 3832 = 111 cm^-1, which is also far from 8.13. See [NIST CCCBDB's H2O frequencies](https://cccbdb.nist.gov/exp2x.asp?casno=7732185&charge=0).

GFN2-xTB is a broadly parametrized semiempirical method, described by its authors as designed for fast molecular structures and noncovalent-interaction calculations. That scope provides context, not a diagnosis of a specific erroneous parameter. See the [authors' GFN2-xTB paper](https://www.cambridge.org/engage/chemrxiv/article-details/60c73f9e567dfe2444ec3a69).

Identifying which approximation or energy contribution creates the discrepancy would require additional evidence, such as force-curvature comparisons at controlled geometries. No such comparison was run here. Any future comparison against another method must record its gap reference, pass mark, method, and settings before running it; the result must be a separate addition.

The original position result remains a pass under its 5% rule. The original gap diagnostic remains 8.13 against 99.14 for H2O, 54.67 against 116.59 for D2O, and 8.13 against the H2O harmonic reference gap of 111 cm^-1. No prior gap pass mark existed, so those comparisons remain diagnostics.

## Calculation details

Let q = (r1, r2, theta), with lengths in bohr and the angle in radians. B contains derivatives of these coordinates with respect to Cartesian positions. From the saved projected Cartesian Hessian H, this addition constructs

- P = pseudoinverse(B);
- F = P^T H P;
- G = B M^-1 B^T.

The reconstruction H ≈ B^T F B is verified. The eigenvalues of F G reproduce the three internal squared frequencies after unit conversion. This construction describes the saved projected harmonic model at its optimized geometry; it does not fit a new force field.

For equivalent bonds, define k as the average diagonal stretch entry of F, c as its cross-bond entry, g as the average diagonal stretch entry of G, and h as its cross-bond entry. In symmetric and antisymmetric stretch coordinates, the difference between the two diagonal stretch products is

(k - c)(g - h) - (k + c)(g + h) = -2 k h - 2 c g.

The table assigns these two terms to mass/geometry and cross-bond force coupling. It adds the difference contributed by the complete stretch–bend calculation. The small effect of asymmetry and reconstruction is recorded separately in the JSON. To express these squared-frequency terms as additive contributions to the gap in cm^-1, each term is multiplied by the common conversion factor squared and divided by the sum of the two full frequencies. No force term was disabled or adjusted.

## Files and provenance

To reproduce the diagnostic from a checkout with NumPy installed:

```sh
python reports/quantum-water-pilot/stretch-gap-mechanism/stretch_gap_mechanism.py
```

The script reads the existing parent study files and writes the new diagnostic JSON to standard output. `--source` can select another copy of the same source files; their recorded hashes must match.

- [Reproducible analysis](stretch_gap_mechanism.py)
- [New diagnostic results](mechanism-results.json)
- [Preservation and scope checks](verification.json)
- [Addition file hashes](addition-checksums.json)
- [Published method and normal-mode code](https://github.com/NousVolition/My-Sources-Project-ALL/blob/f571de5d8ba4caa1b6393b3954ad7cb04dd4f97a/reports/quantum-water-pilot/normal_modes.py)
- [Saved original run records](https://github.com/NousVolition/My-Sources-Project-ALL/blob/f571de5d8ba4caa1b6393b3954ad7cb04dd4f97a/reports/quantum-water-pilot/vibration-final/runs.json)
- [Primary H2O force matrix](https://github.com/NousVolition/My-Sources-Project-ALL/blob/f571de5d8ba4caa1b6393b3954ad7cb04dd4f97a/reports/quantum-water-pilot/vibration-final/gfn2-H2O-start0-step0.00125/hessian)
- [Unchanged published gap diagnostic](https://github.com/NousVolition/My-Sources-Project-ALL/blob/f571de5d8ba4caa1b6393b3954ad7cb04dd4f97a/reports/quantum-water-pilot/stretch-gap-diagnostic.md)

Source revision: f571de5d8ba4caa1b6393b3954ad7cb04dd4f97a. Original study revision: 3469351c1ef9e0fdd0e121179a06d5b7fb0db639.

Generated by Codex under the user's direction, following the diagnostic request supplied from Claude. This addition reanalyzes the saved source revision listed above. Its numerical checks establish data reproduction and the stated decomposition.
