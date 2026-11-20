# The gap between the stretch bands

The band positions pass the 5% rule; the gap between the two stretch bands is **8 against 99 for H₂O** and **55 against 117 for D₂O**, in cm⁻¹. No pass mark was set for the gap beforehand, so this is a **diagnostic**, not a failed test. The original position result is unchanged.

The check subtracts the symmetric stretch from the antisymmetric stretch in the [saved vibration table](vibration_comparison.csv). Using its unrounded numbers:

| Molecule | Predicted gap (cm⁻¹) | Measured infrared gap (cm⁻¹) | Gap underestimated by |
| --- | ---: | ---: | ---: |
| H₂O | 8.1299390581431 | 99.14 | 91.7995% |
| D₂O | 54.6675848803566 | 116.59 | 53.1113% |

The position rule measures each band's error relative to that band's position. It does not measure how accurately their separation is reproduced.

## H₂O compared with the same kind of reference

The model's harmonic stretch frequencies also have a like-for-like comparison with the [NIST harmonic reference](https://cccbdb.nist.gov/exp2x.asp?casno=7732185&charge=0). That reference gives a separation of **3943 − 3832 = 111 cm⁻¹**, against the predicted **8.1299390581431 cm⁻¹**. This remains a gap diagnostic with no pass mark.

The measured infrared gaps above and this harmonic reference gap are separate comparisons. The original [position summary](vibration-summary.json), [protocol](protocol.json), and [148-check verification record](verification.json) remain unchanged.

## Run the check

From the repository root:

```sh
python reports/quantum-water-pilot/stretch_gap_check.py
```

The script prints the [saved diagnostic results](stretch-gap-results.json), including the H₂O harmonic comparison. An optional CSV argument selects the input; the output records only its file name and SHA-256, never its local folder path. The [diagnostic checksums](stretch-gap-checksums.json) bind this script, its results, this explanation and the unchanged input table.

This calculation uses saved predictions. It adds no molecular or fluid simulation and no Lean result.

Review changes supplied by NousVolition from Claude; implementation and arithmetic checks generated with Codex under NousVolition's direction.
