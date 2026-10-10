# Pressure, integral and ratio: current status

Updated October 10, 2026. [Current overview](README.md) · [Historical review](archive/earlier-import/RESULTS-REVIEW.md) · [Calculation record](diagnostic-audit.json) · [Recheck script](diagnostic_audit.py)

## Already repaired in code

The pressure sign and filtered carrying term are corrected. The [separate corrected version](corrected_v1/README.md) also repairs the gate, retained cutoff alias and step scheduling. The four recorded numerical source hashes still match that version's saved verification. No active numerical source or existing simulation was changed by this audit, and no fluid trajectory was rerun.

## Historical pressure values need the original fields

The received [pressure table](results/pressure-at-peak.json) is preserved unchanged. Reversing its pressure column cannot recover the solver-consistent carrying term. Supply either the original velocity arrays at each pressure observation or a complete reproducible original run, including viscosity, domain/grid, solver/filter revision and exact time-step schedule. The corrected diagnostic can then be evaluated on those states. A newly chosen start would be a separate experiment.

## The vorticity integral discrepancy is quantified

The [received table](results/bound-check.json) reports **86.12855895**. Trapezoids using its 12 saved rows give **88.67727603**, a difference of **2.54871708** (2.95920% of the reported value). Every interval contribution is included in the calculation record.

If each printed time was rounded to the nearest four decimal places and the saved vorticity values are held fixed, time rounding can change this trapezoidal result by at most **0.00685941**. Time rounding alone cannot account for the discrepancy under those assumptions. More frequent sampling or a different quadrature rule remains a possible explanation; neither has been recovered for this run.

To recover the reported integral, provide the quadrature formula and every unrounded time and vorticity maximum used in it, or its complete reproducible producer. The saved-row trapezoid is separately labeled and does not replace the historical value. Neither number proves regularity of the continuous equation; [BOUND.md](BOUND.md) explains that requirement.

## The ratio observation is arithmetically consistent

At the maximum saved ratio, **32.20099629 at t=0.1977**, the numerator is **241.45869042**, down **0.104909%** from the start; the recorded denominator is **7.49848508**, down **5.003888%**. The ratio rises **5.157031%**. The maximum saved vorticity, **250.81187785**, occurs later at **t=0.3163**.

This verifies the relationship between the supplied numbers. Reproduction still needs the grid, domain, viscosity, exact initial field or construction, force, solver/filter revision, step schedule and denominator's norm definition. The table's word `average` has not been silently redefined as RMS.

## Search for the missing records

The local search covered **36 supplied hug archives**, including later ZIP and ZIP-formatted AAF deliveries. Each of these three table names had only one distinct content hash in those archives. A keyword scan of script members identified the later `grids_compared.py`, which describes the matched-stretch experiment; it did not identify the producer of these three historical diagnostics. Available velocity files were the separate 48-grid matched-stretch snapshots. No matching original pressure snapshots or complete ratio/integral producer was identified. This is a bounded search of the available deliveries, not proof that the records do not exist elsewhere.

The [archive search inventory](diagnostic-source-search.json) records the inspected filenames, archive hashes, matches and array member names. Historical files remain unchanged. The [earlier intake pages](archive/README.md) are retained for provenance.

## Recheck without rerunning a simulation

From the repository root:

```sh
python math/imports/hug-ns/diagnostic_audit.py --out scratch/historical-diagnostic-audit.json
```

Create `scratch/` first. The script checks saved table arithmetic and versioned source hashes. It performs no time evolution and does not upgrade an unreproduced trajectory to verified evidence.
