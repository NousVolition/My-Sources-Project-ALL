# Review of the supplied result files

[Code review](REVIEW.md) · [Results](results/) · [Verification record](received-results-check.json)

The latest `hug-ns (5).zip` contains the same files, byte for byte, as `hug-ns (4).zip`; only the ZIP container differs. This update adds one copy of eight supplied JSON files and the README's “Later checks” section. The received solver matches the earlier upload. The previously published pressure-diagnostic correction is retained.

Later separate attachments are also included: [ratio-rose.json](results/ratio-rose.json), [stretch-at-max.json](results/stretch-at-max.json), and the revised [step-12.md](step-12.md). The later [stretch-rate.json](results/stretch-rate.json) and [high-fraction.json](results/high-fraction.json) plus [matched-strain.json](results/matched-strain.json) are followed by [matched-continue.json](results/matched-continue.json), then [matched-nostop.json](results/matched-nostop.json) and [matched-finer.json](results/matched-finer.json), bringing the total to sixteen result tables. Their source hashes are recorded separately from the ZIP entries; the superseded step-12 source hash is retained.

## Results at a glance

- All sixteen files pass checks for finite values, increasing saved times and the arithmetic present in each file.
- The 48³ and 64³ starting fields were rebuilt. The later trajectories were not replayed.
- The no-stop table reaches time 0.35 with a largest saved maximum vorticity of 611.139126. The newer N=64 table reaches 0.12 and ends at 200.615406.
- Two calculations need their source records: the historical pressure split and the vorticity time integral. Their specific issues are described below.

## What was checked

- All numeric values are finite and all saved time sequences increase.
- The reported ratios equal `biggest / average`; pressure-table sums equal their supplied three terms; same-spot differences and counts agree.
- Reported start, end, and high-peak summary fields agree with their saved samples where comparable.
- Rebuilding the **starting fields only** at 48³ and 64³ reproduces initial speed, vorticity maximum, and vorticity RMS to numerical roundoff. The 64³ initial energy also matches. No trajectory was rerun.
- Every result file is preserved byte for byte and recorded in [provenance.json](provenance.json).

These checks validate file integrity, arithmetic, and starting values. Full velocity checkpoints, producing scripts for the added diagnostics, exact step schedules, and quadrature details are absent. The later states have not been independently reproduced.

## Recorded results

| Supplied file | Recorded interval | What its saved rows show |
| --- | --- | --- |
| [tube-finer.json](results/tube-finer.json) | 64³, 0–1.6 | Speed 3.59974 → sampled high 6.15678 → 3.34135; saved energy decreases |
| [two-tubes.json](results/two-tubes.json) | 48³, 0–2 | Speed 2.97769 → 2.18623; no saved peak exceeds the start; saved energy decreases |
| [same-spot.json](results/same-spot.json) | 48³, 0–2 | All 17 saved smoothing magnitudes are smaller than their same-spot carrying magnitudes; these are unsigned magnitudes |
| [pressure-at-peak.json](results/pressure-at-peak.json) | 48³ to 2; 64³ to 1.6 | Supplied signed pressure values include both signs; all saved smoothing terms are negative; see the diagnostic warning below |
| [spin-ratio.json](results/spin-ratio.json) | 48³, 0–2 | Maximum/RMS vorticity ratio 21.12167 → 11.18506 |
| [spin-ratio-64.json](results/spin-ratio-64.json) | 64³, 0–1.6 | Ratio 20.02316 → 16.55321 |
| [spin-ratio-low-nu.json](results/spin-ratio-low-nu.json) | Viscosity 0.0005, 0–1 | Ratio 21.12167 → 13.58004; final speed 4.81185 exceeds initial 3.74068 |
| [bound-check.json](results/bound-check.json) | 48³, 0–2 | Reported accumulated maximum-vorticity integral 86.12856; final maximum 29.30359 |
| [stretch-at-max.json](results/stretch-at-max.json) | 0–0.24 | Seven positive and three negative stretching entries; all ten smoothing entries are negative |
| [ratio-rose.json](results/ratio-rose.json) | Tight tube in a stretch, 0–0.4 | Ratio 30.62182 → sampled high 32.20100 → 21.32382; grid, viscosity and exact starting parameters are not recorded |

“Average” matches the root mean square of vorticity magnitude at the verified starting states. It is not the arithmetic mean magnitude or the unnormalized L2 norm. The ratio decreases from first to last in these tables, but each series has intervening increases. The README's phrase “the hole does not open” has no separately defined or recorded gap measurement here.

## Pressure values need recalculation

The new ZIP still supplies the earlier `pressure_at_peak` function with its reversed Fourier pressure sign. The received pressure JSON is kept as historical data alongside that source provenance. It is **not** output from the corrected function on GitHub. Its listed `sum` is arithmetically consistent, but that does not make it the implemented equation's local energy rate.

Current code corrects that sign and uses the same filtered carrying term as the solver. Correcting the saved pressure signs alone is therefore insufficient to reconstruct a fully consistent total from these summary tables. Original velocity snapshots or a reproducible replay are needed to regenerate the corrected split. Neither was supplied. The existing grid-filter and step-scheduling issues also remain documented in [the code review](REVIEW.md).

The saved pressure values change sign, and reversing every sign would still leave a mixture. This observation alone does not validate the magnitude or total rate of the reported split.

## The vorticity integral needs its calculation record

Using only the saved vorticity maxima and rounded times in `bound-check.json`, trapezoidal integration gives **88.67728**, while the file records **86.12856**. This difference could reflect integration at more frequent internal steps or another quadrature rule. The archive does not specify which, so the recorded integral is not independently reproduced.

A finite discrete integral does not prove regularity of the continuous equation. The clarification in [BOUND.md](BOUND.md) remains applicable now that the sample tables are available.

The different files also use different observation times. For example, the saved 64³ peak is 6.15678 in `tube-finer.json` and 6.06010 in the pressure table; their sample times differ. These are separate reported observations, not an established matched-time grid-convergence comparison.

## The ratio-rise attachment

The arithmetic of `ratio-rose.json` is consistent. At its maximum ratio (`32.200996` at time `0.1977`), maximum vorticity is `241.458690`, slightly below its starting `241.712268`; the recorded average fell from `7.893465` to `7.498485`. The ratio rise at that point is driven by the denominator falling. Maximum vorticity later reaches `250.811878` at time `0.3163`. These are distinct observations.

This file supplies no grid size, viscosity, complete initial-field parameters, step size, or average-norm definition. Its numbers can be checked arithmetically, but the earlier 48³/64³ initial-field reproductions do not validate this different start.

## Review of step-12.md

The revised `step-12 (1).md` replaces the earlier note. It now permits a constant depending on the entire starting field. The earlier objection to a universal finite bound using only the current L2 norm does not, by itself, rule out this different claim.

The constants can be checked against the supplied samples:

| Candidate | Smallest fitted constant covering the saved rows |
| --- | --- |
| `B <= C A` | `32.200996292744684` |
| `B <= C A (1 + log(1+A))` | `10.255460414336675` |

The stated `32.2` is a rounded display of the observed maximum; a strict numeric ceiling of exactly 32.2 misses it slightly. The stated `10.3` covers the recorded samples for the logarithmic candidate. Neither is an estimate derived from the equation or a bound between samples or beyond time 0.4. The logarithmic candidate permits larger B as A increases, so it is a looser bound than the linear candidate at the same C.

The revised note correctly observes that `dA/dt <= C A^2` does not rule out finite-time growth without a bound. This is a limitation of that estimate; it is not a demonstrated singularity in the actual flow. A proved finite-time integral bound for the actual solution would be needed for the continuation argument described in BOUND.md.

The revised note defines A as RMS; the producing script must confirm that the ratio-rise table uses that convention. Earlier tables' initial “average” values match RMS; the ratio-rise file does not specify its normalization. On a volume-216 cube, RMS equals the unnormalized L2 norm divided by `sqrt(216)`. A logarithmic expression using `1+A` also needs a dimensionless convention or a reference scale.

## Stretch-at-maximum attachment

[stretch-at-max.json](results/stretch-at-max.json) contains ten finite samples through time 0.24. Stretch is positive in seven and negative in three; smoothing is negative in all ten. The reported maximum goes from `241.712268` to `235.545895`, reaching a saved high of `243.038072` in between.

The file does not include the producing formula, units, grid, viscosity, time-step schedule, or selected point coordinates. These values cannot yet be interpreted as a verified complete balance for the derivative of maximum vorticity. In particular, sampling a grid maximum is different from following one point or locating the exact continuous maximum. No equation or pressure-sign change was made based on this table.

## Signed stretching-rate attachment

[stretch-rate.json](results/stretch-rate.json) supplies ten samples through time 0.4 and a final signed integral of `-0.1882288755`. At time zero, its rate matches the earlier table's `stretch / biggest`. That single match supports the label but does not supply the missing producing formula for every sample.

This is a signed rate integral. It is different from the nonnegative integral of maximum vorticity used in the continuation condition in BOUND.md. A negative value here does not establish that condition. The file omits the quadrature rule and internal samples; a check of its saved rows is retained in the verification record.

## Review of the proposed starting-time bound

[step-12-bound.md](step-12-bound.md) uses `B(0) = 241.7` and the hypothetical inequality `B′ <= B^2`. Under that assumption, scalar comparison gives `B(t) <= B(0)/(1-B(0)t)` for `t < 1/B(0)`, approximately `0.004137` in these units.

The assumption is not established by the supplied files. In Navier–Stokes, vorticity stretching depends on the velocity gradient; replacing that dependence with the stated unit-coefficient bound in B alone needs a derivation. The reciprocal computation is correct as a conditional ODE estimate. It has not been proved to be a lifespan bound for this flow. A finite grid run to time 0.4 does not verify the missing inequality or certify the continuum lifespan.

The original note treated 0.0041 as a proved lower bound. The revised note states the assumption explicitly: the inequality still needs a derivation. This time is not a predicted failure of the flow.

## High-fraction attachment

[high-fraction.json](results/high-fraction.json) contains ten finite samples through time 0.15. Every `fraction` equals `rate / biggest` to the recorded precision. There is one negative entry (`-0.052377` at time zero); the largest positive fraction is `0.027169`. This signed ratio is not a fraction of spatial volume or a probability.

Its maximum-vorticity column starts at `128.907398`, reaches a saved high of `133.675275`, and ends at `106.860564`. This is a different starting value from the preceding tight-tube table. The file does not define the rate formula, units or full run settings, so no common trajectory or universal growth bound is inferred from the ratio.

## Matched-strain attachment

[matched-strain.json](results/matched-strain.json) contains eight finite samples through time 0.08. Its reported maximum increases from `60.198351` to `102.461848`. All fractions match `rate / biggest`; the first is negative and the remaining seven are positive. The largest positive fraction is `0.279297`.

As with the preceding fraction table, the producing rate formula and full run settings are missing. The different initial maximum identifies a distinct numerical record; the file alone does not establish a matched control or a complete maximum-growth balance.

## Matched continuation attachment

[matched-continue.json](results/matched-continue.json) supplies nine finite rows from time 0 to 0.1999. Its initial row exactly matches `matched-strain.json`. The remaining saved times differ, so this is not a row-for-row extension of that earlier table.

The reported maximum increases at every saved row, from `60.198351` to `462.301472`. The file also reports a high of `475.544253` and the stop reason `biggest passed 400`. That high exceeds every saved row; its time and state are not included. All saved fractions agree with `rate / biggest`; two are negative and seven positive.

The table documents numerical growth and a reported stopping threshold. It does not establish breakdown of the continuous equation. The producing script, complete settings, velocity checkpoints and resolution checks are absent, so the cause of the growth and continuity of the underlying run are not independently verified.

## Matched no-stop attachment

[matched-nostop.json](results/matched-nostop.json) contains sixteen finite samples through time 0.35. The initial reported maximum is `60.198351`, exactly matching the preceding two matched tables. Subsequent sample times differ. The file does not contain a rate or fraction column.

The largest saved maximum is `611.139126` at time `0.3239`; the last is `513.808318`. The series includes rises and falls. It supplies observations beyond the preceding file's reported 400 stopping threshold and ends with the reason `time 0.35`.

Crossing that chosen threshold was therefore not, by itself, evidence of a numerical failure. This table still omits the grid, viscosity, step schedule, source revision and full velocity fields. It does not verify that the same trajectory continued, distinguish resolved growth from numerical error, or establish a singularity.

## Matched finer-grid attachment

[matched-finer.json](results/matched-finer.json) records `N=64` and eight finite samples through time 0.12. Maximum vorticity increases at every saved time, from `59.279954` to `200.615406`.

The preceding matched tables start at `60.198351` and use different sample times. The new file provides a grid size, but omits viscosity, initial-field parameters, time-step settings and full velocity fields. The saved tables therefore record growth in both sets of observations; they do not yet establish a controlled resolution comparison or a convergence rate.

## Climbing-start bound

[climb-bound.md](climb-bound.md) compares the no-stop table with a proposed square-growth estimate. The largest forward difference, divided by the square of the value at the interval start, is **0.16394196** on **0.0747–0.0997**. This reproduces the note's rounded 16% observation. It is an interval average, not a measured maximum instantaneous derivative. The final short interval gives 2.13%, so the later values do not all remain below 2%.

The reciprocal of the starting maximum is **0.01661175**. Its use as a lifespan estimate depends on deriving `B′ ≤ B²` for the actual solution. That derivation is absent. The revised note states the assumption and calculation separately; the original source hash and exact wording changes are preserved.

## PDF report

[hug-runs.pdf](hug-runs.pdf) is preserved byte for byte. Both pages were rendered and inspected; the table and plotted curve are readable. The figure has no numeric vertical tick labels, so the JSON tables remain the source for exact values.

The report's rounded no-stop values (60, 611 and 514 at time 0.35) and N=64 final value (201 at time 0.12) agree with the corresponding saved tables. Its rounded 16% growth figure is reproduced by the interval calculation in [climb-bound.md](climb-bound.md).

The following statements need qualification:

- **“Not a coarse-grid effect”:** growth is recorded at N=64, but matched times, full run settings and resolution-error estimates are missing. This does not establish convergence or rule out numerical error.
- **“Also falls with no smoothing”:** no zero-viscosity result table was supplied. `spin-ratio-low-nu.json` records viscosity 0.0005, which is positive.
- **“Earliest blowup” at 0.0166:** the reciprocal is conditional on `B′ ≤ B²`; that inequality has not been derived for this flow.
- **“The hug is the smooth start”:** the finite Fourier field is smooth, but the written gate formula has a corner in its transition region. The [code review](REVIEW.md#1-the-underlying-gate-has-corners-in-its-transition-region) distinguishes the two.
- **Pressure at the fastest point:** the saved pressure table uses the older diagnostic. Its signs and magnitudes need a consistent recalculation before use in the local energy balance.

The ring-to-time-40 and other longer-run summaries repeat claims from the uploaded README. Their full trajectories were not supplied or rerun in this review. The PDF correctly leaves the all-time bound open.

## Reproduce these checks

From the repository root, after installing its pinned requirements:

```sh
python math/imports/hug-ns/check_received_results.py --check-starts --out scratch/hug-ns-received-check.json
python -m pytest -q math/tests/test_hug_ns_received_results.py
```

Create `scratch/` first if needed. `--check-starts` constructs initial fields and advances zero time steps. Omit it to check only the supplied tables and their arithmetic.
