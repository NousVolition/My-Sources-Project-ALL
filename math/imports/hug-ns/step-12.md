# Step 12: a proposed bound on maximum vorticity

[Vorticity estimate](BOUND.md) · [Data and fitted constants](RESULTS-REVIEW.md#review-of-step-12md)

**The proposed formulas fit the saved samples with the constants below. Neither constant has been derived from the equation.**

Let `B` be maximum vorticity and `A` its RMS over the fixed periodic volume. The normalization of `average` must be confirmed in the producing script before applying this definition to every table.

## Candidate estimates

The first proposal is

```text
B(t) ≤ C A(t)
```

Here `C` may depend on the full starting field but must remain fixed in time. In [ratio-rose.json](results/ratio-rose.json), `B/A` starts at 30.62182 and reaches **32.2009963**. The starting ratio is too small to bound the later samples. A fitted constant at least as large as the saved maximum covers those samples through time 0.4.

The second proposal is

```text
B(t) ≤ C A(t) (1 + log(1 + A(t)))
```

With the file's numerical convention, the smallest fitted constant is **10.2554604**; 10.3 covers the saved rows. A dimensionless convention or reference scale for `A` must be specified inside the logarithm. At the same `C`, this formula allows larger `B` than the linear formula.

## Why this does not finish the argument

Dropping the nonnegative viscous term from the vorticity energy estimate gives

```text
d/dt A² ≤ 2 B A²
```

Even if `B ≤ C A` were proved, it would give only `A′ ≤ C A²` when `A > 0`. The scalar equation `a′ = C a²` can become unbounded in finite time. This estimate therefore does not rule out such growth; it also does not show that the actual flow becomes unbounded.

The next mathematical requirement is control of the time integral of `B`, as explained in [BOUND.md](BOUND.md). Fitting constants to saved samples does not supply that control between samples or after the recorded interval.
