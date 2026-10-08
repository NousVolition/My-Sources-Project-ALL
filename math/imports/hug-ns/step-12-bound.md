# The proposed short-time estimate

[Full argument](BOUND.md) · [Review of this assumption](RESULTS-REVIEW.md#review-of-the-proposed-starting-time-bound)

**The reciprocal calculation is correct under the stated assumption. That assumption has not been established for this flow.**

Take a starting maximum vorticity `B(0) = 241.7`. If a nonnegative quantity `B` satisfies

```text
B′ ≤ B²
```

then scalar comparison gives

```text
B(t) ≤ B(0) / (1 − B(0)t), for 0 ≤ t < 1/B(0)
1/B(0) ≈ 0.004137
```

This is a conditional calculation. Navier–Stokes vorticity stretching depends on the velocity gradient. A derivation is needed before replacing that dependence with the unit-coefficient inequality `B′ ≤ B²`.

The saved run reaches time 0.4 with finite reported values. That observation does not prove the assumed inequality. The number 0.004137 is therefore not an established lifespan bound for this flow, and it is not a predicted failure time.
