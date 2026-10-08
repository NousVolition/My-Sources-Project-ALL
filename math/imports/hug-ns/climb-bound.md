# Climbing flow: observed rates and a conditional bound

[Saved observations](results/matched-nostop.json) · [Vorticity argument](BOUND.md) · [Result review](RESULTS-REVIEW.md#climbing-start-bound)

**The saved growth-rate calculation is reproducible. The proposed differential bound still needs a derivation.**

## What the saved samples show

Maximum vorticity starts at `B(0) = 60.198351`. The no-stop table reaches time 0.35, with a largest saved value of `611.139126` and a final value of `513.808318`.

For each adjacent pair of saved times, define

```text
q_i = (B_(i+1) − B_i) / ((t_(i+1) − t_i) B_i²)
```

This compares the average increase over that interval with the square of the value at its start. The largest recorded `q_i` is **0.16394196**, or **16.39%**, on the interval **0.0747–0.0997**.

Later, the interval 0.1993–0.2242 gives **1.12%**. Three saved intervals have decreases. The final short interval gives **2.13%**, so the later rates do not all stay below 2%.

These are finite differences between observations. They do not measure the largest instantaneous growth rate inside each interval.

## What the proposed bound assumes

If the actual maximum satisfied `B′ ≤ B²`, scalar comparison would give

```text
B(t) ≤ B(0) / (1 − B(0)t), for t < 1/B(0)
1/B(0) ≈ 0.01661175
```

The reciprocal uses only the starting value. Its application to this flow requires proving the assumed differential inequality. Finite observations after 0.0166 do not prove that assumption. This calculation does not predict a failure at 0.0166 or establish smoothness for all time.
