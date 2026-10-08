# What the vorticity bound requires

[Experiment overview](README.md) · [Numerical observations](RESULTS-REVIEW.md#the-vorticity-integral-needs-its-calculation-record) · [Mathematical review](REVIEW.md#review-of-boundmd)

**The estimate below is valid for a smooth, periodic, divergence-free solution. The missing step is a bound on the time integral of its maximum vorticity.**

## Starting from the equation

The gate prepares the initial data. The subsequent equation is

```text
∂t u + (u · ∇)u = −∇p + ν Δu
∇ · u = 0
```

Let `ω = curl u` be vorticity. Taking the curl gives

```text
∂t ω + (u · ∇)ω = (ω · ∇)u + ν Δω
```

The pressure term disappears because the curl of a gradient is zero. The two terms on the right describe stretching and viscous diffusion.

For the periodic domain, the vorticity energy estimate is

```text
d/dt ||ω||₂² + 2ν ||∇ω||₂² ≤ 2 ||ω||∞ ||ω||₂²
```

`||ω||∞` is the maximum vorticity magnitude. `||ω||₂` is its volume-integrated root-square norm. The estimate follows from the vorticity energy identity and `||∇u||₂ = ||curl u||₂` for periodic, divergence-free fields.

## The step still required

For a smooth solution approaching a finite time `T`, the continuation criterion uses

```text
∫₀ᵀ ||ω(t)||∞ dt < ∞
```

A bound on this integral allows smooth continuation past `T`. Applying the argument for all time requires control for every finite `T`. A result covering every admissible starting field must establish that control for each such field. See the [Navier–Stokes continuation criterion discussed by Chemin and Zhang, p. 133](https://smf.emath.fr/sites/default/files/2024-05/ens_ann-sc_49_131-167__sample.pdf).

## What the saved table establishes

For the reported 48³ pulled-tube run, maximum vorticity starts at about 71.7 and ends at about 29.3 at time 2. The file records an integral of **86.12856**. Integrating only its saved rows by the trapezoidal rule gives **88.67728**.

The file does not include the internal integration samples or quadrature rule. That calculation record is needed to reconcile the two values. Even a fully reproduced finite grid integral is a numerical observation, not the required bound on the continuous solution.

The estimate is established. The bound needed to complete this argument remains open.
