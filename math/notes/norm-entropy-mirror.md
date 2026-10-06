# Norm and Entropy: two mirrored halves

[Math guide](../README.md) · [Experiment](../mirror_family.py) · [Saved results](../results/mirror-family.json)

**Contributor's idea:** give the test two named personalities, Norm and Entropy, and look for two halves of one data family.

Here those names identify the original curve and its flipped neighbor. We test a correspondence between numerical samples; we have not measured personality or mathematical entropy.

## What we built

We reuse the maintained smooth field. At `y = z = 0`, take 65 equally spaced samples from `x = -3` to `x = 3`. Each point has two features: position `x` and velocity `u_x`, both nondimensional.

- **Norm:** the original samples.
- **Entropy:** the same samples shifted right by one box length and reflected vertically using the existing endpoint-matched construction.
- **Alignment:** undo that known shift and reflection before comparing each pair.

The mirror relationship is built into the example. This is a controlled check, not a discovery of an unknown relationship in independent observations.

## What happened

| Check | Result |
| --- | --- |
| Largest paired difference after alignment | 0 in the saved calculation |
| Largest change in distances within either half | 0 in this two-feature slice |
| Corresponding pairs in the same K-means group, before alignment | 26 of 65 |
| Corresponding pairs in the same K-means group, after alignment | 65 of 65 |

For the illustrative K-means comparison, we explicitly request **two groups**. We standardize features once using the placed data, retain those scales after alignment, and run seeds 0–4 with 20 initializations each. The half names and pair IDs are not clustering features.

Before alignment, one group contains 65 Norm points and 26 Entropy points; the other contains 39 Entropy points. After alignment, both groups contain equal numbers from each named half: 33 + 33 and 32 + 32. All five seeds give those compositions; numeric group IDs can swap. Thus the grouping does not simply reproduce the supplied names. The requested two-group split is not evidence for two natural personalities.

**What this lets us study:** whether a grouping changes because of how related data are positioned or represented. The next useful extension would be to add small, recorded differences to one half and measure how well alignment still identifies its corresponding pairs.

## What still happens at the physical join

The [reflection probe](../reflection_probe.py) checks the earlier boundary idea separately. Over a 33 × 33 sample of the joining face:

| Copy operation | Largest velocity mismatch | One-sided center-line slopes | Interior divergence at (4.2, 0, 0) |
| --- | ---: | --- | ---: |
| Repeat | 0.0537094 | +0.164579 / +0.164579 | about 0 |
| Reverse signs | 0.135246 | +0.164579 / −0.164579 | about 0 |
| Match endpoints, then flip | 2.17 × 10⁻¹⁹ | +0.164579 / −0.164579 | −0.329145 |

These are different quantities: velocity mismatch, velocity derivative, and divergence. The saved slopes and divergence use finite differences with step `1e-5`, with refinement and an analytic checkpoint recorded in [the results](../results/reflection-probe.json).

Matching the endpoint values closes the gap, but leaves a corner and introduces nonzero divergence in the three-dimensional copy. Alignment for comparing data does not fix this physical boundary issue. Away from the chosen center slice, the reflection offset depends on position; global rigid-shape equivalence is not established. Neither experiment steps Navier–Stokes or establishes a singularity result.

## Reproduce

From the repository root, after installing `requirements.txt`:

```sh
python math/reflection_probe.py --json scratch/reflection-probe.json
python math/mirror_family.py --json scratch/mirror-family.json
python -m pytest math/tests/test_reflection_probe.py math/tests/test_mirror_family.py -q
```

The nine checks cover face matching, slope reversal, divergence, exact paired correspondence, and a negative control that changes one copied sample and breaks the match.

Method references: [K-means](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html) and [StandardScaler](https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html). These synthetic field samples are separate from the project's [conversation clustering study](../../language/clustering/README.md).
