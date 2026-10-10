# Current code review status

[Current overview](README.md) · [Full earlier code review](archive/earlier-import/REVIEW.md) · [Historical archive](archive/README.md)

The earlier detailed review is archived. Its findings remain evidence about the preserved imported implementation.

| Finding | Current disposition |
| --- | --- |
| Reversed pressure sign and inconsistent carrying term | Corrected diagnostic is available; the earlier pressure JSON still needs its original velocity fields for recalculation |
| Corner in the initial gate | A smooth product gate is implemented in `corrected_v1`; this explicitly changes the starting field |
| Aliased retained cutoff | `corrected_v1` uses the strict integer cutoff and has convolution checks |
| Rounded or stale time-step estimate | `corrected_v1` monitors the estimate at each step and logs accepted/rejected trials |

Read the [separate repair and its saved checks](corrected_v1/README.md). Its 28 regression tests passed. The three expected failures in the earlier implementation remain deliberate records of those earlier defects. The repair's bounded checks do not establish long-time or spatial convergence.

[Current status of the pressure, integral and ratio tables](DIAGNOSTIC-STATUS.md).

<details>
<summary>Links to sections in the archived page</summary>

<a id="review-of-the-supplied-hug-ns-files"></a>

[Review of the supplied hug-ns files](archive/earlier-import/REVIEW.md#review-of-the-supplied-hug-ns-files)

<a id="what-was-corrected-in-the-same-point-calculation"></a>

[What was corrected in the same-point calculation](archive/earlier-import/REVIEW.md#what-was-corrected-in-the-same-point-calculation)

<a id="results-supplied-during-the-review"></a>

[Results supplied during the review](archive/earlier-import/REVIEW.md#results-supplied-during-the-review)

<a id="the-invisible-boundary"></a>

[The invisible boundary](archive/earlier-import/REVIEW.md#the-invisible-boundary)

<a id="how-this-code-differs-from-the-existing-project"></a>

[How this code differs from the existing project](archive/earlier-import/REVIEW.md#how-this-code-differs-from-the-existing-project)

<a id="findings-that-need-correction-before-relying-on-larger-runs"></a>

[Findings that need correction before relying on larger runs](archive/earlier-import/REVIEW.md#findings-that-need-correction-before-relying-on-larger-runs)

<a id="1-the-underlying-gate-has-corners-in-its-transition-region"></a>

[1. The underlying gate has corners in its transition region](archive/earlier-import/REVIEW.md#1-the-underlying-gate-has-corners-in-its-transition-region)

<a id="2-the-fourier-cutoff-retains-an-alias-at-n-divisible-by-three"></a>

[2. The Fourier cutoff retains an alias at N divisible by three](archive/earlier-import/REVIEW.md#2-the-fourier-cutoff-retains-an-alias-at-n-divisible-by-three)

<a id="3-the-runner-can-exceed-its-own-step-estimate"></a>

[3. The runner can exceed its own step estimate](archive/earlier-import/REVIEW.md#3-the-runner-can-exceed-its-own-step-estimate)

<a id="additional-limits-in-the-measurements"></a>

[Additional limits in the measurements](archive/earlier-import/REVIEW.md#additional-limits-in-the-measurements)

<a id="review-of-boundmd"></a>

[Review of BOUND.md](archive/earlier-import/REVIEW.md#review-of-boundmd)

<a id="checks-performed-on-the-reviewed-version"></a>

[Checks performed on the reviewed version](archive/earlier-import/REVIEW.md#checks-performed-on-the-reviewed-version)

<a id="included-files"></a>

[Included files](archive/earlier-import/REVIEW.md#included-files)

</details>
