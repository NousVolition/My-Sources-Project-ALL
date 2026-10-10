# Hug flow: current reading guide

[Math index](../../README.md) · [Reports](../../../reports/README.md) · [Historical archive](archive/README.md)

## Start here

- **Current code repair:** [the separately versioned corrected solver](corrected_v1/README.md), its recorded checks and its limitations.
- **The three older diagnostic questions:** [pressure, vorticity integral and ratio-rise status](DIAGNOSTIC-STATUS.md).
- **Completed mirror work:** [four supplied fluid tests](mirror-fluid/README.md) and [20-run signed calibration](mirror-calibration/README.md), each with its own starting field and evidence.
- **Larger studies:** [matched-stretch results](../../../reports/matched-stretch/README.md) and [completed vortex matrix](../../../reports/study/final-matrix/README.md).

## Current code status

The pressure diagnostic uses the corrected sign and the solver's filtered carrying term. The separately versioned `corrected_v1` also repairs the smooth-gate construction, retained cutoff alias and step scheduler. Its 28 repair-specific regression tests passed; the saved verification states the short-run scope. The replacement gate changes the initial field, so results from the earlier and repaired versions must retain their version labels.

## Earlier files are preserved

The older overview and two detailed reviews now live in the [historical archive](archive/README.md). Their previous page addresses continue to work. Raw [result tables](results/), imported source, [provenance](provenance.json), verification records and original PDF remain at their recorded paths. The root `solver.py` and `run.py` belong to the earlier imported version; use the versioned repair documentation when reviewing the corrected implementation.

The old pressure table cannot be regenerated from its summary alone. The reported vorticity integral lacks its internal quadrature record, and the ratio-rise table lacks its producing configuration. These gaps are listed precisely in the [diagnostic status](DIAGNOSTIC-STATUS.md).

The conditional mathematical estimate remains in [BOUND.md](BOUND.md). A finite sampled-grid integral does not establish a continuum regularity bound.

<details>
<summary>Links to sections in the archived page</summary>

<a id="hug-flow-code-and-results"></a>

[Hug flow: code and results](archive/earlier-import/README.md#hug-flow-code-and-results)

<a id="completed-mirror-tests"></a>

[Completed mirror tests](archive/earlier-import/README.md#completed-mirror-tests)

<a id="earlier-controls"></a>

[Earlier controls](archive/earlier-import/README.md#earlier-controls)

<a id="earlier-results"></a>

[Earlier results](archive/earlier-import/README.md#earlier-results)

<a id="what-the-boundary-does"></a>

[What the boundary does](archive/earlier-import/README.md#what-the-boundary-does)

<a id="equation-and-measurements"></a>

[Equation and measurements](archive/earlier-import/README.md#equation-and-measurements)

<a id="run-the-code"></a>

[Run the code](archive/earlier-import/README.md#run-the-code)

<a id="what-is-established"></a>

[What is established](archive/earlier-import/README.md#what-is-established)

</details>
