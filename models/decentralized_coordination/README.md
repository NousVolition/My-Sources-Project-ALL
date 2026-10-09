# Coordination with the freedom to refuse

**Simulation study and optional volunteer protocol. No human data were collected.**

Can a group coordinate through local interactions while its members retain the ability to refuse? Does apparent leadership follow a person or their position? This study makes those questions testable without identifying a behavioral rule with actual water or claiming to prove or disprove Freud.

[Open the saved report](results/report.html) · [Methods](METHODS.md) · [Volunteer protocol](PROTOCOL.md) · [Sources](SOURCES.md) · [Machine-readable results](results/summary.json)

Download or clone the repository, then open `results/report.html` in a browser. GitHub's HTML file view displays source. The report, plots, tables, and code work locally; no service or account is needed.

## Conditions

| Condition | Nine-agent implementation |
| --- | --- |
| A: fixed hierarchy | The leader commands their initial color. Followers act every 10 seconds; primary compliance is 100%. Leadership identity is fixed within a round and randomized across blocks. |
| B: deterministic local organization | Nine-seat ring; every 10 seconds, simultaneously switch only if both immediate neighbors have the opposite color. Otherwise hold. No leader or communication. |
| C: voluntary local coordination model | Same ring. Agents can imitate a neighbor, decline a differing proposed color, or initiate a switch. Timing and probabilities are specified assumptions, varied in sensitivity analyses. |

All primary rounds use independent fair red/blue starting colors and stop at first unanimity or 120 seconds. Starting states, seating, and potential designated leader are matched within each three-condition block. Already-unanimous starts are retained at time zero. There are four separate pilot blocks (four rounds per condition), followed by 4,000 main blocks.

**B is a water-inspired analogy, not molecular dynamics or a claim about actual water. C is a stochastic toy model of choices, not a model or measurement of metaphysical free will.**

## Reproduce

Use Python 3.12. From this directory:

```console
python -m venv .venv
# Activate with .venv\Scripts\activate on Windows, or source .venv/bin/activate on macOS/Linux.
python -m pip install -r requirements.txt
python run_study.py --output reproduced
python verify_results.py --actual reproduced --expected results
```

The runner refuses to overwrite a nonempty output directory. Its defaults reproduce seed `20261008`, 12 pilot rounds, 12,000 main rounds, 25,200 sensitivity rounds (21 scenarios × 1,200 starts), and 8,640 crossed influence rounds, plus exhaustive enumeration of 512 B starts. A smaller exploratory run:

```console
python run_study.py --output quick-check --blocks 60 --sensitivity-blocks 30 --influence-replicates 10
```

`--no-plots` avoids importing Matplotlib and omits the HTML/figures; NumPy is still required for analysis. The dynamics alone in `simulation.py` use only the standard library. Seeds are derived through SHA-256 and never depend on Python's randomized object hashes. Exact data reproduction targets the saved Python/minor version and direct dependency versions; plot bytes can differ across rendering environments. `manifest.json` records the actual versions and hashes of sources and data.

From the repository root, run the study tests using only NumPy and the standard library:

```console
python -m unittest discover -s models/tests -p test_coordination.py -v
```

Set `PYTHONPATH=models` for that unittest command (PowerShell: `$env:PYTHONPATH='models'`; POSIX: `PYTHONPATH=models python ...`). Alternatively, install the root test dependencies and run `python -m pytest -q`; the existing root configuration includes these tests. The tests cover an independent majority-rule oracle for all 512 states, symmetries, pinned-pair behavior, censoring, tied initiators, credit conservation, crossed seating, recovery of known effects, and a small end-to-end run.

## Read the outputs

| File | Purpose |
| --- | --- |
| `results/report.html` | Interpreted report with four figures and comparison tables. |
| `results/rounds.csv` | Pilot and main records; filter `phase` before summarizing. Arrays use JSON inside CSV fields. |
| `results/pilot-traces.json` | Every simulated pilot decision, proposal, refusal and observed state. |
| `results/*-schedule.csv` | Starting states, identity-to-position mappings, leader assignment and actual randomized condition order. |
| `results/exact-B.csv` | Attractor and period for every one of the 512 starting patterns. |
| `results/sensitivity-*.csv` | Round data and summaries for changed assumptions, including lower A compliance. |
| `results/influence-observations.csv` | One row per agent per crossed round; nine agent rows are not independent replications. |
| `results/influence-summary.json` | Identity/position effects, held-out prediction and top-copy-source persistence. |
| `results/paired-contrasts.csv` | Paired differences and bootstrap intervals resampling matched blocks. |
| `results/manifest.json` | Sizes, seed, parameters, versions, repository base, source and data hashes. |
| `forms/` | Blank volunteer event, round and private self-report templates, without participant data. |

Identity numbers and seat numbers are zero-based (0–8); red is 0 and blue is 1. Any screenshots/plots marked “pilot” also contain simulated data. No real identities, choices or participation records belong in the public repository.

## Scope of the conclusion

The exact B result is 152/512 consensus starts (29.6875%), with 360 frozen disagreements and no persistent cycles on this odd ring. The primary C result and its sensitivity range are reported from the saved run. The ring itself has no structurally special seat. A separate **C-only star-network extension** tests an unequal topology; it does not silently redefine condition B. Heterogeneous-identity runs insert a particular identity advantage deliberately.

Thus position can explain apparent leadership in one constructed mechanism, and identity can matter in another. No estimate here decides between those explanations in people. First movement, being copied, authority, and intentional leadership are distinct outcomes. A 120-second consensus task also does not test survival, welfare, coercion, psychoanalytic theory, water physics, or long-term stability of agreement.
