# Intermediate obstacle-only contact stage — self-crossing found

[Final repaired model and results](../README.md). [Exact stage summary](data/stage-summary.json).

This stage replaced the former penetration spring with object and square-wall barriers, but did not include material self-contact. It completed **34 configurations** (28 from its main protocol and 6 from its sensitivity protocol). A broader audit found proper crossings of nonadjacent material edges in **11** of those configurations, despite positive object clearance, positive cell areas and passing energy budgets. Object exclusion alone was insufficient.

For example, in `amplitude4-one.npz`, crossing is already present at model time 4.6. The final model uses this recording as a rejection regression check. The final model also reruns the case with self-contact and compares it against independent DOP853 integration.

The main stage was stopped when these crossings were discovered. Its remaining planned cases have no completed result; they must not be counted as successes or numerical failures. The frozen plans remain in `data/protocol.json` and `data/sensitivity-protocol.json`. Individual `.json` files record actually completed cases; matching `.npz` files contain full states. `data/audits` gives recomputed errors, dense edge probes and the proper crossing count summed across saved frames. That count is not a count of distinct physical collision events.

To reproduce the failed audits from this directory:

```text
python analyze_results.py --audit-ready
```

It is expected to return a failing status and list 11 configurations. The old report builder requires passing geometry, so no successful report was generated for this stage. `MODEL-NOTES.md` retains working notes and the originally planned workflow; its uncompleted placeholders are not completed results.

`data/initial-attempts` retains two recording-time roundoff stops and two completed early trials. `data/low-barrier-attempts` retains eight completed weaker-barrier trials, including an energy-budget miss. The independent DOP853 recordings here were completed, but the queued matching fine-step tests had not finished when the stage was stopped, so no completed independent-comparison result is claimed for them.
