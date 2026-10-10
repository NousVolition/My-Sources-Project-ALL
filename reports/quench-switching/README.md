# Quench and switching-identity diagnostics

Real-data analysis of peak-vorticity timing and material-marker identity changes, using published fields from the central-response study.

## Peak vorticity (1.5× band)

![Peak vorticity from real saved fields](wmax_real.png)

| Case | Source | W₀ | Final W | First arrival | Sustained (0.02) |
|------|--------|----|---------|---------------|------------------|
| aligned | fd4-n112-base | 80 | 135.35 | **0.36** | **0.36** |
| compressive | fd4-n112-base | 80 | 111.80 | null | null |

Aligned crosses the 1.5× threshold at t = 0.36 and remains above it. Compressive ends below the threshold.

## Switching identity

![Marker index holding the global peak](switching_identity.png)

Global-max marker identity changes on nearly every recorded sample across the six Fourier runs.

| Run | Sequence | Changes |
|-----|----------|---------|
| aligned n112 base | 32 → 18 → 53 → 59 → 62 | 4 |
| aligned n112 half | 32 → 46 → 53 → 5 → 62 | 4 |
| aligned n160 base | 32 → 46 → 53 → 59 → 3 | 4 |
| compressive n112 base | 32 → 15 → 6 → 0 → 1 | 4 |
| compressive n112 half | 32 → 15 → 6 → 0 → 1 | 4 |
| compressive n160 base | 32 → 49 → 6 → 0 → 0 | 3 |

## Files

- [real_data_diagnostics.json](real_data_diagnostics.json) — numbers
- [PROTOCOL.md](PROTOCOL.md) — locked definitions
- [CONTROLS.md](CONTROLS.md) — required gates
- [report.html](report.html) — offline report

All values come from published saved fields. No synthetic series.
