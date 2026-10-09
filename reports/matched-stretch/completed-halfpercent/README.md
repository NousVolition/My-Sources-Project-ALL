# Three completed 0.5% perturbation runs

All three 64-grid runs reached time 0.40. Their saved final fields were independently remeasured before publication. The original supplied field, viscosity 0.001 and zero force were retained.

| Seed | Largest saved spin | Accumulated peak | Final relative separation | Endpoint log growth rate |
| --- | ---: | ---: | ---: | ---: |
| 101 | 681.66294 | 166.72859 | 0.797910 | 12.68139 |
| 202 | 715.33026 | 166.63841 | 0.804008 | 12.70043 |
| 303 | 756.43654 | 168.09001 | 0.799290 | 12.68571 |

![Three-seed comparison](seed-comparison.png)

Relative separation is the L2 difference from the baseline, divided by the initial baseline L2 norm. Each initial separation is 0.005. The endpoint log growth rate is log(D(0.4)/D(0))/0.4. It describes a finite window and finite perturbation.

**Spatial resolution fails:** about 67% of endpoint enstrophy is in the upper retained band. This prevents a physical sensitivity or asymptotic Lyapunov conclusion. The original raw-strain periodic-join limitation remains.

[Direct field verification](verification.json) · [Main controls](../README.md) · [Original protocol](../protocol.json)
