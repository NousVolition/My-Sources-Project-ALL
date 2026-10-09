# Test design

| Stage | Measurement | Control |
| --- | --- | --- |
| Amplitude selection | 3 μ values × 3 initial amplitudes; positive-peak events through μt=16. | Compare analytic averaged envelope; last peaks within 0.01 of 2; spread below 0.002. |
| Fast–slow separation | μ=5,10,20,40; final full cycle after t up to 12μ+40. | Refine Radau tolerance; compare DOP853 at μ=10; use fixed traversal thresholds. |
| Conservation and damping | Duffing ε=0.1; cubic velocity damping ε=0.1,0.2,2. | Energy quadrature or accumulated loss; no assumed agreement outside small ε. |
| Small-angle prediction | Pendulum a=0.1,0.2,0.4. | Exact elliptic period; omitted frequency term scales as a⁴. |
| Parametric growth | Swing ε=0.1,0.05,0.025; γ=−0.75,−0.25,0,0.25,0.75. | One-period linear map; determinant=1; nonlinear zero/seed controls. |
| Bifurcation mechanism | Explicit pitchfork and supercritical Hopf normal forms. | Eigenvalues, determinant, exact radius evolution; do not infer missing equations. |

See data/protocol.json and the saved exercise/bifurcation results for exact parameters, formulas and acceptance criteria. These are mathematical benchmark tests, not water or biological experiments.
