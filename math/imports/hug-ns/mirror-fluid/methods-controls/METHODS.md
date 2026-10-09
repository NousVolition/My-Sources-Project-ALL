# Measurements and controls

The reports distinguish **what was measured**, **which controls passed**, and **what is still running**. A reproducible number can still come from an unresolved grid.

## Three different experiments

| Experiment | What changes the field after time zero? | What the result checks |
| --- | --- | --- |
| Initial asymmetry | Unforced Navier–Stokes evolution | Response to a supplied initial difference |
| Mirror averaging | At the stated time, replace `u` by `(u + M[u])/2` | An explicit deletion of the mirror-odd part |
| Breathing rerun | Before every fluid step, add `0.15 sin(2πt/0.4) breath_u` | Response to prescribed velocity kicks |

The breathing increment has no `dt` factor. Changing its step count changes its driving. Its scalar coefficient is not a measured geometric opening. Mirror averaging deliberately sets the asymmetry to zero; it does not measure spontaneous recovery.

## Definitions

On an `N³` uniform grid in the cube of side `L=6`, use

```text
<a,b> = (dx³/L³) sum_grid sum_components a_i b_i = mean_grid(a · b)
||a|| = sqrt(<a,a>)
M[u](x,y,z) = (-u_x(-x,y,z), u_y(-x,y,z), u_z(-x,y,z))
index of reflected x sample i = (-i) modulo N
a = u - M[u]
U = ||u||; A = ||a||
E = A/U
D = <a,phi>, with ||phi|| = 1
B = ||a - D phi||
A² = D² + B²
```

`E` is dimensionless and unsigned. `D` has the velocity unit used in the simulation and retains a direction. `U`, `A`, and `B` use the same velocity unit. No conversion to physical metres or seconds is specified by the supplied model.

For the breathing, supplied-source and dashboard cases, construct `C` by projecting the Gaussian velocity `(0, 0.25 exp[-8((x−0.55)²+(y−0.4)²+z²)], 0)`. Form `phi=(C−M[C])/2`, then divide by its norm. The template is divergence-free, mirror-odd, and fixed throughout each run. **Positive means alignment with this declared template; negative means its mirror.** It does not label a person as positive or negative.

The 20-run calibration used a different documented Gaussian template and the half-difference `(u−M[u])/2`; its signed coordinate therefore has a factor of two relative to `D` here. The dynamic-q report defines `q_fluid=<u,phi>/U_AB(0)=D/[2 U_AB(0)]`. Its earlier `q_record` follows the imposed diagnostic equation `dq/dt=−0.2q+0.8D`. These are distinct measurements, not interchangeable labels.

`E` cannot generally be recovered from `|D|`: the remainder `B` can grow while the projection falls. In the verified driven small-lean run, 38.24% of the final squared mirror difference lies outside its original template.

For translated or rotated tests, transform the whole field, the mirror plane, and the template together. The new controls implement this by pulling the field back to the original coordinates before measuring `E` and `D`.

## Fluid quantities

| Quantity | Definition | Simulation units |
| --- | --- | --- |
| Speed | `|u|` | length/time |
| Peak spin `W` | `max_grid |ω|`, `ω=curl u` | 1/time |
| Accumulated peak spin `I` | Time integral of `W`; trapezoidal quadrature at every step in the corrected dashboard | dimensionless |
| Kinetic energy `K` | `0.5 dx³ Σ |u|²`, density set to one | length⁵/time² |
| Enstrophy `Z` | `0.5 dx³ Σ |ω|²` | length³/time² |
| Spin ratio | `max|ω| / mean_grid|ω|` | dimensionless |
| Fixed-threshold volume | `dx³ count(|ω| ≥ c)` for `c=50,100,200,400` | length³ |
| Half-peak volume | `dx³ count(|ω| ≥ W/2)` | length³ |

The old uploaded dashboard used right-endpoint rectangles for `I`. Both the original and the corrected values are retained under their stated definitions.

The equivalent-sphere diameter derived from half-peak volume is an aggregate size, **not a tube radius or minimum core thickness**. Several disconnected peaks can contribute to it. Physical narrowing needs spatial structure and grid checks as well as this number. The fixed thresholds can all be empty in low-spin runs; that is reported as zero volume, not evidence of adequate resolution.

## Terms at one location

At every saved output, select the current global vorticity maximum. Save its index, position, change of index, and shortest periodic displacement from the previous selected maximum. This is a sequence of global maxima; it is not tracking the same fluid parcel or proving that the same vortex remains dominant.

At that **same location**, record two separate budgets:

```text
Velocity:  −(u·∇)u  −∇p  +ν Δu
Vorticity: −(u·∇)ω  +(ω·∇)u  +ν Δω
```

Save each vector and its signed projection along the local velocity or vorticity, respectively. A positive projection increases that local magnitude at that instant; a negative projection decreases it. Pressure contributes to the velocity budget; its curl vanishes and it is not a separate term in this constant-density vorticity equation. Terms evaluated at other maxima cannot be combined into this budget.

The corrected dashboard also compares the vorticity terms with the curl of the actual discrete velocity RHS. Their difference is reported as a discrete product-rule residual; it is not silently assigned to physical stretching.

## Energy, enstrophy and forcing

For unforced periodic flow, the continuum energy balance is `K(t)−K(0)+2ν∫Z dt=0`. Report its numerical residual relative to `K(0)`. Also save the measured nonlinear energy rate; do not assume the discrete advection contribution cancels exactly.

For each breathing step, first measure the exact energy and enstrophy changes caused by its kick. Then integrate the fluid contribution between the kicked state and the saved next state. The reported balance subtracts accumulated kick input. Omitting that input would mislabel prescribed driving as a numerical energy error.

Enstrophy rates use the actual curl of the discrete velocity RHS. The physical stretching/advection/viscosity decomposition is recorded separately. Time quadrature has its own error and must be checked under timestep refinement.

## Numerics and spectra

The supplied `clay_hug.py` uses Fourier differentiation and two-stage explicit Heun stepping. Each RHS first projects velocity to divergence-free Fourier coefficients, forms advection on the physical grid, filters its transform, adds viscosity, and projects the RHS. The final Heun field is projected again. Viscosity is `0.01` in these mirror controls. The matched-stretch parameter matrix has its separately recorded viscosities.

NumPy's forward FFT is unnormalized; its inverse divides by `N³`. Modal integrated energy is `0.5 L³ Σ_components |û|²/N⁶`; multiply by `|k|²` for modal enstrophy. All three velocity components contribute.

The supplied mask retains component indices satisfying `|j_x|,|j_y|,|j_z| ≤ floor(N/3)` as evaluated by its existing code. This inclusive rule is retained and audited; it should not be silently replaced by a stricter mask. Boundary-mode aliasing must be tested for grids divisible by three. Passing a symmetry check on particular fields does not prove that the mask removes every possible quadratic alias.

The reported high band is the **retained** set where at least one component index reaches 80% of that cutoff. Shell plots sum modes by rounded radial index and mark the component cutoff and outer corner of retained support. The mask is cubical, so there is no single spherical cutoff.

At every full diagnostic output, record maximum and RMS divergence, relative projection residual, energy and enstrophy spectra, and high-band fractions. The original uploaded band `k²>0.6 max_grid(k²)` contains no retained modes on its 33³ grid. Its reported zero could not diagnose resolution.

## Numerical controls and thresholds

The separate [control protocol](protocol.json) fixes 22 runs before execution: even and opposite starts on 33³, 49³ and 65³; 65³ half-step repeats; two additional equal-amplitude odd shapes; and translated and rotated pairs. It retains the original equation and solver. It does not replace the separate 64³–128³–256³ matched-stretch matrix.

The new runs end at `0.4`. The base step count is the ceiling of `0.4/[0.04 dx/max|u₀|]`; the actual timestep is `0.4/steps`. Paired and geometric controls share the reference timestep. Half-step runs double the count. `W`, `I`, energy and enstrophy are sampled every step; full diagnostics and fields are saved at approximately 0.02 intervals. Every result records its actual times, grid, viscosity, step count and timestep.

Mirror tests compare `F_dt(Mu)` with `M(F_dt(u))` for one complete Heun step, on all three grids at base and half dt. They test both an even and an asymmetric input. A preserved even field alone is insufficient to establish equivariance.

Use `1e-10` as a conservative relative roundoff ceiling for these double-precision mirror checks. Also report the measured residual. This numerical audit threshold is not a continuum approximation error. Signed projection residuals are scaled by the initial velocity norm when the control itself is zero.

Curve comparisons use relative L2 error over their common available interval, with interpolation onto shared times. Near-zero signed curves use the baseline velocity norm as the reference instead of dividing by roundoff. The predefined practical curve threshold is 1%, alongside decreasing refinement error and spectral checks. High-band energy above 0.1% or enstrophy above 1% raises a resolution flag. These are numerical screening criteria, not a theorem.

The original 48/64/80 replay has a common available interval only through `0.12`. The matched-stretch and mirror studies use different starting fields. A common available interval is labelled **resolved** only if its grid, timestep and spectral checks pass together.

## Initial-source comparison

The existing no-C, C-only, reflected-C, 50/50 C–D and 70/30 C–D starts are remeasured for energy, enstrophy, spectra and circulation. Their energies are **measured rather than assumed equal**. Local circulation is defined on specified counterclockwise rectangular loops in the plane `z=0`, using the retained Fourier interpolant. Zero mean vorticity over a periodic box is not a substitute for a local circulation measurement.

The new shape controls hold the initial odd-field L2 amplitude and total energy fixed. Enstrophy and spectral shape can still differ; those differences are reported. Each shape has a declared, fixed signed template, and its opposite run uses the full mirror of its initial velocity.

## Marble comparison

The supplied ball code retains both position and velocity and includes damping, so overshoot is possible. Its single-well and double-well forces are explicitly chosen. The coefficient multiplying the linear term controls whether the center restores or destabilizes a lean; the cubic term specifies saturation. A rotating-hoop model would require a separate rotation/gravity parameter and a stated approximation before an overdamped one-coordinate reduction could be used. No such hoop parameter has been fitted to the fluid.

The useful conceptual correspondence is signed position with signed projection, and distance from symmetry with an unsigned mismatch. It is not an exact numerical mapping `|x|=E`: the fluid has the additional remainder `B` and a changing normalization `U`. The ball's two stable wells do not demonstrate two fluid attractors.

## Reporting categories

| Category | Evidence required |
| --- | --- |
| Verified mirror equivariance | One-step commutator residuals, including asymmetric inputs, below the declared numerical tolerance |
| Mirror-paired response | Opposite runs have equal `E`, opposite `D`, and reflected full fields within the declared tolerance |
| Maintained asymmetry | An imposed initial difference remains measurable during the stated finite window |
| Persistent asymmetric state | Converged long-time behavior without repeated driving or manual resets; not established here |
| Resolved local growth | Local growth with agreeing grid/time curves and adequately controlled spectra and structure |
| Unresolved local growth | Growth with material grid dependence or concentration at the cutoff; the current matched-stretch case belongs here |

For a future **finite-window plateau candidate**, the predefined observation window is `[T/2,T]` with length at least one simulation time unit; `|D|` must exceed 100 times the measured numerical symmetry floor, vary by at most 1%, and agree under three grids and a halved timestep. No run in this 0.4 control batch can meet that duration. Even a plateau pass would not by itself establish an attractor.

Use **unfinished** for jobs still running, **numerically unresolved** for completed calculations that fail numerical controls, and **finite over the tested window** for the observed finite values. These labels answer different questions.
