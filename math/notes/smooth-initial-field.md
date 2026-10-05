# A smooth initial velocity field with swirl

This construction is smooth through the axis, divergence-free everywhere, and rapidly decreasing. Its energy can be calculated exactly and adjusted to a chosen value.

The [editable LaTeX note](smooth-initial-field.tex) contains the full derivation. [The reference implementation](../initial_field.py) uses the Cartesian formula directly.

## 1. Define the shape and the flow

Write $q=x^2+y^2$ and $r=\sqrt q$. Choose $R\geq0$, $a,b>0$, and real $A,B$.

| Parameter | Meaning |
| --- | --- |
| $R$ | Radius of the ring where the envelope peaks |
| $a$ | Radial concentration control |
| $b$ | Vertical width |
| $A$ | Strength and direction of meridional circulation |
| $B$ | Strength and direction of swirl |

If coordinates have units of length, $R,a,b$ are lengths and $A,B$ are inverse times. The examples use nondimensional variables.

Define the envelope, stream function, and swirl:

$$
\phi=\exp\left[-\frac{(r^2-R^2)^2}{a^4}-\frac{z^2}{b^2}\right],
\qquad \psi=A r^2z\phi,\qquad u_\theta=B r\phi.
$$

Differentiating the stream function gives, for $r>0$,

$$
u_r=-A r\left(1-\frac{2z^2}{b^2}\right)\phi,\qquad
u_z=A z\left(2-\frac{4r^2(r^2-R^2)}{a^4}\right)\phi.
$$

For evaluation on all of space, including the axis, use the Cartesian expression. Set

$$
H=1-\frac{2z^2}{b^2},\qquad C=2-\frac{4q(q-R^2)}{a^4}.
$$

Then

$$
\boxed{\mathbf u_0=
\bigl(-A xH-B y,\;-A yH+B x,\;A zC\bigr)\phi.}
$$

There is no division by $r$. On the axis,

$$
\mathbf u_0(0,0,z)=
\left(0,0,2Az\,e^{-R^4/a^4-z^2/b^2}\right).
$$

The previous correction is recovered with $a=b=\alpha$ and $A=B=1$.

## 2. Why it is smooth and incompressible

Every component is a polynomial in $x,y,z$ times the exponential of a polynomial. Hence the field is real analytic, including at $x=y=0$.

For decay, the inequality

$$
(q-R^2)^2\geq q^2/2-R^4
$$

bounds the envelope by

$$
\phi\leq e^{R^4/a^4}\exp\left[-\frac{q^2}{2a^4}-\frac{z^2}{b^2}\right].
$$

Each derivative is another polynomial times this rapidly decreasing envelope. Thus $\mathbf u_0$ belongs to the Schwartz space: all derivatives decay faster than every inverse power of distance. In particular, kinetic energy and the squared integrals of all derivatives are finite.

Incompressibility can be checked directly in Cartesian coordinates:

$$
\partial_xu_x+\partial_yu_y=-AHC\phi,\qquad
\partial_zu_z=AHC\phi.
$$

Their sum is exactly zero. This calculation also holds on the axis.

## 3. Exact energy and a controllable example

Define energy per unit constant density by

$$
E_0=\frac12\int_{\mathbb R^3}|\mathbf u_0|^2\,dV.
$$

With

$$
Z_0=b\sqrt{\pi/2},
$$

$$
M_0=\frac{a^2\sqrt\pi}{2\sqrt2}
\left[1+\operatorname{erf}\left(\frac{\sqrt2 R^2}{a^2}\right)\right],
\qquad
M_1=R^2M_0+\frac{a^4}{4}e^{-2R^4/a^4},
$$

the exact result is

$$
\boxed{
E_0=\frac{\pi Z_0}{2}
\left[
\left(B^2+\frac{3A^2}{4}+\frac{A^2b^2R^2}{a^4}\right)M_1
+\frac{3A^2b^2}{4}M_0
\right].
}
$$

Here $\operatorname{erf}$ is the error function. The derivation integrates the squared velocity over the angle and height, then reduces the remaining radial integral to $M_0$ and $M_1$; details are in the LaTeX note.

For $R=1.5$, $a=b=1$, and $A=B=1$:

$$
E_0=24.05715786588763.
$$

To give a nonzero base field a target energy $E_*>0$, multiply both $A$ and $B$ by $\sqrt{E_*/E_0}$. This preserves the ratio of swirl to meridional flow. For the example, choosing

$$
A=B=0.20388150977663805
$$

gives unit energy.

## 4. The axis correction changes the shape too

The old envelope uses $(r-R)^2$, with $r=\sqrt{x^2+y^2}$. Along $y=0$, that becomes $(|x|-R)^2$. For $R>0$ its left and right derivatives at the axis disagree.

This also affects the original velocity. At fixed $z$, its axial component has the expansion

$$
u_{z,\mathrm{old}}(x,0,z)=
Az\phi_{\mathrm{old}}(0,z)
\left[2+\frac{6R}{\alpha^2}|x|+O(x^2)\right].
$$

It has a cusp when $AzR\ne0$. The corrected field uses a polynomial in $r^2$.

![Comparison of the original radial Gaussian and corrected smooth ring envelope](../figures/smooth-initial-profile.png)

The new radial width near $r=R>0$ is approximately $a^2/(2R)$, in the convention $\exp[-((r-R)/w)^2]$. Equal scale values therefore do not give equal widths. To match the old width $\alpha$ to quadratic order at the ring, choose $a=\sqrt{2R\alpha}$ and $b=\alpha$. This matches local curvature at the peak, not the whole original shape.

## 5. What has been checked

Seven code tests cover axis values, stream-function differentiation, numerical Cartesian divergence, rotational symmetry, energy scaling, the zero field, and invalid geometric parameters. They pass using the Python standard library.

Direct numerical integration of the squared Cartesian velocity agrees with the exact energy formula to relative error below $1.3\times10^{-14}$ in four parameter choices, including zero ring radius and unequal widths. The calculation used a 180-point Gauss-Legendre rule in both radius and height, on a large truncated domain.

These calculations check the implementation. The analytic identities above establish smoothness, incompressibility, decay, and finite energy for the full parameter family.

If a scalar initial field is needed, $S_0=S_*\phi$ is also Schwartz. Smoothness of a force containing $|\nabla S|$ at critical points remains a separate issue. The editable note now examines a discrete projection and short box experiments. Forcing, boundary treatment and convergence remain separate questions.

For the precise initial-data and forcing requirements in the standard Navier-Stokes problem, see [Fefferman's official Clay Mathematics Institute problem description](https://www.claymath.org/library/monographs/MPPc.pdf).

## 6. Discrete follow-up

The existing LaTeX note includes the compatible centered Poisson stencil, the missing-spacing and oscillation defects in the submitted iteration, and measurements from length-6 and length-18 boxes. For the original bump, the length-6 raw divergence maxima are 0.46351, 0.91375, and 1.71323 at 8, 16, and 32 points. The finer-grid maxima lie on the periodic boundary. The corrected projection reduces measured divergence to floating-point roundoff, without proving that the boundary model or evolution is accurate. The smooth comparison starts both branches from the same projected field.
