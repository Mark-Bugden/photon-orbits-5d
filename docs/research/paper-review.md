# Review of Bugden (2018), "Spherical photon orbits around a 5D Myers-Perry black hole"

*arXiv:1801.03389 (v1, 10 Jan 2018, the only arXiv version); Gen. Rel. Grav. 50 (2018), doi:10.1007/s10714-018-2352-x.*
*This review covers the arXiv v1 LaTeX source and PDF. The journal version is paywalled, so I did not check it. Page numbers below refer to the arXiv PDF.*

The review has two jobs:

- **Part A** is an equation-by-equation errata check, with symbolic and numerical verification.
- **Part B** is everything an implementer needs for the interactive 5D demo: the parameter ranges, the solver for the controls, equations of motion suited to JavaScript, the projection, presets for the paper's figures, the extra time-frequency maths for resonances, and unit-test values.

All checks were done in Python (run with `uv`: sympy, numpy, scipy, mpmath, matplotlib, pillow, pymupdf). The scripts are in the session scratch folder, not in the repo, and are listed in [Appendix C](#appendix-c-verification-scripts).

---

## 0. Executive summary

**Verdict.** The core physics of the paper is correct. I re-derived the following independently and they agree exactly:

- the inverse metric (3.3);
- the Hamilton-Jacobi separation;
- the first-order equations (2.3);
- the Hamiltonian formulation (3.4)–(3.8);
- the cubic coefficients $a_0,\dots,a_3$;
- the Fourier-Budan analysis;
- the equatorial and polar radii (3.9)–(3.15).

All the figure parameters reproduce the figures. For Figs 4.9 and 4.11 the curves match to pixel accuracy, about 0.005 rad.

**Confirmed errors.** None of these affects the conclusions.

1. **Eq. (2.1):** the last term reads $(x+b)^2\cos^2\theta\,d\psi^2$. It should be $(x+b^2)\cos^2\theta\,d\psi^2$.
2. **Mass formula (p. 2):** the paper gives $M = 3\mu/(8\sqrt{\pi}G)$. It should be $M = 3\pi\mu/(8G)$.
3. **§3.1.1, $a=b$ (p. 4):** the maximum of $\Theta$ is $a^2+K-(|\Phi|+|\Psi|)^2$, not $a^2+K-(\Phi+\Psi)^2$. It sits at $\tan^2\theta = |\Phi|/|\Psi|$. The paper's criterion $K\ge(\Phi+\Psi)^2-a^2$ is wrong whenever $\Phi\Psi<0$; there is a counterexample in §A4.
4. **§3.4, the display of $\chi$ for equatorial orbits (p. 7):** the constant term is printed as a copy of the $x^1$ coefficient. It should be $b^2\big[(\Phi+a)^2 + b^2(a^2-\Phi^2)\big]$. The solutions (3.9)–(3.12) are nevertheless correct, because they solve the correct cubic and not the printed one.
5. **Figs 4.9 and 4.11:** the curve labelled "x" actually plots $r=\sqrt{x}$. It sits at 1.378 and 1.126, whereas $x$ is 1.899 and 1.267.

There is also a short list of imprecise statements (horizon condition, sign of $\Omega_{a,b}$, the Kerr-Schild remark), notational slips ($\mathcal X$ vs $\chi$, the $\theta$ range $[0,\pi]$ vs $[0,\pi/2]$) and typos, including two in the references. See §A4.

**Key results for the demo (Part B).**

- For fixed non-extremal $(a,b)$, the pair $(\Phi/E,\Psi/E)$ always determines a **unique** spherical photon orbit $(x,K)$, when one exists. Equivalently, $V(x)=P(x)/\Delta(x)$ has exactly one critical point outside the horizon. I found no exception in $2\times10^5$ random samples plus dense grids.
- That orbit is always radially unstable.
- The set of $(\Phi,\Psi)$ that admit a spherical orbit is a **compact convex "curved diamond"** containing the origin. Its four vertices lie on the axes, at the paper's equatorial and polar circular orbits. Its edges are the constant-$\theta$ circular orbits, as in Fig 4.12.
- **$(\Phi,\Psi)$ are therefore good controls.** Store them internally as (direction, fraction of the way to the boundary) so they survive changes of $a$ and $b$.
- **Equations of motion:** use Mino time and an angle variable $\chi$ with $\sin^2\theta = \bar y - \Delta y\cos\chi$. With this choice:
  - there are no turning points and no sign flips;
  - the axis spikes are handled analytically;
  - one $\theta$-oscillation is a fixed-step Simpson quadrature with about 256 steps;
  - any length of history follows by periodic extension.
- **Do not integrate the radial equation.** A full $(x,\theta,p_x,p_\theta)$ integration falls off the orbit after about 6 oscillations, even at rtol $10^{-12}$.
- **Projection:** use the paper's projection, $(x_1,y_1,z_1)$ with a cyclic colour for $\psi$. It is lossless and it matches the figures.

---

# Part A — careful read and errata check

## A1. Summary of the paper

### Setup (§2)

The metric is the 5D Myers-Perry metric in Boyer-Lindquist-type coordinates $(t,x,\theta,\phi,\psi)$, with $x=r^2$, mass parameter $\mu$ and rotation parameters $a$ (in the $\phi$-plane) and $b$ (in the $\psi$-plane). With the typo in (2.1) corrected, it is

$$
ds^2=-dt^2+\frac{\mu}{\rho^2}\big(dt+a\sin^2\theta\,d\phi+b\cos^2\theta\,d\psi\big)^2+\frac{\rho^2}{4\Delta}dx^2+\rho^2d\theta^2+(x+a^2)\sin^2\theta\,d\phi^2+(x+b^2)\cos^2\theta\,d\psi^2,
$$

where

$$
\rho^2=x+a^2\cos^2\theta+b^2\sin^2\theta,\qquad \Delta=(x+a^2)(x+b^2)-\mu x .
$$

This is exactly the form used by Frolov & Stojković (FS03, gr-qc/0301016), by Kagramanova & Reimers (KR12, arXiv:1208.3686) and by Diemer et al. (DKLR14, arXiv:1404.3865). In those papers $\theta\in[0,\pi/2]$ and $\phi,\psi\in[0,2\pi)$.

**Horizons.**
$x_\pm=\tfrac12\big(\mu-a^2-b^2\pm\sqrt{(\mu-a^2-b^2)^2-4a^2b^2}\big)$. The exterior is $x>x_+$.

**Discrete symmetry (2.2).**
$a\leftrightarrow b$, $\theta\leftrightarrow\pi/2-\theta$, $\phi\leftrightarrow\psi$, and therefore $\Phi\leftrightarrow\Psi$.

**Conserved quantities.**
$p_t=-E$, $p_\phi=\Phi$, $p_\psi=\Psi$, plus a Carter-like constant $K$, which comes from separating the Hamilton-Jacobi equation. The paper normalises by $E$: in effect $E=1$, and $\Phi,\Psi,K$ stand for $\Phi/E$, $\Psi/E$, $K/E^2$.

**First-order equations (2.3)** (from FS03). With
$\mathcal E = E+\frac{a\Phi}{x+a^2}+\frac{b\Psi}{x+b^2}$:

$$
\begin{aligned}
\rho^2\dot t&=E\rho^2+\frac{\mu(x+a^2)(x+b^2)}{\Delta}\mathcal E,\\
(\rho^2\dot x)^2&=4\mathcal X,\qquad \mathcal X=\Delta\Big[xE^2+(a^2-b^2)\Big(\tfrac{\Phi^2}{x+a^2}-\tfrac{\Psi^2}{x+b^2}\Big)-K\Big]+\mu(x+a^2)(x+b^2)\mathcal E^2,\\
(\rho^2\dot\theta)^2&=\Theta=E^2(a^2\cos^2\theta+b^2\sin^2\theta)+K-\frac{\Phi^2}{\sin^2\theta}-\frac{\Psi^2}{\cos^2\theta},\\
\rho^2\dot\phi&=\frac{\Phi}{\sin^2\theta}-\frac{\mu a(x+b^2)}{\Delta}\mathcal E-\frac{(a^2-b^2)\Phi}{x+a^2},\\
\rho^2\dot\psi&=\frac{\Psi}{\cos^2\theta}-\frac{\mu b(x+a^2)}{\Delta}\mathcal E+\frac{(a^2-b^2)\Psi}{x+b^2}.
\end{aligned}
$$

### Conditions for spherical orbits (§3)

**§3.1, the $\theta$-equation.**

- Nonzero $\Phi$ excludes the $\phi$-axis ($\theta=0$). Nonzero $\Psi$ excludes $\theta=\pi/2$.
- For $a=b$ the paper finds the maximum of $\Theta$ (it contains error 3 above).
- For $a\neq b$, substituting $y=\sin^2\theta$ gives $\Theta=Q(y)/[y(1-y)]$, with $Q$ cubic. The Fourier-Budan theorem then yields necessary conditions on $K$: (3.2), plus its $a^2<b^2$ analogue.

**§3.2, the radial equation.**

- $\mathcal X$ is a cubic in $x$, written $\chi=x^3+a_2x^2+a_1x+a_0$.
- Spherical orbits need $\chi=\chi'=0$ at $x>x_+$. The paper suggests solving for $(x,K)$ with Mathematica or Maple and checking with the cubic discriminant.
- Every example in the paper has $\chi''>0$, so the orbits are unstable to radial perturbations.

**§3.3.** The equations are rewritten in Hamiltonian form, (3.3)–(3.8), so that the integrator has no turning-point sign flips. This follows Liu, Ding & Jing (2017).

**§3.4, explicit solutions.** These are the equatorial ($\theta=\pi/2$, $\Psi=0$, $K=\Phi^2-b^2$) circular photon orbits with $\mu=1$:

$$
x=1\pm a-b^2+\sqrt{(1\pm a)^2-b^2},\qquad \Phi=\pm\Big(1+\sqrt{(1\pm a)^2-b^2}\Big).
$$

- For $b=0$ this reduces to $x=2(1\pm a)$, in agreement with Cardoso et al. 2009 (CMBWZ09): $r_c^2=2(\mu\pm a\sqrt\mu)$.
- The polar ($\theta=0$) orbits follow from the symmetry: $x=1\pm b-a^2+\sqrt{(1\pm b)^2-a^2}$.

### Special cases

| Case | Content |
|---|---|
| $b=0$ (§4.2.1) | Zero-momentum orbit at $x=2\mu-a^2=x_++\mu$. |
| $a=b$ | The $\theta$-equation is analysed in §3.1.1. |
| $a=b=0$ (Schwarzschild-Tangherlini) | Not discussed explicitly. It is the $a,b\to0$ limit: $x=2\mu$, $K=4\mu$. |

### Stability, quasiperiodicity, closed orbits

- All the examples are radially unstable.
- The abstract and Fig 4.3 mention "quasiperiodic behaviour".
- There is **no** discussion of frequencies, resonances or closed orbits.

### How the figures were made (§4.1)

The figures were evidently made in Mathematica: default plot styling, `0.` tick labels, `Hue` colours.

**Coordinates.** The paper uses the "pseudo-Cartesian" Myers-Perry coordinates

$$
x_1+iy_1=\sqrt{r^2+a^2}\,\sin\theta\,e^{i[\phi-\arctan(a/r)]},\qquad x_2+iy_2=\sqrt{r^2+b^2}\,\cos\theta\,e^{i[\psi-\arctan(b/r)]}.
$$

It then projects to 3D with $z_1=\sqrt{x_2^2+y_2^2}$, written in the paper as $z_1=\sqrt{r^2+b^2}\cos\theta$.

**What is drawn.**

- The orbit lies on the spheroid $\frac{x_1^2+y_1^2}{r^2+a^2}+\frac{z_1^2}{r^2+b^2}=1$, drawn as an opaque black surface with a white ring at $z_1=0$.
- The fourth dimension is shown by colouring the path by $\psi$ on a hue wheel (Fig 4.1).
- The side panels (Figs 4.9 and 4.11) plot $\theta,\phi,\psi$ (and "x") against the affine parameter $\lambda$.

## A2. Equation-by-equation verification

Script: `check_symbolic.py`, output in `check_symbolic.out`.

**Method.** Rational-function identities were tested by exact evaluation at 8 random rational points. A nonzero rational function vanishes at such points with probability 0. To keep sympy fast, the angular dependence is written with $S=\sin^2\theta$.

| Item (page) | Status | How verified |
|---|---|---|
| (2.1) metric (p. 2) | **Typo** $(x+b)^2\to(x+b^2)$ | Compared with FS03, KR12 and DKLR14. Also, (3.3) is the inverse of the corrected metric. |
| $\rho^2,\Delta$, $x_\pm$ (p. 2) | ✓ | $\Delta(x_+)=0$ symbolically. |
| $M=3\mu/(8\sqrt\pi G)$ (p. 2) | **Error** | Standard Myers-Perry: $M=\frac{(D-2)\Omega_{D-2}}{16\pi G}\mu=\frac{3\pi\mu}{8G}$ with $\Omega_3=2\pi^2$. DKLR14 eq. (1) gives the same. |
| $\Omega_{a}=a/(x_++a^2)$, $\Omega_b$ (p. 2) | ✓ in magnitude, sign convention | The horizon generator for this metric is $\partial_t-\frac{a}{x_++a^2}\partial_\phi-\frac{b}{x_++b^2}\partial_\psi$ (checked: $g(\xi,\xi)=0$ on $x_+$ for all $\theta$; the "+" choice fails). Frame dragging $\omega=-g_{t\phi}/g_{\phi\phi}<0$ for $a>0$. FS03 writes the same formula. See §A4(iii). |
| Horizon condition $\mu\ge(a+b)^2$ (p. 2) | Imprecise | See §A4(ii). |
| (2.2) symmetry (p. 3) | ✓ | $\mathcal X$, $\Theta$, $\Delta$, $\rho^2$ are all invariant. |
| (2.3a–e), $\mathcal X$, $\Theta$, $\mathcal E$ (p. 3) | ✓ | $\rho^2\cdot 2H = 4\Delta p_x^2+p_\theta^2-(\mathcal X+\Delta\Theta)/\Delta$ exactly, with the inverse metric from sympy. $\dot q^\mu=g^{\mu\nu}p_\nu$ reproduces (2.3a,d,e). They also match FS03's first-order equations term by term. |
| $E$-normalisation (p. 3) | ✓ | $\mathcal X$ and $\Theta$ are homogeneous of degree 2 in $(E,\Phi,\Psi,\sqrt K)$. |
| (3.1) (p. 4) | ✓ | |
| §3.1.1, max of $\Theta$ for $a=b$ (p. 4) | **Error** (signs) | Symbolic maximisation gives $-(|\Phi|+|\Psi|)^2$; see §A4(i). |
| $Q(y)$ (p. 4) | ✓ | $y(1-y)\Theta-Q\equiv0$. |
| $\mathcal F_Q(0)$, $\mathcal F_Q(1)$ (p. 4) | ✓ | Derivatives recomputed. |
| Fourier-Budan logic and $K$ constraints, incl. (3.2) (p. 5) | ✓ (necessary, not sufficient, as the paper says) | Logic re-checked; the $a^2>b^2$ case follows by symmetry. |
| Cubic $\chi$, coefficients $a_3..a_0$ (p. 5) | ✓ | Coefficient-by-coefficient match. Also: $\chi=P(x)-K\Delta(x)$, see Part B. |
| Discriminant $D$ (p. 5) | ✓ | Standard cubic discriminant. |
| "$\chi''>0$ ⇒ unstable" (p. 5) | ✓ | Since $\chi\ge0$ on both sides of the double root. |
| (3.3a–h) inverse metric (p. 6) | ✓ | All 8 components, exact. Identical to FS03. |
| (3.4)–(3.7) (p. 6) | ✓ | $f=-(\mathcal X+\Delta\Theta)/(2\rho^2\Delta)$ is $K$-independent, and $H$ in (3.7) equals $\tfrac12g^{\mu\nu}p_\mu p_\nu$ exactly. |
| (3.8a–g) (p. 7) | ✓ | (3.8c–e) equal $\dot q^\mu=g^{\mu\nu}p_\nu$ exactly. |
| §3.4, equatorial $\chi$, $x^2$ and $x^1$ coefficients (p. 7) | ✓ | |
| §3.4, equatorial $\chi$, constant term (p. 7) | **Typo** | Printed as a copy of the $x^1$ coefficient. |
| (3.9)–(3.12) (p. 7) | ✓ | All four branches satisfy $\chi=\chi'=0$ exactly. `sympy.solve` finds exactly 5 solutions, the 5th being $x=-b^2$, $\Phi=-a$, as the paper says. |
| "+ roots $\ge x_+$, − roots $\le x_+$" (p. 8) | ✓ | 14,520 $(a,b,\pm)$ combinations on a grid with $a,b\ge0$, $a+b<1$; no violations. Equality only at extremality. |
| (3.13)–(3.15) (p. 8) | ✓ | (3.15) via symmetry. Also confirmed numerically as the axis intercepts of the allowed region (§B2). |
| §4.1 projection, spheroid footnote (pp. 8–9) | ✓, with remarks | See §B4. Oblate iff $a^2>b^2$ ✓. |
| $x=2\mu-a^2\equiv x_++\mu$ for the zero-momentum orbit, $b=0$ (p. 9) | ✓ | Derived symbolically: double root at $x=2\mu-a^2$ (or $-a^2$), with $K=4\mu-a^2$. |
| Figure parameters and claims (pp. 10–15) | ✓ | See §B5. |
| Fig 4.12: $\Phi\approx-0.9292$, $\Psi\approx0.9036$ (p. 15) | ✓ | Exact: $\Phi=-0.929241$, $\Psi=+0.903667$, $x=1.455284$, $K=3.203956$. The paper's $\Psi$ is truncated rather than rounded, which is trivial. |

**Comparisons with the literature.**

| Comparison | Result |
|---|---|
| FS03 | The metric, inverse metric and first-order equations are identical; the paper copies them correctly. |
| KR12 | They use the same equations, with their $\mathcal E_{KR}=(x+a^2)(x+b^2)\mathcal E$, and Mino time $\rho^2d\gamma=d\tau$, which is what I recommend in §B3. |
| DKLR14 | Same metric and mass formula $M=3\pi\mu/8G$. They also treat $\theta\in[0,\pi/2]$. |
| CMBWZ09 | 5D singly spinning: $r_c^2=2(\mu\pm a\sqrt\mu)$ and $D_c=-a\mp2\sqrt\mu$. These match (3.14) and $\Phi=\pm(2\pm a)$, i.e. $|\Phi|=2+a$ counter-rotating and $2-a$ co-rotating. |
| Teo 2003 (Kerr) | The 5D analogue of Teo's one-parameter family $r\in[r_{ph}^{pro},r_{ph}^{retro}]$ is the two-parameter family of §B2. The qualitative types in Figs 4.2 and 4.3 are Teo's. |

## A3. Independent numerical integration

Scripts: `integrate_check.py`, output in `integrate_check.out`; also `lyapunov.py`.

**Method.**

- Hamilton's equations were built directly from the sympy inverse of the corrected metric. This deliberately does **not** use the paper's $\mathcal X$ or $\Theta$.
- Each run starts on the predicted spherical orbit at a $\theta$ turning point, with $p_x=0$ and $p_\theta=0$.
- Integrator: `solve_ivp(DOP853, rtol=1e-12, atol=1e-13)` over 3 $\theta$-oscillations.
- The result is compared with the Mino-time angle-variable integration recommended in §B3.
- $\theta$ is folded into $[0,\pi/2]$ for the comparison, because the BL integration runs through the axes into the paper's $[0,\pi]$ range.

The $\theta$ bands are in degrees, as predicted from the roots of $Q$ and as reached by the Boyer-Lindquist integration. The last column is the largest difference between the Hamiltonian and Mino-χ results at equal $\lambda$.

| Case | $a$ | $b$ | $\Phi$ | $\Psi$ | $x_0$ | θ band, predicted | θ band, BL Hamiltonian | max $\lvert x-x_0\rvert$ | max $\lvert H\rvert$ | max diff in $\theta$, $\phi$, $\psi$, $t$ |
|---|---|---|---|---|---|---|---|---|---|---|
| Tangherlini | 0 | 0 | 0.7 | 0.9 | 2.000000 | 23.6955–60.5654 | 23.6955–60.5653 | 7e-9 | 3e-12 | 4e-8, 7e-8, 5e-8, 5e-9 |
| Fig 4.3 | 0.5 | 0 | 2 | 0 | 2.750000 | 55.9375–90 | 55.9375–89.962 | 1e-8 | 3e-12 | 2e-8, 1e-8, 0, 4e-9 |
| Fig 4.7 | 0.6 | 0.3 | 1 | −1 | 1.754557 | 37.5908–51.5256 | 37.5908–51.5256 | 3e-9 | 2e-12 | 1e-8, 1e-8, 2e-8, 2e-9 |
| Fig 4.8 (extremal) | 0.9 | 0.1 | 0.8 | 0.8 | 1.899039 | 22.2319–66.0957 | 22.2319–66.0956 | 2e-8 | 1e-11 | 6e-8, 9e-8, 9e-8, 1e-8 |
| Fig 4.10 | 0.9 | 0.09 | 0.1 | 0.01 | 1.267423 | 2.8074–89.6874 | 2.8074–89.676 | 7e-9 | 2e-11 | 4e-6, 7e-6, 5e-4 (see note), 2e-8 |
| generic $a>b$ | 0.55 | 0.2 | −0.6 | 0.4 | 1.406699 | 19.4371–76.2465 | 19.4371–76.2463 | 1.5e-8 | 2e-11 | 1e-7, 1e-7, 3e-7, 2e-8 |
| generic $b>a$ | 0.15 | 0.7 | 0.9 | −0.6 | 1.196022 | 35.5266–66.1852 | 35.5266–66.1852 | 4e-10 | 2e-12 | 4e-8, 3e-8, 6e-8, 5e-9 |
| $a=b$ | 0.4 | 0.4 | 1.2 | −0.5 | 1.967354 | 36.8582–72.6529 | 36.8582–72.6528 | 1.5e-8 | 2e-11 | 5e-8, 3e-8, 1e-7, 9e-9 |
| negative $b$ | 0.5 | −0.3 | 0.8 | 1.0 | 1.692776 | 29.4015–54.2913 | 29.4015–54.2913 | 2e-9 | 3e-12 | 3e-8, 4e-8, 3e-8, 2e-9 |
| 15-oscillation closed orbit (§B6) | 0.6 | 0.3 | 0.1051349 | −0.1354416 | 1.534087 | 3.0511–85.9296 | 3.0511–85.9290 | 2e-8 | 2e-11 | 5e-7, 6e-6, 4e-6, 2e-8 |
| $\Phi=0$ (through the axis) | 0.8 | 0 | 0 | 0.5 | 1.316830 | 0–73.9964 | 0.057–73.9964 | 1e-8 | 1e-11 | 1e-7, 7e-9, 2e-7, 1e-8 |

**Notes on the table.**

- The few-$10^{-4}$ $\psi$ difference for Fig 4.10 is an artefact of comparing at equal $\lambda$ by interpolation: near $\theta_{\max}$, where $\Psi=0.01$, $\psi$ jumps by about $\pi$ over $\Delta\lambda\approx0.02$. The per-period $\Delta\psi$ agrees to $10^{-8}$.
- The shortfalls in the BL bands (89.676 vs 89.687, and 0.057 vs 0) are output sampling at sharp turning points or axis passages.

**Conclusions.**

- $r$ stays constant to about $10^{-8}$ over 3 oscillations.
- The $\theta$ range matches the paper's claims and the roots of $Q(y)$.
- The two formulations agree to $10^{-7}$ or better.
- The Mino period matches the closed form $2K(m)/\sqrt{G_{\max}}$ (§B3) to $10^{-10}$.

**Radial instability caveat** (`lyapunov.py`). Perturbations grow by about $e^{4}$ per $\theta$-oscillation; the growth rate is $\sqrt{2\chi''}$ in Mino time. As a result, a full integration of $(x,\theta,p_x,p_\theta)$ leaves $r=$ const by 0.1%:

- after **3.6–3.9 oscillations at rtol $10^{-8}$**;
- after **5.8–6.4 oscillations at rtol $10^{-12}$**, in every case tested.

Figs 4.3, 4.5 and 4.7 show many more oscillations than that (Fig 4.5 has about 20 loops). So the paper's integration must either have used higher working precision or held $x$ fixed; the paper does not say which. This matters for the demo (§B3).

## A4. List of mistakes and issues

### (i) Confirmed errors

1. **Eq. (2.1), p. 2.**
   - *Says:* $\dots+(x+b)^2\cos^2\theta\,d\psi^2$.
   - *Should be:* $(x+b^2)\cos^2\theta\,d\psi^2$.
   - *Verified:* FS03, KR12 and DKLR14 all have $(x+b^2)$. The paper's own inverse metric (3.3), and every later equation, is consistent only with $(x+b^2)$; I inverted the corrected metric and recovered (3.3) exactly.

2. **Mass, p. 2.**
   - *Says:* $M=\dfrac{3\mu}{8\sqrt\pi\,G}$.
   - *Should be:* $M=\dfrac{3\pi\mu}{8G}$ (5D Myers-Perry with $\Omega_3=2\pi^2$; DKLR14 eq. 1).
   - Not used anywhere else in the paper.

3. **§3.1.1 ($a=b$), p. 4.**
   - *Says:* the maximum of $\Theta$ is at $\theta=\tan^{-1}\sqrt{\Phi/\Psi}$ with $\Theta_{\max}=a^2+K-(\Phi+\Psi)^2$, so a root exists iff $K\ge(\Phi+\Psi)^2-a^2$.
   - *Should be:* the maximum is at $\tan^2\theta=|\Phi|/|\Psi|$, with $\Theta_{\max}=a^2+K-(|\Phi|+|\Psi|)^2$, so a root exists iff $K\ge(|\Phi|+|\Psi|)^2-a^2$.
   - *Verified:* symbolic maximisation. Counterexample: $a=b=0.4$, $\Phi=1.5$, $\Psi=-1.5$. The radial equations give $K=3.6646$, which satisfies the paper's bound ($K\ge-0.16$). But $\max_\theta\Theta=-5.175<0$, so no $\theta$-motion is possible and there is no orbit. The paper itself uses opposite-sign momenta (Fig 4.7, $\Phi=1$, $\Psi=-1$, although that one has $a\neq b$), so the sign restriction is not implicit.

4. **§3.4, display of $\chi$ and $d\chi/dx$ for equatorial orbits, p. 7.**
   - *Says:* the constant term of $\chi$ is $(a^2+\Phi^2+2a\Phi+b^4+2a^2b^2-2\Phi^2b^2)$, identical to the $x^1$ coefficient.
   - *Should be:* $b^2\big[(\Phi+a)^2+b^2(a^2-\Phi^2)\big]=a^2b^2+2ab^2\Phi+b^2\Phi^2+a^2b^4-b^4\Phi^2$.
   - The $x^3$, $x^2$ and $x^1$ coefficients and $d\chi/dx$ are correct.
   - *Verified:* sympy, from the general cubic with $\mu=1$, $\Psi=0$, $K=\Phi^2-b^2$. The printed cubic does **not** vanish at the paper's own solution: for $a=1/2$, $b=1/4$ it gives 7.61. The solutions (3.9)–(3.12) satisfy the corrected cubic exactly, so this is a transcription typo.

5. **Figs 4.9 and 4.11, pp. 14–15.**
   - *Says:* the blue curve is labelled "x".
   - *Actually:* it plots $r=\sqrt x$.
   - *Verified:* the rasters were digitised and calibrated against the tick marks. The blue line sits at $1.378\pm0.001$ (Fig 4.9), against $r=1.37806$ and $x=1.89904$; and at $1.124$ (Fig 4.11), against $r=1.12580$ and $x=1.26742$. In the same plots the $\theta$ curves match my $\theta(\lambda)$ to an rms of 0.004–0.006 rad, so the parameters are certainly those in the captions.

### (ii) Imprecise statements and suspected issues

1. **Horizon condition, p. 2.**
   - *Says:* "This horizon exists and is positive provided $\mu\ge(a+b)^2$".
   - *Should say:* $\mu\ge(|a|+|b|)^2$, since negative $a$ or $b$ are allowed.
   - When $ab=0$ and $\mu=a^2+b^2$ one gets $x_+=0$, which is not positive (the singular extremal limit). §4.2.1 itself correctly requires $a^2<\mu$ strictly.

2. **Kerr-Schild remark, p. 8.**
   - *Says:* "In these coordinates, the metric (2.1) takes the Kerr-Schild form".
   - *Issue:* this is true only if $\phi,\psi$ (and $t$) are the Kerr-Schild angles, which differ from the Boyer-Lindquist ones by $r$-dependent shifts.
   - The sign $\phi-\arctan(a/r)$ is consistent with this metric's sign convention: mapping $a\to-a$ from Kerr's $(r+ia)$ gives $(r-ia)$.
   - For constant-$r$ orbits the difference is a rigid rotation in each plane, so **the pictures are unaffected**. The statement is just loose. Also "h is a function of $\{r,x_i,y_i\}$" is redundant, since $r$ is determined by $x_i,y_i$.

3. **Section 3.1.2 is correct but weak.**
   - The exact criterion for $\theta$-motion is a one-dimensional convex minimisation:
     $K\ge K_\theta(\Phi,\Psi)=\min_{0<y<1}\big[\Phi^2/y+\Psi^2/(1-y)-a^2(1-y)-b^2y\big]$.
   - It satisfies the elementary bounds $(|\Phi|+|\Psi|)^2-\max(a^2,b^2)\le K_\theta\le(|\Phi|+|\Psi|)^2-\min(a^2,b^2)$.
   - It also follows that because $Q(0)=-\Phi^2\le0$ and $Q(1)=-\Psi^2\le0$, a cubic $Q$ has either 0 or 2 roots in $(0,1)$. So **the allowed $\theta$-band is always a single interval**. This is not an error, only a cleaner statement that the implementation uses.

4. **"Cusps" (Fig 4.8, p. 12).**
   - A true cusp in the $(x_1,y_1,z_1)$ projection needs $\dot\theta=\dot\phi=0$ at the same time.
   - For the stated parameters, at $\theta_{\max}$: $d\phi/d\sigma=0.034$, against $d\psi/d\sigma=5.10$.
   - So these are extremely sharp turns rather than exact cusps. The paper's explanation, motion almost entirely in $\psi$, is right.

5. **Stability, p. 5.** The paper could only check stability case by case. In fact every spherical photon orbit is radially unstable: $V=P/\Delta$ has a single critical point, which is its global minimum (§B2). This rests on numerical evidence, not a proof.

### (iii) Notational inconsistencies and ambiguities

1. **$\mathcal X$ versus $\chi$.** §2 defines $\mathcal X$ (`\mathcal X`), but §3 onwards writes $\chi$ ("the equation $\chi(x)$", "zeros of $\Theta$ and $\chi$"). The two are the same function.

2. **Range of $\theta$.**
   - FS03 and DKLR14 use $\theta\in[0,\pi/2]$, which covers $S^3$ once.
   - The paper speaks of "$\theta$ approaches 0 or $\pi$", "$(\pi/2,\pi)$" and "south pole", i.e. a double cover in which $(\pi-\theta,\psi+\pi)$ is the same point as $(\theta,\psi)$.
   - Correspondingly, $z_1$ is defined as $\sqrt{x_2^2+y_2^2}\ge0$ (p. 8) but used as $\sqrt{r^2+b^2}\cos\theta$, which can be negative (p. 9).
   - This is harmless for the figures, and it is in fact what makes the $\Psi=0$ orbits in Figs 4.2, 4.3 and 4.6 pass smoothly through the "equator". It should be stated, though.

3. **Sign of rotation.**
   - With this metric's $+a\sin^2\theta\,d\phi$, a black hole with $a>0$ rotates towards **negative** $\phi$. Co-rotating photons therefore have $\Phi<0$: the smaller radius in (3.11) has $\Phi=-1-\sqrt{(1-a)^2-b^2}$.
   - "$\Omega_a=a/(x_++a^2)$" is the magnitude.
   - The discussion of Fig 4.10 (φ decreasing near the equator because of frame dragging) is consistent with this, but the reader has to infer the sign.

4. **"Quantities related to the mass $M$ and angular velocities", p. 2.** $a$ and $b$ are more naturally related to the angular momenta: $J_a=\tfrac23Ma$ and $J_b=\tfrac23Mb$.

5. **Inconsistent terms.** "Zero momentum orbit", "zero momenta orbit" and "zero-angular-momentum orbit" are used interchangeably. "Lense Thirring" (p. 14) also appears as "Lense-Thirring" (p. 16), and "Myers Perry" (p. 2) as "Myers-Perry".

### (iv) Typos and reference glitches

| Page | Says | Should say |
|---|---|---|
| 1 | "a fixed cordinates radius" | "coordinate radius" |
| 3 | "the physically allowed range for θ occur" | "occurs" |
| 5 | "$a^2>b^2$,we" | missing space |
| 7 | "the radii of such orbits was computed" | "were" |
| 8, 11 | "agreeance" | "agreement" |
| 14 | "intertial frames" | "inertial" |

References:

- **[12]** "Mashroom, B." should be **Mashhoon, B.** The title has "Stabi1ity" (digit 1) for "Stability".
- **[15]** The em-dash is lost in the PDF: "black holesshadows" should be "black holes—shadows".
- **[3], [5], [11], [15]** Old-style category prefixes are attached to new-style identifiers, e.g. "arXiv:hep-th/0812.1806". The standard form is "arXiv:0812.1806 [hep-th]".

### Sections that check out completely

§2 apart from the two items on p. 2, §3.1.2, §3.2, §3.3, the solutions (3.9)–(3.15), §4.2 claims and parameters (every figure reproduced; see §B5), and §5.

---

# Part B — what the implementer needs

## Conventions used throughout Part B

| Quantity | Convention |
|---|---|
| Metric | The corrected metric above |
| Coordinates | $x=r^2$; $\theta\in[0,\pi/2]$; $y\equiv\sin^2\theta\in[0,1]$ |
| Units | **$\mu=1$**: lengths in units of $\sqrt\mu$, times in units of $\sqrt\mu$; $M=3\pi\mu/8G$ |
| Normalisation | **$E=1$**: $\Phi,\Psi$ mean $\Phi/E,\Psi/E$ and $K$ means $K/E^2$. Every spherical orbit found has $dt/d\lambda>0$ (checked on all valid grid points), so $E>0$ and the normalisation loses nothing. |
| Sign of rotation | For $a>0$ ($b>0$) the hole rotates towards $-\phi$ ($-\psi$). Co-rotating photons have $\Phi<0$ ($\Psi<0$). |

## B1. Allowed black holes and the horizon

The horizon is

$$
x_+=\tfrac12\Big(1-a^2-b^2+\sqrt{(1-a^2-b^2)^2-4a^2b^2}\Big),\qquad r_+=\sqrt{x_+}.
$$

It exists iff $|a|+|b|\le1$. Recommended slider domain: $a,b\in[0,0.99]$ with $a+b\le0.99$, with the region clamped dynamically.

**Negative spins.** $a\to-a$ is equivalent to $\phi\to-\phi$ together with $\Phi\to-\Phi$, so negative spins add nothing visually. If you allow them, everything below works with $|a|,|b|$ in the axis-intercept formulas and the signs mirrored.

**Values to avoid or treat carefully.**

- **Exact extremality, $|a|+|b|=1$.**
  - The horizons merge.
  - The co-rotating circular orbits on the axes reach the horizon. For example $x=1-a-b^2+\sqrt{(1-a)^2-b^2}\to x_+$, so a vertex of the allowed region touches the horizon.
  - $\Delta'(x_+)=0$, so the bracket of the radial solver must start strictly above $x_+$.
  - The paper's Fig 4.8 uses $a=0.9$, $b=0.1$, which is exactly extremal. That orbit, at $x=1.899$, is perfectly fine. Allow it if you start the bracket at $x_++10^{-9}$, but keep the slider maximum below 1 so users don't sit at the degenerate corner.
- **$ab=0$ with $|a|$ or $|b|=1$.** Here $x_+=0$, the singular limit with no horizon. Exclude it.
- **$a=b$.** Allowed and harmless in the recommended formulation. It has enhanced $U(2)$ symmetry: $Q(y)$ becomes quadratic and $x,K$ depend only on $\Phi+\Psi$. Avoid any formula that divides by $a^2-b^2$ or uses the third root of $Q$.
- **$a=b=0$ (Tangherlini).** Allowed. Every orbit has $x=2$ and $K=4$, so a radius slider would be meaningless here.

To show the horizon in the same projection, draw the spheroid with semi-axes $\sqrt{x_++a^2}$ (equatorial) and $\sqrt{x_++b^2}$ (polar).

## B2. Parametrising spherical orbits; the controls question

### Reduction to a one-dimensional "potential"

Because $K$ enters $\mathcal X$ only through $-K\Delta$,

$$
\mathcal X(x)=P(x)-K\,\Delta(x),
$$

where, for general $\mu$ ($\mu=1$ in the demo),

$$
\begin{aligned}
P(x)={}&x^3+(a^2+b^2)x^2+\big[(a^2-b^2)(\Phi^2-\Psi^2)+a^2b^2+\mu(a^2+b^2+2a\Phi+2b\Psi)\big]x\\
&+\mu(ab+b\Phi+a\Psi)^2+(a^2-b^2)(b^2\Phi^2-a^2\Psi^2),
\end{aligned}
$$

$$
\Delta(x)=x^2+(a^2+b^2-\mu)x+a^2b^2 .
$$

Both are verified symbolically against the paper's $\mathcal X$; see `check_symbolic.out` and `check_F_grouped.py`.

So $\mathcal X=\Delta\,(V-K)$ with $V(x)=P(x)/\Delta(x)$. Spherical orbits satisfy $\mathcal X=\mathcal X'=0$, which is equivalent to

$$
V'(x_0)=0,\qquad K=V(x_0).
$$

Equivalently, $x_0$ is a root of the quartic $F=P'\Delta-P\Delta'$, which has the sign of $V'$. Written out ($\mu$ general, verified in `check_F_grouped.py`):

$$
\begin{aligned}
F(x)={}&x^4+2(a^2+b^2-\mu)x^3\\
&+\big[(b^2-a^2)(\Phi^2-\Psi^2)-2\mu(a\Phi+b\Psi)+a^4+4a^2b^2+b^4-2\mu(a^2+b^2)\big]x^2\\
&+2\big[b^2(b^2-a^2-\mu)\Phi^2+a^2(a^2-b^2-\mu)\Psi^2-2\mu ab(\Phi\Psi+b\Phi+a\Psi)+a^2b^2(a^2+b^2-\mu)\big]x\\
&+b^2[(b^2-\mu)^2-a^2b^2]\Phi^2+a^2[(a^2-\mu)^2-a^2b^2]\Psi^2+2ab\mu(\mu-a^2-b^2)\Phi\Psi\\
&+2ab^2\mu(\mu-b^2)\Phi+2a^2b\mu(\mu-a^2)\Psi+a^2b^2(a^2b^2+\mu^2).
\end{aligned}
$$

**Why $V$ has a minimum.** At the horizon, $P(x_+)=\mu(x_++a^2)(x_++b^2)\mathcal E(x_+)^2\ge0$ while $\Delta\to0^+$, so $V\to+\infty$ as $x\to x_+$, unless $\mathcal E(x_+)=0$ (the superradiant threshold, non-generic). Also $V\sim x$ as $x\to\infty$. **So $V$ always has at least one minimum outside the horizon.**

### Uniqueness, an answer to "is $(x,K)$ determined by $(\Phi,\Psi)$?"

**Yes.** In every test, $F$ had **exactly one** real root with $x>x_+$:

- `uniqueness_random.py`: $2\times10^5$ random samples with $(a,b)$ uniform in $|a|+|b|<1$, biased towards extremality and including negative spins, and $(\Phi,\Psi)\in[-6,6]^2$, both inside and outside the allowed region. The result was 200,000 out of 200,000.
- `explore_region.py`: 16 black holes, including the extremal $a=b=0.5$ and negative $b$, on a $241^2$ grid of $(\Phi,\Psi)\in[-4.6,4.6]^2$ each. At most one valid orbit per grid point. On every valid orbit $dt/d\sigma\ge2.03>0$, so $E>0$.

That root is therefore the global minimum of $V$. So for every $(\Phi,\Psi)$:

- there is exactly one candidate radius $x_0$;
- $K=V(x_0)$;
- $\mathcal X''=\Delta V''>0$, i.e. **the orbit is radially unstable**. The minimum $\mathcal X''$ over all grids was 0.066, at $a=0.99$.

**Whether the candidate is physical** depends only on the $\theta$ equation:

$$
\text{orbit exists}\iff K=V(x_0)\ \ge\ K_\theta(\Phi,\Psi)=\min_{0<y<1}\Big[\frac{\Phi^2}{y}+\frac{\Psi^2}{1-y}-a^2(1-y)-b^2y\Big].
$$

In code, test whether $\Theta\ge0$ somewhere using the cubic (§B3, step 3).

### Shape of the allowed region $\mathcal A(a,b)$

Scripts: `region.py`, `region_maps.png`.

**General shape.**

- $\mathcal A$ is a **compact, convex, curvilinear quadrilateral** ("curved diamond") containing the origin.
- Tests on 8 black holes: 0 star-shapedness violations along 60 rays × 80 points; 0 convexity violations in 300 random chords.
- The origin is the zero-momentum orbit, which sweeps all latitudes (Figs 4.2 and 4.6).

**Vertices on the axes.** These are the paper's circular orbits (3.9)–(3.15), confirmed by ray bisection to $10^{-8}$:

| Vertex | Orbit type | $\Phi$ or $\Psi$ | $x$ |
|---|---|---|---|
| $\Psi=0$, $\Phi_{\max}$ | equatorial counter-rotating | $1+\sqrt{(1+a)^2-b^2}$ | $1+a-b^2+\sqrt{(1+a)^2-b^2}$ |
| $\Psi=0$, $\Phi_{\min}$ | equatorial co-rotating | $-1-\sqrt{(1-a)^2-b^2}$ | $1-a-b^2+\sqrt{(1-a)^2-b^2}$ |
| $\Phi=0$, $\Psi_{\max}$ | polar counter-rotating | $1+\sqrt{(1+b)^2-a^2}$ | $1+b-a^2+\sqrt{(1+b)^2-a^2}$ |
| $\Phi=0$, $\Psi_{\min}$ | polar co-rotating | $-1-\sqrt{(1-b)^2-a^2}$ | $1-b-a^2+\sqrt{(1-b)^2-a^2}$ |

These are written for $a,b\ge0$. The vertices are axis intercepts, not necessarily the extreme $\Phi$ or $\Psi$ of the region. For example, for $a=b=0.49$ the minimum $\Phi$ is $-1.209$, against $\Phi_{\min}=-1.141$.

**Edges.** The edges are the **circular orbits of constant $\theta_0$**: the $\theta$-band width on the boundary is 0 to within the bisection tolerance. Along each quadrant's edge, $\theta_0$ runs from $\pi/2$ (at the $\Phi$-axis vertex) to $0$ (at the $\Psi$-axis vertex). Fig 4.12 is one such edge point. For $a=0.5$, $b=0.25$, $\theta_0=\pi/4$, the four edge points are:

| $\Phi$ | $\Psi$ | $x$ | $K$ |
|---|---|---|---|
| −0.764662 | −0.733371 | 0.961862 | 2.088834 |
| −0.929241 | +0.903667 | 1.455284 | 3.203956 |
| +1.048754 | −1.026162 | 1.906766 | 4.149537 |
| +1.168941 | +1.148716 | 2.519041 | 5.215697 |

The second row is the paper's Fig 4.12.

**Interior.**

- $\theta$ oscillates in $[\theta_-,\theta_+]$.
- On the line $\Phi=0$ the orbit passes through the $\phi$-axis ($\theta=0$).
- On the line $\Psi=0$ it passes through $\theta=\pi/2$.

**Radius across the region.** $x$ varies smoothly over $\mathcal A$. Its minimum and maximum are attained at the co-rotating and counter-rotating vertices of the plane with the larger spin; for $a=b$ they are attained along the straight edges. I checked this with 720 boundary points per black hole (`xextreme.py`). So $x_{\min}$ and $x_{\max}$ are the vertex radii in the table above. Examples: for $(0.6,0.3)$, $x\in[0.574575,3.081623]$; for $(0.9,0.09)$, $x\in[0.135489,3.789767]$.

**Special cases.**

- **Tangherlini:** $\mathcal A=\{|\Phi|+|\Psi|\le2\}$, with $x\equiv2$ and $K\equiv4$.
- **$a=b$:** $x$ and $K$ depend only on $\Phi+\Psi$. The edges in the $(+,+)$ and $(-,-)$ quadrants are **straight lines** $\Phi+\Psi=$ const.
- **$b=0$:** $\mathcal A$ is symmetric under $\Psi\to-\Psi$.

**Outside $\mathcal A$.** The radial solution still exists, but $\Theta<0$ for all $\theta$, so there is **no orbit**.

### Computing the boundary

**Option 1: ray bisection (robust; recommended).**

1. For a direction $\alpha$, find $s^*(\alpha)=\sup\{s:\ (s\cos\alpha,s\sin\alpha)\in\mathcal A\}$ by bisection. This works because $\mathcal A$ is star-shaped about the origin.
2. Each test is one radial solve plus one band test, about 40 bisection steps per ray.
3. 64–128 rays draw the region on a 2D pad in a few milliseconds.

**Option 2: exact parametrisation by $\theta_0$.** Write $y_0=\sin^2\theta_0$. The conditions $\Theta(y_0)=\Theta'(y_0)=0$ give

$$
\Phi^2=y_0^2\big[K+2a^2(1-y_0)+b^2(2y_0-1)\big],\qquad \Psi^2=(1-y_0)^2\big[K+2b^2y_0+a^2(1-2y_0)\big]
$$

(verified). Insert these into $\mathcal X=\mathcal X'=0$ with the quadrant's signs, and solve the two equations for $(x,K)$ by Newton's method, continuing in $y_0$ from the $\Phi$-axis vertex ($y_0=1$).

### Recommended solver

Inputs are $(a,b,\Phi,\Psi)$; the outputs are $x$ (or $r$), $K$, and the band $[y_1,y_2]$. The reference implementation is `reference_impl.py`, written in plain floats so it ports directly to JavaScript.

1. Compute $x_+$. Set $lo=x_+(1+10^{-12})+10^{-14}$.
   - $F(lo)<0$ always holds off extremality and away from $\mathcal E(x_+)=0$.
   - If $F(lo)\ge0$, fall back to the polynomial roots of $F$.
2. Set $hi=\max(4,2\,lo)$ and double it until $F(hi)>0$.
3. Bisect on $F$ for 60–200 iterations, or use Brent. This gives $x_0$.
4. Set $K=P(x_0)/\Delta(x_0)$.
5. Run the band test (§B3, step 3). If there is no band, the point $(\Phi,\Psi)$ is outside $\mathcal A$.
6. **No solution:** project back to the boundary along the ray, $(\Phi,\Psi)\to(1-\epsilon)\,s^*(\alpha)(\cos\alpha,\sin\alpha)$, or refuse the move and show "no spherical orbit: θ-motion forbidden".
7. **Multiple solutions:** never observed. Code defensively anyway: if the quartic shows more than one root with $x>x_+$, keep those that are minima of $V$ and pass the band test, and log a warning.

### Which controls to show the user

**Recommendation: $(\Phi/E,\Psi/E)$, presented on a 2D pad.** The pad draws $\mathcal A(a,b)$ (for example via 96 rays) and shows the current point as a draggable dot; optionally it is coloured by $r$, as in `region_maps.png`. Two linked sliders with live clamping are the fallback.

Reasons:

- every point of $\mathcal A$ is exactly one orbit, and every orbit is a point;
- the region is convex and contains the origin, so clamping and interpolating stay inside it;
- the axes and edges have clean meanings: equatorial, polar and constant-$\theta$ circular orbits;
- it works uniformly, including Tangherlini.

**Persist the orbit "type" when $a$ or $b$ moves.** Store $(\alpha,f)$ with $\alpha=\operatorname{atan2}(\Psi,\Phi)$ and $f=s/s^*(\alpha)\in[0,1]$. Recompute $(\Phi,\Psi)=f\,s^*(\alpha)(\cos\alpha,\sin\alpha)$ when $a$ or $b$ changes. In this scheme:

- $f=0$ is the zero-momentum orbit;
- $f=1$ is a circular constant-$\theta$ orbit;
- $\alpha$ sets the momentum ratio $\tan\alpha=\Psi/\Phi$.

These two numbers would also make good slider labels ("direction" and "fraction of maximum").

**Alternatives considered.**

*Radius plus a mixing parameter.*

- At fixed $x$, $F=0$ is a **conic** in $(\Phi,\Psi)$, because $F$ is quadratic in the momenta.
- The orbits at that radius form the arc of the conic inside $\mathcal A$, ending at two constant-$\theta$ circular orbits.
- This works, but:
  - the radius range $[x_{\min},x_{\max}]$ shrinks to a point as $a,b\to0$ (Tangherlini: $x\equiv2$);
  - the arc parameter has no intuitive meaning;
  - it needs a conic-arc solver.
- It is worse on every count except "radius is directly controllable". You can still *display* $r$ as a readout, or colour the pad by $r$.

*Ratio $\Psi/\Phi$ plus magnitude.* This is essentially the $(\alpha,f)$ scheme above, which is the best normalised choice.

## B3. Equations of motion for a fast and robust JavaScript integration

### Recommended formulation: Mino time plus an angle variable

**Do not integrate $x$.** It is exactly constant on the orbit, and integrating it makes it drift away within about 6 oscillations (§A3).

With Mino time $d\sigma=d\lambda/\rho^2$ (as in KR12) and $y=\sin^2\theta$, everything decouples and depends only on $y(\sigma)$:

$$
\Big(\frac{dy}{d\sigma}\Big)^2=4Q(y),\qquad Q(y)=(a^2-b^2)y^3+(b^2-2a^2-K)y^2+(\Phi^2-\Psi^2+a^2+K)y-\Phi^2,
$$

$$
\frac{d\phi}{d\sigma}=\frac{\Phi}{y}-A_\phi,\quad \frac{d\psi}{d\sigma}=\frac{\Psi}{1-y}-A_\psi,\quad \frac{dt}{d\sigma}=\rho^2+A_t,\quad \frac{d\lambda}{d\sigma}=\rho^2=x+a^2+(b^2-a^2)\,y,
$$

with constants on the orbit

$$
\mathcal E=1+\frac{a\Phi}{x+a^2}+\frac{b\Psi}{x+b^2},\quad
A_\phi=\frac{\mu a(x+b^2)\mathcal E}{\Delta}+\frac{(a^2-b^2)\Phi}{x+a^2},\quad
A_\psi=\frac{\mu b(x+a^2)\mathcal E}{\Delta}-\frac{(a^2-b^2)\Psi}{x+b^2},\quad
A_t=\frac{\mu(x+a^2)(x+b^2)\mathcal E}{\Delta}.
$$

**The angle variable removes the turning points.** Let $[y_1,y_2]$ be the band, $\bar y=(y_1+y_2)/2$ and $\Delta y=(y_2-y_1)/2$. Then

$$
Q(y)=(y-y_1)(y_2-y)\,G(y),\qquad G(y)=(b^2-a^2)(y+y_1+y_2)-(b^2-2a^2-K)>0 \text{ on the band}
$$

(verified, valid also for $a=b$). Set

$$
y(\chi)=\bar y-\Delta y\cos\chi,\qquad \frac{d\sigma}{d\chi}=\frac{1}{2\sqrt{G(y(\chi))}} .
$$

Then $\chi$ increases monotonically with no turning points and no sign flips. $\chi=0$ is $\theta_{\min}$ and $\chi=\pi$ is $\theta_{\max}$. **Every right-hand side depends on $\chi$ only**, so the "integration" is a quadrature:

$$
\frac{d(\sigma,\lambda,t,\phi,\psi)}{d\chi}=\frac{1}{2\sqrt G}\Big(1,\ \rho^2,\ \rho^2+A_t,\ \frac{\Phi}{y}-A_\phi,\ \frac{\Psi}{1-y}-A_\psi\Big).
$$

### The algorithm

1. Solve for the orbit as in §B2.

2. Compute the constants $A_\phi$, $A_\psi$, $A_t$ and $G$.

3. **Band test.**
   - Form $R(y)=Q(y)$. If $\Phi=0$, divide by $y$; if $\Psi=0$, divide by $(1-y)$. This matters because then $\Theta(0)$ or $\Theta(1)$ is finite.
   - Candidates for the maximum are the critical points of $R$ in $(0,1)$ (roots of the quadratic $R'$, solved exactly), plus $y=0$ only if $\Phi=0$, plus $y=1$ only if $\Psi=0$.
   - If $\max R<0$: no orbit.
   - If $\max R\approx0$: a circular orbit with $y_1=y_2$.
   - Otherwise, bisect for $y_1$ and $y_2$.
   - *Pitfall I hit:* testing $Q\ge0$ at the endpoints with a tolerance gives false positives for tiny nonzero $\Phi$, because $Q(0)=-\Phi^2\approx0$.

4. **Tabulate one $\chi$-period** with fixed step $h=2\pi/N$ using Simpson's rule, which is the same as RK4 since the right-hand side depends on $\chi$ only.
   - **Subtract the axis spikes analytically.** Write
     $\frac{\Phi}{2y\sqrt G}=\frac{\Phi}{2y\sqrt{G(y_1)}}+\frac{\Phi}{2y}\big(\tfrac{1}{\sqrt{G(y)}}-\tfrac{1}{\sqrt{G(y_1)}}\big)$.
     The second term is smooth.
   - The first term integrates in closed form:
     $$\int_0^\chi\frac{d\chi'}{\bar y-\Delta y\cos\chi'}=\frac{2}{\sqrt{y_1y_2}}\Big[\tfrac{\chi}{2}+\operatorname{atan2}\big((k-1)\sin\tfrac\chi2\cos\tfrac\chi2,\ \cos^2\tfrac\chi2+k\sin^2\tfrac\chi2\big)\Big],\qquad k=\sqrt{y_2/y_1}.$$
     This is continuous in $\chi$.
   - Do the same for $\Psi/(1-y)$, subtracting with $G(y_2)$, using $(1-y_2,1-y_1)$ in place of $(y_1,y_2)$, and evaluating the closed form at $\chi+\pi$ minus its value at $\pi$. This works because $1-y=(1-\bar y)-\Delta y\cos(\chi+\pi)$.

5. **Periodic extension.** $(\sigma,\lambda,t,\phi,\psi)(\chi+2\pi n)=\text{table}(\chi)+n\cdot(\Delta\sigma,\Delta\lambda,\Delta t,\Delta\phi,\Delta\psi)$. Any history length then costs nothing, and the error does not accumulate.

**Accuracy** (`reference_impl.out`). The error is against DOP853 at rtol $10^{-12}$; "naive" means without the subtraction.

| Case | $N$ | error in $\Delta\phi$ (subtracted) | error in $\Delta\phi$ (naive) | error in $\Delta\psi$ (subtracted) | error in $\Delta\psi$ (naive) |
|---|---|---|---|---|---|
| $\Phi=10^{-3}$, $\Psi=0.5$, $a=0.6$, $b=0.3$ | 128 | $1\times10^{-7}$ | **13 rad** | $10^{-11}$ | $10^{-11}$ |
| same | 512 | $1\times10^{-8}$ | 1.4 rad | | |
| Fig 4.10 ($\Phi=0.1$, $\Psi=0.01$) | 128 | $2\times10^{-9}$ | $7\times10^{-6}$ | $5\times10^{-7}$ | **0.14 rad** |
| Generic cases | 32 | $\le10^{-11}$ | | | |

**Recommendation:** $N=256$ with subtraction, and $N=512$ if $\min(y_1,1-y_2)<10^{-4}$. One period costs about 1,500 square roots, i.e. microseconds.

### Special handling

**$\Phi=0$ (or $\Psi=0$).**

- The axis term is absent.
- The orbit passes *through* the axis, so use **signed** coordinates for drawing:
  $\sin\theta_s=\sqrt{y_2}\sin(\chi/2)$ when $\Phi=0$, and $\cos\theta_s=\sqrt{1-y_1}\cos(\chi/2)$ when $\Psi=0$.
- This is exactly what the paper's $[0,\pi]$ convention does in effect. In terms of the embedding, it is a phase jump of $\pi$ in $\phi$ (or $\psi$) per $\chi$-period.

**Circular orbits** ($y_1=y_2$, on the edge of $\mathcal A$): the same code works, with $\Delta y=0$.

**Initial data.** To start at $\theta_0$ moving towards larger $\theta$, set $\chi_0=\arccos\big((\bar y-\sin^2\theta_0)/\Delta y\big)$. For decreasing $\theta$, use $2\pi-\chi_0$.

### Closed forms (for the record and for frequencies)

**$\theta$ motion.** Let $G_{\max}$ and $G_{\min}$ be the larger and smaller of $G(y_1),G(y_2)$, and set $m=1-G_{\min}/G_{\max}\in[0,1)$. Then

$$
y(\sigma)=y_e\pm(y_2-y_1)\,\mathrm{sn}^2\!\big(\sqrt{G_{\max}}\,\sigma\,\big|\,m\big),
$$

where $y_e$ is the endpoint at which $G=G_{\max}$, and $\sigma=0$ there. The $\sigma$-period is

$$
\Lambda_\sigma=\frac{2K(m)}{\sqrt{G_{\max}}}.
$$

Both are verified to $10^{-14}$ (`check_jacobi.py`). $m$ is small in practice (at most 0.2 in all the tests), so $y$ is nearly sinusoidal in Mino time.

**$\phi,\psi$** are incomplete elliptic integrals of the third kind, and **$t,\lambda$** of the second kind.

**For $a=b$, everything is elementary.** $G=a^2+K$ is constant, $\chi=2\sqrt{a^2+K}\,\sigma$, and per $\chi$-period

$$
\Delta\phi=\pi\Big[\operatorname{sgn}\Phi-\frac{A_\phi}{\sqrt{a^2+K}}\Big],\qquad \Delta\psi=\pi\Big[\operatorname{sgn}\Psi-\frac{A_\psi}{\sqrt{a^2+K}}\Big],\qquad A_\phi=A_\psi,
$$

verified to $10^{-10}$.

**Is the closed form worth using?** Not for drawing: the quadrature above is simpler and just as accurate. The complete integrals are worth having for frequency readouts and resonance search (§B6). $K(m)$ is one AGM: $K(m)=\pi/[2\,\mathrm{AGM}(1,\sqrt{1-m})]$.

### Which "time" for the history slider and the animation

| Variable | Pros | Cons |
|---|---|---|
| Coordinate time $t$ | The distant observer's clock, so physically meaningful. Strictly monotonic ($dt/d\chi>0$ on every valid orbit). | |
| Affine parameter $\lambda$ | What the paper plots. | Its normalisation depends on $E$. |
| Mino time $\sigma$ | Makes the $\theta$ motion periodic with a fixed period. | No direct physical meaning. |

All three are smooth and monotonic in $\chi$. In practice $\Delta t$ per $\theta$-oscillation lies between 5.9 and 6.4 for every case tested (exactly $2\pi$ for Tangherlini), and $t\approx2\lambda$.

**Recommendation.**

- **History slider:** coordinate time $t$, displayed alongside the equivalent number of $\theta$-oscillations ($\approx t/\Delta t$).
- **Animation:** advance $t$ at a constant rate.
- Tabulate in $\chi$ and invert $t(\chi)$ by interpolation. It is monotonic, so binary search works.

### Integrator: the brief's choice (item 4)

The brief proposes adaptive RK45 in a Web Worker instead of fixed-step RK4, because Leo Stein's Kerr page notes fixed-step inaccuracy near the axis. My assessment:

- **That inaccuracy belongs to the Boyer-Lindquist $\theta$ equation**: turning points, and the $\Phi/\sin^2\theta$ spike near the axis.
- **The $\chi$ formulation with spike subtraction removes both.** Fixed-step Simpson with $N=256$ is accurate to $10^{-7}$ even at $\Phi=10^{-3}$.
- One period is a sub-millisecond computation, and periodic extension makes history free, so **no Web Worker is needed**.
- Adaptive RK45 on the $\chi$ system *without* subtraction would also work. It is acceptable, just unnecessary.
- **What must not be done** is adaptive or fixed-step integration of the full Boyer-Lindquist system including $x$ and $p_x$. That system leaves the orbit after about 6 oscillations, whatever the tolerance (§A3).

### Side panel (item 4)

Plotting $\theta,\phi,\psi$ against $\lambda$ reproduces Figs 4.9 and 4.11 one-for-one, which is handy for validating presets. My suggestions:

- default to the **same variable as the history slider** ($t$), with a toggle to $\lambda$;
- plot $\phi,\psi$ **unwrapped**, as the paper does (the jumps of about $\pi$ near the axes are the point of Fig 4.11);
- plot $\theta\in[0,\pi/2]$, or signed for orbits that cross an axis;
- draw dashed lines at $\theta_\pm$;
- drop the paper's "x" line, which is really $r$ and constant. Show $r$ as a number instead.

## B4. Embedding and projection

### The paper's projection, exactly as used for its figures

$$
x_1=\sqrt{r^2+a^2}\sin\theta\cos\!\big[\phi-\arctan\tfrac ar\big],\quad
y_1=\sqrt{r^2+a^2}\sin\theta\sin\!\big[\phi-\arctan\tfrac ar\big],\quad
z_1=\sqrt{r^2+b^2}\cos\theta,
$$

with $\theta$ in the paper's double-covered range, so $z_1$ can be negative.

**Colour.** The colour is $\mathrm{Hue}[\psi/2\pi]$ (HSV hue $=\psi/2\pi \bmod 1$, saturation 1, value 1). I verified this by sampling the Fig 4.1 raster:

| $\psi$ | Colour |
|---|---|
| 0 | red |
| $\pi/3$ | yellow |
| $2\pi/3$ | green |
| $\pi$ | cyan |
| $4\pi/3$ | blue |
| $5\pi/3$ | magenta |

The hue increases counter-clockwise on the wheel. The colours of the reproduced Fig 4.10 arcs ($\psi\approx2.9$ cyan, $\psi\approx5.85$ crimson) match the paper.

**Check.** The $\arctan(a/r)$ offset is a constant rotation for constant $r$; drop it. The underlying map is the exact flat-space ($\mu\to0$) spheroidal embedding of $\mathbb R^4$: I verified symbolically that $dX_1^2+\dots+dY_2^2$ pulls back to the spatial part of the metric at $\mu=0$. In that embedding:

- $X_1+iY_1=\sqrt{x+a^2}\sin\theta\,e^{i\phi}$, so **$a$ pairs with $\phi$ and $\sin\theta$**;
- $X_2+iY_2=\sqrt{x+b^2}\cos\theta\,e^{i\psi}$, so **$b$ pairs with $\psi$ and $\cos\theta$**;
- the constant-$r$ surface is the ellipsoid $\frac{X_1^2+Y_1^2}{x+a^2}+\frac{X_2^2+Y_2^2}{x+b^2}=1$;
- its image in $(x_1,y_1,z_1\ge0)$ is half of a spheroid with equatorial semi-axis $\sqrt{x+a^2}$ and polar semi-axis $\sqrt{x+b^2}$. **It is oblate iff $a^2>b^2$, prolate iff $a^2<b^2$, and a sphere iff $a^2=b^2$.** ✓

### Comparison with "three of $(X_1,Y_1,X_2,Y_2)$ drawn plus the fourth as colour"

The paper's map is the $\mathbb R^4$ embedding written in **polar form in the second plane**: height $z_1=|X_2+iY_2|$ and hue $\psi=\arg(X_2+iY_2)$. Height and hue together recover $(X_2,Y_2)$ exactly, so the picture is **lossless** apart from $z_1=0$, where $\psi$ is undefined anyway. The orbit lies on a 2D surface that can be drawn opaque, as the paper does. And $\psi$ is cyclic, so a cyclic colour map is the natural choice.

The alternative, drawing $(X_1,Y_1,X_2)$ and colouring by $Y_2$, is an orthogonal projection:

- the curve then wanders *inside* a solid ellipsoid, with no surface to anchor it, which is much harder to read;
- the colour mixes the amplitude $\cos\theta$ with the angle $\psi$;
- it needs a diverging, non-cyclic map.

**Recommendation: use the paper's projection.**

- Signed $z_1=\sqrt{x+b^2}\cos\theta_s$ (§B3), so that $\Psi=0$ orbits cross the white "equator" smoothly, as in Figs 4.2, 4.3 and 4.6. Remember that the lower half then shows the same points as the upper half with $\psi\to\psi+\pi$.
- Hue $=\psi/2\pi$ to match the paper. Optionally offer a perceptually uniform cyclic map (CET-C2, cmocean "phase" or "twilight"); HSV has bright and dark bands, at yellow and blue.
- Draw the orbit's spheroid translucent or dark, and the horizon spheroid ($x=x_+$) inside it.

**Optional "swap planes" toggle.** Draw $(X_2,Y_2,\sqrt{x+a^2}\sin\theta_s)$ and colour by $\phi$. This makes the $a\leftrightarrow b$ symmetry visible: the paper notes that Figs 4.2 and 4.6 are the same orbit type, which looks different only because the projection privileges the $a$-plane.

**Stereographic projection of $S^3$** (after normalising each plane) shows all four dimensions as position, with constant-$\theta$ surfaces becoming nested Clifford tori. It is attractive as an "advanced view", but it is unbounded near the projection point, unfamiliar, and not comparable with the paper or the Kerr tool. I would not make it the default.

### Coordinate singularities

$\theta=0$ is the $\phi$-axis. In the picture it is the top pole $X_1=Y_1=0$, where $\phi$ is undefined:

- only $\Phi=0$ orbits reach it, and they pass straight through when the signed $\sin\theta_s$ is used;
- small $\Phi\neq0$ produces a fast swing of $\phi$ by about $\pi$ near the pole. In the drawing this is a tiny loop, handled exactly by the subtraction.

$\theta=\pi/2$ is the $\psi$-axis: the white ring $z_1=0$, where the hue is undefined.

- Only $\Psi=0$ orbits reach it.
- Small $\Psi\neq0$ gives a rapid change of hue by about $\pi$. That is the "cusp" and "sharp colour change" effect of Figs 4.8 and 4.10, and it is real.

**Neither singularity affects the $\chi$ equations.**

## B5. Presets reproducing Figs 4.2–4.12

Scripts: `figures.py`, `fit_figs.py`, `digitize.py`; outputs `figures.out`, `fit_figs.out`, `figures_repro.png`, `figs_fitted_3d.png`.

**Conventions for the table.**

- All presets have $\mu=1$ and $E=1$.
- $\theta$ ranges are given in the 5D convention, $[0^\circ,90^\circ]$. The paper's figures show $\Psi=0$ orbits mirrored through $z_1=0$.
- "per osc." means per $\chi$-period, i.e. one $y$-oscillation $\theta_{\min}\to\theta_{\max}\to\theta_{\min}$. For orbits through an axis, the paper's "one latitudinal oscillation" (pole → pole → pole) is **two** $\chi$-periods.
- The paper states no initial conditions or $\lambda$ ranges except implicitly in Figs 4.9 and 4.11. Those were fitted from the digitised curves, with rms 0.004–0.006 rad. The other ranges are my estimates from the figures and are marked "est."

| Fig | $a$ | $b$ | $\Phi$ | $\Psi$ | $x$ (=$r^2$) | $r$ | $K$ | $\theta$ range | per osc.: $\Delta\lambda$, $\Delta t$, $\Delta\phi$, $\Delta\psi$ | Shown, and suggested start |
|---|---|---|---|---|---|---|---|---|---|---|
| 4.2 | 0.5 | 0 | 0 | 0 | 1.750000 | 1.322876 | 3.750000 | all (0–90°, through both axes) | 2.9913, 6.1838, −0.7981, 0 | "one latitudinal oscillation" = 2 osc. (λ≈0–6); start $\theta=0$, φ=ψ=0 |
| 4.3 | 0.5 | 0 | 2 | 0 | 2.750000 | 1.658312 | 5.750000 | 55.94°–90° (paper: 55.94°–124.06°) | 3.6970, 6.3480, 2.5362, 0 | "many": est. 20+ osc. |
| 4.4 | 0.8 | 0 | 0 | 0.5 | 1.316830 | 1.147532 | 3.240444 | 0–74.00° | 2.7313, 6.1061, −1.3797 (+π at each pole passage), 3.4643 | est. 2 osc. |
| 4.5 | 0.8 | 0 | 0 | 0.5 | same | | | same | same | est. ~20 loops (≈20 osc.) |
| 4.6 | 0 | 0.99 | 0 | 0 | 1.019900 | 1.009901 | 3.019900 | all | 2.5124, 5.8785, 0, −1.6662 | 2 osc.; $\phi$ constant |
| 4.7 | 0.6 | 0.3 | 1 | −1 | 1.754557 | 1.324597 | 4.020430 | 37.59°–51.53° | 3.0191, 6.1240, 2.1133, −3.8200 | "many": est. ~20 osc. The orbit is close to the $\nu_\phi=1/3$ resonance ($\nu_\phi=0.33635$), which is why the figure shows about 7 slowly precessing strands, each of 3 oscillations. |
| 4.8/4.9 | 0.9 | 0.1 (extremal) | 0.8 | 0.8 | 1.899039 | 1.378056 | 4.522172 | 22.23°–66.10° | 3.2663, 6.1839, 1.9249, 3.3649 | Fig 4.9: λ∈[0,5]. Fitted at λ=0: θ=0.3924 (increasing, just past θ_min=0.3880), φ=1.748, ψ=0.449. Fig 4.8 shows ≈6 arcs, i.e. ≈6 osc. (λ≈0–20, est.) |
| 4.10/4.11 | 0.9 | 0.09 | 0.1 | 0.01 | 1.267423 | 1.125799 | 3.360570 | 2.81°–89.69° | 2.7007, 5.9469, 1.6886, 2.9209 | λ∈[1,5] shown. Fitted at λ=1: θ=1.021 (decreasing), φ=−0.980, ψ=2.957 (at λ=0: 0.781, −0.416, −0.072) |
| 4.12 | 0.5 | 0.25 | −0.929241 | +0.903667 | 1.455284 | 1.206351 | 3.203956 | 45° (circular) | dφ/dt=−0.64847, dψ/dt=+0.43978 | one turn; paper's rounded values give θ∈[44.69°, 45.32°] |

**Feature checks, all reproduced** (`figures_repro.png`).

| Fig | Feature in the paper | What my integration gives |
|---|---|---|
| 4.2 | All latitudes, non-planar | $\phi$ drifts by −0.80 per osc. through frame dragging; ψ constant (red) |
| 4.3 | Band about the equator | Band edges at $z/R=\cos55.9^\circ=0.56$, matching the figure; quasiperiodic, with $\Delta\phi/2\pi=0.4037$ |
| 4.4/4.5 | Returns to the pole, has a minimum latitude, ψ increasing | All reproduced: flower pattern, θ_max = 74° |
| 4.6 | Prolate "football", $\phi$ constant | Reproduced, with axis ratio $\sqrt{2}/1.01=1.40$ |
| 4.7 | Band in (0, π/2) | Narrow band 37.6°–51.5° with a wave pattern |
| 4.8 | Cusps | At $\theta_{\max}$, $d\phi/d\sigma=0.034$ while $d\psi/d\sigma=5.1$, giving near-cusps with sharp hue changes. The 3D rendering shows the sharp corners. |
| 4.10 | φ decreases near the equator and increases near the pole; ψ the opposite | $d\phi/d\sigma=+40.8$ at $\theta_{\min}$ and $-0.81$ at $\theta_{\max}$; $d\psi/d\sigma=-0.125$ at $\theta_{\min}$ and $+335.7$ at $\theta_{\max}$. Exactly the paper's description. Arc colours (cyan, crimson) match. |
| 4.12 | Constant θ | θ ≡ π/4 exactly at the exact values |

**Unresolved figure detail: scale.** In Mathematica's perspective 3D plots I could not reliably read the absolute scale of the drawn spheroid from the tick labels. The spheroid in some figures looks closer to unit radius than to $\sqrt{x+a^2}$, so the author may have normalised. This does not affect the shapes or $\theta$ ranges.

## B6. Frequencies and resonances (stretch goal)

**The paper only says "quasiperiodic".** It has no frequency analysis.

### Definitions

Per $\chi$-period (one $y$-oscillation), define

$$
\nu_\phi=\frac{\Delta\phi_{\rm eff}}{2\pi},\qquad \nu_\psi=\frac{\Delta\psi_{\rm eff}}{2\pi},
$$

with $\Delta\phi_{\rm eff}=\Delta\phi+\pi$ if $\Phi=0$, and likewise for $\psi$. The extra $\pi$ is the phase jump at an axis passage. With it, $\nu$ is continuous mod 1 across $\Phi=0$ (verified: at $\Phi=-10^{-4}$, $0$, $+10^{-4}$ the values are $-0.6435$, $0.3565$, $0.3565$).

The coordinate-time frequencies are $\Omega_\theta=2\pi/\Delta t$, $\Omega_\phi=\Delta\phi_{\rm eff}/\Delta t$ and $\Omega_\psi=\Delta\psi_{\rm eff}/\Delta t$.

### Closure condition

The spatial orbit closes in $\mathbb R^4$ after $n$ oscillations iff $n\nu_\phi\in\mathbb Z$ and $n\nu_\psi\in\mathbb Z$. **Both ratios must be rational.** The time $t$ does not need to close.

| Case | Closure structure |
|---|---|
| Tangherlini | $\nu_\phi=\nu_\psi=\tfrac12$ for **every** orbit (great circles; verified to $10^{-10}$). Everything closes after 2 oscillations. Frequency-map Jacobian: rank 0. |
| $a=b$ | $\nu_\phi-\tfrac12\operatorname{sgn}\Phi=\nu_\psi-\tfrac12\operatorname{sgn}\Psi=-A/(2\sqrt{a^2+K})$, a function of $\Phi+\Psi$ only. Take $\operatorname{sgn}0=+1$, consistent with the $+\pi$ rule. Closure is a single condition, so closed orbits lie on lines $\Phi+\Psi=$ const. Jacobian: rank 1. |
| Generic $a\neq b$ | The Jacobian of $(\Phi,\Psi)\mapsto(\nu_\phi,\nu_\psi)$ has rank 2 (e.g. det $=6.3\times10^{-4}$ at Fig 4.7's parameters, singular values 0.128 and 0.0049), so closed orbits are **isolated points** of $\mathcal A$. The second singular value is small because the map becomes degenerate as $a\to b$, so undamped Newton steps can be large; use damping. |

### Resonance hunt recipe

1. Compute $\nu(\Phi,\Psi)$ from the one-period table ($N\approx1024$).
2. Pick a target $(p/n,\ q/n)$ near the current $\nu$.
3. Use damped Newton with a finite-difference Jacobian, backtracking whenever a step leaves $\mathcal A$.

**Example** (`resonance.py`). For $a=0.6$, $b=0.3$ the orbit with $\nu=(1/3,-3/5)$ is $\Phi=0.1051348934$, $\Psi=-0.1354415960$, with $x=1.53408704$ and $K=3.56062587$. It closes after 15 oscillations: $15\,\Delta\phi/2\pi=5.000000000$ and $15\,\Delta\psi/2\pi=-9.000000000$.

## B7. Unit-test values

These come from `tables.py`, with the reference implementation at $N=2048$. They agree with the independent DOP853 and Hamiltonian checks to at least $10^{-8}$.

- $\mu=E=1$; $\theta$ in degrees.
- Increments are per $\chi$-period (one $\sin^2\theta$ oscillation), and $\nu=\Delta_{\rm eff}/2\pi$.
- The input $\Phi=-0.929241$, $\Psi=0.903667$ is Fig 4.12 rounded to 6 digits, hence the $\pm0.013^\circ$ band.

| a | b | Φ | Ψ | x₊ | x = r² | r | K | θ_min | θ_max | Δσ | Δλ | Δt | Δφ | Δψ | ν_φ | ν_ψ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 1.000000 | 2.000000 | 1.414214 | 4.000000 | 0.0000 | 90.0000 | 1.570796 | 3.141593 | 6.283185 | 0 | 0 | 0.5 | 0.5 |
| 0 | 0 | 0.7 | 0.9 | 1.000000 | 2.000000 | 1.414214 | 4.000000 | 23.6955 | 60.5654 | 1.570796 | 3.141593 | 6.283185 | 3.141593 | 3.141593 | 0.5 | 0.5 |
| 0.5 | 0 | 0 | 0 | 0.750000 | 1.750000 | 1.322876 | 3.750000 | 0.0000 | 90.0000 | 1.596242 | 2.991345 | 6.183829 | −0.798121 | 0 | 0.372975 | 0.5 |
| 0.5 | 0 | 2 | 0 | 0.750000 | 2.750000 | 1.658312 | 5.750000 | 55.9375 | 90.0000 | 1.325491 | 3.696987 | 6.347970 | 2.536226 | 0 | 0.403653 | 0.5 |
| 0.8 | 0 | 0 | 0.5 | 0.360000 | 1.316830 | 1.147532 | 3.240444 | 0.0000 | 73.9964 | 1.650188 | 2.731270 | 6.106098 | −1.379712 | 3.464257 | 0.280412 | 0.551354 |
| 0 | 0.99 | 0 | 0 | 0.019900 | 1.019900 | 1.009901 | 3.019900 | 0.0000 | 90.0000 | 1.683065 | 2.512396 | 5.878526 | 0 | −1.666234 | 0.5 | 0.234811 |
| 0.4 | 0.4 | 1.2 | −0.5 | 0.640000 | 1.967354 | 1.402624 | 4.232685 | 36.8582 | 72.6529 | 1.498942 | 3.188781 | 6.189437 | 2.577388 | −3.705797 | 0.410204 | −0.589796 |
| 0.6 | 0.3 | 1 | −1 | 0.482906 | 1.754557 | 1.324597 | 4.020430 | 37.5908 | 51.5256 | 1.523646 | 3.019127 | 6.123974 | 2.113327 | −3.820010 | 0.336346 | −0.607973 |
| 0.9 | 0.1 | 0.8 | 0.8 | 0.090000 | 1.899039 | 1.378056 | 4.522172 | 22.2319 | 66.0957 | 1.411811 | 3.266340 | 6.183856 | 1.924896 | 3.364853 | 0.306357 | 0.535533 |
| 0.9 | 0.09 | 0.1 | 0.01 | 0.132313 | 1.267423 | 1.125799 | 3.360570 | 2.8074 | 89.6874 | 1.622139 | 2.700653 | 5.946941 | 1.688610 | 2.920855 | 0.268751 | 0.464868 |
| 0.55 | 0.2 | −0.6 | 0.4 | 0.638551 | 1.406699 | 1.186043 | 3.157429 | 19.4371 | 76.2465 | 1.727826 | 2.712565 | 6.168261 | −4.134277 | 2.761033 | −0.657991 | 0.439432 |
| 0.15 | 0.7 | 0.9 | −0.6 | 0.463725 | 1.196022 | 1.093628 | 2.762068 | 35.5266 | 66.1852 | 1.773146 | 2.645541 | 6.302181 | 3.176527 | −4.440477 | 0.505560 | −0.706724 |
| 0.5 | 0.25 | −0.929241 | 0.903667 | 0.663967 | 1.455284 | 1.206351 | 3.203956 | 44.9869 | 45.0131 | 1.713827 | 2.761890 | 6.193134 | −4.016074 | 2.723607 | −0.639178 | 0.433475 |
| 0.6 | 0.3 | 0 | 0 | 0.482906 | 1.521687 | 1.233567 | 3.528757 | 0.0000 | 90.0000 | 1.621893 | 2.830970 | 6.086223 | −1.037979 | −0.605934 | 0.334800 | 0.403563 |
| 0.6 | −0.3 | 1 | 1 | 0.482906 | 1.754557 | 1.324597 | 4.020430 | 37.5908 | 51.5256 | 1.523646 | 3.019127 | 6.123974 | 2.113327 | 3.820010 | 0.336346 | 0.607973 |

The last row is a **symmetry test**: $b\to-b$ with $\Psi\to-\Psi$ must reproduce the row for $(0.6,0.3,1,-1)$ with $\Delta\psi\to-\Delta\psi$.

**Further tests.**

**Vertices of $\mathcal A$**, i.e. ray bisection must hit these:

| a | b | Φ_max | x | Φ_min | x | Ψ_max | x | Ψ_min | x |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 2 | 2 | −2 | 2 | 2 | 2 | −2 | 2 |
| 0.5 | 0 | 2.5 | 3 | −1.5 | 1 | 1.866025 | 1.616025 | −1.866025 | 1.616025 |
| 0.6 | 0.3 | 2.571623 | 3.081623 | −1.264575 | 0.574575 | 2.153256 | 2.093256 | −1.360555 | 0.700555 |
| 0.4 | 0.4 | 2.341641 | 2.581641 | −1.447214 | 0.887214 | 2.341641 | 2.581641 | −1.447214 | 0.887214 |
| 0.9 | 0.09 | 2.897867 | 3.789767 | −1.043589 | 0.135489 | 1.614898 | 0.894898 | −1.134536 | 0.234536 |

**Out of region:** $(a,b,\Phi,\Psi)=(0.4,0.4,1.5,-1.5)$ gives $x=1.659921$ and $K=3.664599$, but there is no $\theta$-band, so the solver must report no orbit. Similarly $(0.55,0.2,-1.3,0.6)$ is invalid.

**Invariances.**

- $a=b=0.4$ with $\Phi+\Psi=0.7$ and any split gives $x=1.967354452441$ and $K=4.232684973306$.
- Tangherlini: any $|\Phi|+|\Psi|<2$ gives $x=2$, $K=4$, $\Delta t=2\pi$ and $\nu_\phi=\nu_\psi=\tfrac12$.

**Closed orbit:** $(0.6,0.3,0.1051348934,-0.1354415960)$ gives $\nu=(1/3,-3/5)$, closing after 15 oscillations.

**Hamiltonian cross-check:** start the full Boyer-Lindquist Hamiltonian at a turning point. $|x-x_0|$ should stay below $10^{-7}$ for 3 oscillations and exceed $10^{-3}x_0$ within about 4–7 oscillations. This is expected; it is the instability.

---

## Appendix A. Pseudocode summary for the implementation

```
// mu = 1, E = 1, a,b with |a|+|b| < 1, controls (Phi, Psi)
xp = 0.5*(1 - a2 - b2 + sqrt((1 - a2 - b2)^2 - 4 a2 b2))
P(x) = x^3 + (a2+b2) x^2 + ((a2-b2)(Phi^2-Psi^2) + a2 b2 + a2 + b2 + 2aPhi + 2bPsi) x
       + (ab + bPhi + aPsi)^2 + (a2-b2)(b2 Phi^2 - a2 Psi^2)
D(x) = x^2 + (a2+b2-1) x + a2 b2
F(x) = P'(x) D(x) - P(x) D'(x)              // < 0 below the orbit, > 0 above
x0   = bisect F on (xp*(1+1e-12), hi)       // hi doubled from 4 until F(hi) > 0
K    = P(x0)/D(x0)
band = thetaBand(Q, Phi, Psi)               // R = Q / y^[Phi==0] / (1-y)^[Psi==0]
if (!band) -> outside region: clamp along ray (bisection on s)
E_   = 1 + a Phi/(x0+a2) + b Psi/(x0+b2); Dl = D(x0)
Aphi = a (x0+b2) E_/Dl + (a2-b2) Phi/(x0+a2)
Apsi = b (x0+a2) E_/Dl - (a2-b2) Psi/(x0+b2)
At   = (x0+a2)(x0+b2) E_/Dl
G(y) = (b2-a2)(y+y1+y2) - (b2 - 2a2 - K)
y(chi) = ybar - dy cos(chi)
d/dchi [sigma, lambda, t, phi, psi] = [1, rho2, rho2+At, Phi/y - Aphi, Psi/(1-y) - Apsi] / (2 sqrt(G(y)))
   rho2 = x0 + a2 + (b2-a2) y
   (subtract the Phi/y, Psi/(1-y) spikes analytically, see B3)
tabulate chi in [0, 2pi] with Simpson, N = 256..512; extend periodically
draw: X1 = sqrt(x0+a2) s cos(phi), Y1 = sqrt(x0+a2) s sin(phi), Z = sqrt(x0+b2) c, hue = psi/2pi
   s = sqrt(y) or sqrt(y2) sin(chi/2) if Phi==0; c = sqrt(1-y) or sqrt(1-y1) cos(chi/2) if Psi==0
```

## Appendix B. Formulas for a general $\mu$

The cleanest route is to work in units of $\mu=1$ and rescale: $x\to\mu x$, $(a,b,\Phi,\Psi,r,t,\lambda)\to\sqrt\mu\,(\cdot)$, $K\to\mu K$, $\sigma\to\sigma/\sqrt\mu$. For example, $\Phi_{\max}=\sqrt\mu+\sqrt{(\sqrt\mu+a)^2-b^2}$.

## Appendix C. Verification scripts

The scripts, their outputs, `region_maps.png` and `figures_repro.png` are copied into [`scripts/`](scripts/) beside this report. They were written and run in a scratch folder that no longer exists. The paper's source and PDF and the reference papers (FS03, KR12, DKLR14, CMBWZ09) were not copied; download them from arXiv. The figure-digitising scripts (`axes.py`, `calib.py`, `digitize.py`, `extract_imgs.py`, `fit_figs.py`, `inspect_figs.py`) need the paper's PDF and its extracted figures. Python was run via `uv`, e.g. `uv run --with sympy --with numpy --with scipy --with mpmath python tables.py`.

| Script | Purpose |
|---|---|
| `check_symbolic.py` (→ `.out`) | Inverse metric, Hamilton-Jacobi separation, (2.3), (3.3)–(3.8), cubic coefficients, horizon and $\Omega$ sign, flat-space embedding, symmetry, $Q(y)$, Fourier sequences, the $a=b$ maximum, §3.4 coefficients and solutions |
| `check_F_grouped.py` | Closed form of the quartic $F$; the boundary formulas for $\Phi^2$, $\Psi^2$ |
| `mp5d.py` | Core numerics: horizon, $P$, $\Delta$, $F$, $Q$, band, solver, Mino constants, elliptic period |
| `integrate_check.py` (→ `.out`) | Full Boyer-Lindquist Hamiltonian integration against the Mino-χ formulation |
| `lyapunov.py` | How long the full integration stays on the orbit |
| `explore_region.py`, `negE.py`, `uniqueness_random.py` | Uniqueness of $(x,K)$; $dt/d\sigma>0$; random test with $2\times10^5$ samples |
| `region.py` (→ `region.out`, `region_maps.png`), `xextreme.py` | Region shape, vertices, star-shapedness and convexity, maps; where $x$ is extremal |
| `figures.py` (→ `figures.out`, `figures_repro.png`) | Reproduce Figs 4.2–4.12; solve Fig 4.12 |
| `fit_figs.py`, `digitize.py`, `axes.py`, `calib.py`, `extract_imgs.py`, `inspect_figs.py` (→ `fit_figs.out`, `figs_fitted_3d.png`, `figs/`) | Raster extraction, digitisation and initial-condition fits for Figs 4.9 and 4.11; colour-wheel sampling |
| `reference_impl.py` (→ `.out`) | JavaScript-portable reference implementation plus accuracy table |
| `tables.py` (→ `.out`) | Unit-test tables, frequencies, $a=b$ closed forms, Fig 4.12 quadrants, §3.4 inequality |
| `resonance.py` | Rank of the frequency map; closed-orbit Newton search |
| `check_jacobi.py` | Jacobi $\mathrm{sn}^2$ closed form and period |
