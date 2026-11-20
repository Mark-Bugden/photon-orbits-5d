# Equations, as implemented

This is every equation the demo's physics uses, in the form the code evaluates it, for review before any drawing work starts. The reference implementation is the Python package in [`src/photon_orbits/`](../src/photon_orbits/):

- [`orbit.py`](../src/photon_orbits/orbit.py): the horizon, the radius, $K$ and the latitude band;
- [`oscillation.py`](../src/photon_orbits/oscillation.py): one latitude oscillation, periodic extension and the picture;
- [`region.py`](../src/photon_orbits/region.py): the allowed region of momenta.

The demo's JavaScript, [`demo/src/physics/`](../demo/src/physics/), is a line-by-line port. It reproduces the Python results in [`tests/fixtures/orbits.json`](../tests/fixtures/orbits.json) to a relative $3\times10^{-16}$.

The derivations and the errata are in the [paper review](research/paper-review.md) (Part A for the errata, Part B for the method). This document follows Part B, and says where the implementation departs from it.

## What to review

Most of this is the paper's own content, re-expressed. Five things are new, and need the closest look:

1. **The radius is solved in $u = x - x_+$** (§3), from a decomposition of $P$ that keeps full precision next to the horizon. The research solver worked in $x$. At an extremal hole that cut the allowed region short by about $3\times10^{-4}$, and at $a = b = \tfrac12$ it could divide by zero.
2. **At exact extremality the co-rotating edge of the allowed region is the segment $\Phi + \Psi = -1$** (§8). Along that segment the orbits close on the horizon. Neither the paper nor the review has this. The plan said only that "the co-rotating corner touches the horizon".
3. **The upper edge of the θ-band is found as $\cos^2\theta_{\max}$** (§4), by the symmetry (2.2). This keeps it exact when $\Psi$ is tiny.
4. **The axis-spike subtraction is rewritten without cancellation** (§6), and the closed-form axis integral is evaluated from each node's index, so its values at $\chi = 0, \pi, 2\pi$ are exact.
5. **Coordinate time $t$ is not computed.** The clock is $\lambda$ throughout, as the plan decided.

**Added after the review**, for the drawing (milestone 2): the state between two nodes, `between` (§6, "Between the nodes").

## 1. Conventions

**Units.** $\mu = 1$ and $E = 1$, so $\Phi$, $\Psi$ and $K$ stand for $\Phi/E$, $\Psi/E$ and $K/E^2$.

**Coordinates.**

- $x = r^2$.
- $\theta \in [0, \pi/2]$, with $y = \sin^2\theta$ and $w = 1 - y = \cos^2\theta$.
- $\phi$ turns in the plane of $a$, and $\psi$ in the plane of $b$.

**The metric** is eq. (2.1) with its typo corrected ($(x+b)^2 \to x + b^2$):

$$
ds^2 = -dt^2 + \frac{1}{\rho^2}\big(dt + a\sin^2\theta\,d\phi + b\cos^2\theta\,d\psi\big)^2 + \frac{\rho^2}{4\Delta}dx^2 + \rho^2 d\theta^2 + (x+a^2)\sin^2\theta\,d\phi^2 + (x+b^2)\cos^2\theta\,d\psi^2,
$$

$$
\rho^2 = x + a^2\cos^2\theta + b^2\sin^2\theta = x + a^2 + (b^2 - a^2)\,y, \qquad \Delta = (x+a^2)(x+b^2) - x.
$$

**Constants of motion.** $p_t = -1$, $p_\phi = \Phi$, $p_\psi = \Psi$, and the Carter-like $K$. Equations (2.3), with

$$
\mathcal E = 1 + \frac{a\Phi}{x+a^2} + \frac{b\Psi}{x+b^2},
$$

give the radial function

$$
\mathcal X = \Delta\Big[x + (a^2-b^2)\Big(\frac{\Phi^2}{x+a^2} - \frac{\Psi^2}{x+b^2}\Big) - K\Big] + (x+a^2)(x+b^2)\,\mathcal E^2
$$

and the latitude function $\Theta = a^2\cos^2\theta + b^2\sin^2\theta + K - \Phi^2/\sin^2\theta - \Psi^2/\cos^2\theta$.

**Signs.** For $a > 0$ the hole turns towards $-\phi$, so co-rotating photons have $\Phi < 0$. Likewise for $b$ and $\Psi$.

**Spins.** The module accepts any $a, b$ with $|a| + |b| \le 1$; the demo uses $a, b \in [0, 0.99]$. The corners in §8 assume $a, b \ge 0$.

## 2. The horizon

`horizons(a, b)` computes

$$
s = 1 - a^2 - b^2,\qquad d = s^2 - 4a^2b^2 = \big(1-(a+b)^2\big)\big(1-(a-b)^2\big),\qquad e = \sqrt{d},\qquad x_\pm = \tfrac12(s \pm e).
$$

- There is no horizon if $s \le 0$ or $d < -10^{-12}$.
- A $d$ in $[-10^{-12}, 0)$ is an extremal hole after rounding, and is clamped to $0$.
- $e = x_+ - x_-$, so $\Delta(x_+ + u) = u(u + e)$. $e = 0$ exactly when the hole is extremal.
- $a = 1, b = 0$ gives $s = 0$, i.e. $x_+ = 0$, and counts as no horizon.

## 3. The radius and $K$

### The equations

Since $K$ enters $\mathcal X$ only through $-K\Delta$,

$$
\mathcal X = P - K\Delta,\qquad P = \Delta W + (x+a^2)(x+b^2)\,\mathcal E^2,\qquad W = x + (a^2 - b^2)\Big(\frac{\Phi^2}{x+a^2} - \frac{\Psi^2}{x+b^2}\Big).
$$

$P$ is the monic cubic of the review (§B2). `p_coefficients` keeps it in that form for the tests:

$$
P = x^3 + (a^2+b^2)x^2 + \big[(a^2-b^2)(\Phi^2-\Psi^2) + a^2b^2 + a^2 + b^2 + 2a\Phi + 2b\Psi\big]x + (ab + b\Phi + a\Psi)^2 + (a^2-b^2)(b^2\Phi^2 - a^2\Psi^2).
$$

A spherical orbit needs $\mathcal X = \mathcal X' = 0$. With $V = P/\Delta$, that is $V'(x_0) = 0$ and $K = V(x_0)$. Equivalently, $x_0$ is a root of the quartic $F = P'\Delta - P\Delta' = \Delta^2 V'$.

### In $u = x - x_+$

`shifted_coefficients` writes $P$ and $F$ in $u$. At $x_+$, $\Delta = 0$, so the Taylor coefficients of $P$ there come straight from the decomposition above:

$$
\begin{aligned}
P_0 &= (x_++a^2)(x_++b^2)\,\mathcal E_+^2,\\
P_1 &= e\,W_+ + (2x_+ + a^2 + b^2)\,\mathcal E_+^2 - 2\,\mathcal E_+\Big[a\Phi\frac{x_++b^2}{x_++a^2} + b\Psi\frac{x_++a^2}{x_++b^2}\Big],\\
P_2 &= 3x_+ + a^2 + b^2,\qquad P_3 = 1,
\end{aligned}
$$

where $\mathcal E_+$ and $W_+$ are $\mathcal E$ and $W$ at $x_+$. The bracket in $P_1$ is $-(x+a^2)(x+b^2)\,\mathcal E'$ at $x_+$. With $\Delta = u^2 + eu$, this gives

$$
F = u^4 + 2e\,u^3 + (eP_2 - P_1)\,u^2 - 2P_0\,u - eP_0 .
$$

At an extremal hole ($e = 0$) this is $F = u\,(u^3 - P_1 u - 2P_0)$, a depressed cubic times $u$.

### Why $u$, not $x$

At an extremal hole $x_+$ is a double root of $\Delta = x^2 + (a^2+b^2-1)x + a^2b^2$. Rounding blurs the position of a double root to about $\sqrt{\epsilon_{\rm mach}} \approx 10^{-8}$. Near $\mathcal E_+ = 0$, $F$ also has a near-quadruple root at $x_+$, so its sign is rounding noise within about $3\times10^{-5}$ of the horizon.

The research solver started its bracket at $x_+(1 + 10^{-12}) + 10^{-14}$, inside that noise. In $u$, every coefficient carries exactly the power of $\mathcal E_+$ it should, and the noise goes away.

### The solver

`radial` works as follows.

1. If $P_0 = 0$ ($\mathcal E_+ = 0$), there is **no orbit ("at horizon")**. $P$ and $\Delta$ both vanish at the horizon, and $V$ is $0/0$ there. This is the plan's first guard.
2. Otherwise $F < 0$ just above the horizon ($V \to +\infty$ there) and $F \to +\infty$. So bisect $F$ on $u \in (0, h]$, with $h = 4$ doubled until $F(h) > 0$. Stop when the interval's width is below $10^{-15}h$, or after 200 steps.
3. $x_0 = x_+ + u$, $\Delta(x_0) = u(u + e)$, and $K = P(u)/\Delta(x_0)$.
4. **The second-root check.** Divide $F$ by $(u - u_0)$. The cubic $C$ that remains is positive just above $u = 0$ and at $+\infty$. Count its sign changes between its turning points in $u > 0$; these come in pairs, a maximum of $V$ and a second minimum. A nonzero count raises a `SecondRootWarning` (a `console.warn` in JS) and is stored on the orbit. In 3,000 random samples (`test_one_orbit_per_pair_of_momenta`), numpy's roots of $F$ agree that there is exactly one root above the horizon.

## 4. The latitude band

### The test

$$
Q(y) = y(1-y)\,\Theta = (a^2-b^2)\,y^3 + (b^2-2a^2-K)\,y^2 + (\Phi^2-\Psi^2+a^2+K)\,y - \Phi^2 .
$$

$Q(0) = -\Phi^2 \le 0$ and $Q(1) = -\Psi^2 \le 0$, so a cubic $Q$ is non-negative on at most one interval of $(0,1)$: the band $[y_1, y_2]$.

`band_peak` builds $R$, which has the sign of $\Theta$:

- $R = Q/y$ if $\Phi = 0$;
- $R = Q/(1-y)$ if $\Psi = 0$.

These keep $\Theta$ finite on an axis the orbit can then reach. The division by $(1-y)$ works up from the constant term, $(1-y)(s_0 + s_1 y + \dots) = c_0 + c_1y + \dots$, so the constant term of the quotient stays exactly $-\Phi^2$.

The largest value of $R$ on $[0,1]$ is taken over three kinds of candidate:

- the roots of $R'$ in $(0,1)$, by the stable quadratic formula (for $a = b$, $R'$ is linear);
- $y = 0$ if $\Phi = 0$;
- $y = 1$ if $\Psi = 0$.

Call it $R_m$, at $y_m$, and let $\tau = 10^{-12}(1 + |K| + \Phi^2 + \Psi^2)$. Then:

- $R_m < -\tau$: **no orbit ("no band")**, since $\theta$-motion is forbidden;
- $|R_m| \le \tau$: a circular orbit at constant $\theta$, with $y_1 = y_2 = y_m$;
- otherwise, find the edges.

This replaces §3.1 of the paper. In particular it does not use the $a = b$ criterion with the sign error (erratum 3).

### The edges

- $y_1 = \sin^2\theta_{\min}$ is the root of $R$ in $(0, y_m)$, or $0$ if $R(0) \ge 0$, which happens only when $\Phi = 0$.
- $w_2 = \cos^2\theta_{\max}$ is the root in $(0, 1-y_m)$ of $\tilde R(w) = R(1-w)$, or $0$ if $\tilde R(0) \ge 0$, which happens only when $\Psi = 0$. By the paper's symmetry (2.2), $Q(1-w)$ is $Q$ with $(a, \Phi) \leftrightarrow (b, \Psi)$, so $\tilde R$ is built by the same code with the arguments swapped.

Both are found by bisection, until the interval's width is below $10^{-15}$ of its upper end. Then $y_2 = 1 - w_2$.

Finding $w_2$ directly matters when $\Psi$ is small. Then $\cos^2\theta_{\max} \approx \Psi^2/(b^2+K)$, which for $\Psi = 10^{-10}$ is far below the spacing of floats near 1 ($1.1\times10^{-16}$); as $1 - y_2$ it would round to 0, and the orbit would divide by zero at the equator.

## 5. Mino time and the angle variable

With Mino time $d\sigma = d\lambda/\rho^2$, equations (2.3) become

$$
\Big(\frac{dy}{d\sigma}\Big)^2 = 4Q(y),\quad \frac{d\phi}{d\sigma} = \frac{\Phi}{y} - A_\phi,\quad \frac{d\psi}{d\sigma} = \frac{\Psi}{1-y} - A_\psi,\quad \frac{d\lambda}{d\sigma} = \rho^2,
$$

with constants of the orbit, computed in `solve_orbit`, and $\Delta = \Delta(x_0) = u(u+e)$:

$$
A_\phi = \frac{a(x+b^2)\,\mathcal E}{\Delta} + \frac{(a^2-b^2)\,\Phi}{x+a^2},\qquad
A_\psi = \frac{b(x+a^2)\,\mathcal E}{\Delta} - \frac{(a^2-b^2)\,\Psi}{x+b^2}.
$$

These are (2.3d, e) divided through by $\rho^2\dot\sigma = 1$.

**Factorisation.** On the band,

$$
Q(y) = (y-y_1)(y_2-y)\,G(y),\qquad G(y) = (b^2-a^2)(y + y_1 + y_2) - (b^2 - 2a^2 - K) > 0 .
$$

$G$ is linear: $-(b^2-a^2)$ times $(y - y_3)$, with $y_3$ the third root of $Q$. It holds for $a = b$ too, where $G = a^2 + K$. `solve_orbit` stores $G_1 = G(y_1)$ and $G_2 = G(y_2)$.

**The angle.** With $\bar y = \tfrac12(y_1+y_2)$ and $\Delta y = \tfrac12(y_2-y_1)$, set

$$
y = \bar y - \Delta y\cos\chi = y_1 + (y_2 - y_1)\sin^2\tfrac\chi2,\qquad 1 - y = w_2 + (y_2-y_1)\cos^2\tfrac\chi2 .
$$

The code always uses the right-hand forms, so $y$ near the pole and $1-y$ near the equator keep full relative precision.

Then $(y-y_1)(y_2-y) = \Delta y^2\sin^2\chi$ and $dy/d\chi = \Delta y\sin\chi$, so $d\sigma/d\chi = 1/(2\sqrt{G})$, and

$$
\frac{d(\sigma, \lambda, \phi, \psi)}{d\chi} = \frac{1}{2\sqrt{G(y)}}\Big(1,\ \rho^2,\ \frac{\Phi}{y} - A_\phi,\ \frac{\Psi}{1-y} - A_\psi\Big).
$$

$\chi$ increases monotonically, with no turning points. $\chi = 0$ is $\theta_{\min}$ and $\chi = \pi$ is $\theta_{\max}$.

## 6. One oscillation

### The quadrature

`tabulate` integrates $\chi$ from $0$ to $2\pi$ by Simpson's rule on $N$ steps, from $0$ at $\theta_{\min}$:

$$
Z_{i+1} = Z_i + \tfrac{h}{6}\big[f(\chi_i) + 4f(\chi_i + \tfrac h2) + f(\chi_{i+1})\big],\qquad h = 2\pi/N,\quad \chi_i = ih .
$$

This is the same as RK4, since the right-hand side depends on $\chi$ alone. $N = 256$, or $512$ when $\min(y_1, w_2) < 10^{-4}$, i.e. when the band comes within about 0.6° of an axis.

### The axis terms

$\Phi/(2y\sqrt G)$ is sharply peaked at $\chi = 0$ when $y_1$ is small, and $\Psi/(2(1-y)\sqrt G)$ is sharply peaked at $\chi = \pi$ when $w_2$ is small. Each splits into a part integrated exactly and a smooth remainder. Because $G - G_1 = (b^2-a^2)(y - y_1)$, the remainder has a form with no cancellation:

$$
\frac{\Phi}{2y\sqrt G} = \frac{\Phi}{2y\sqrt{G_1}} - \frac{\Phi\,(b^2-a^2)\,(y-y_1)}{2y\sqrt G\sqrt{G_1}\,\big(\sqrt G + \sqrt{G_1}\big)},
$$

$$
\frac{\Psi}{2(1-y)\sqrt G} = \frac{\Psi}{2(1-y)\sqrt{G_2}} + \frac{\Psi\,(b^2-a^2)\,(y_2-y)}{2(1-y)\sqrt G\sqrt{G_2}\,\big(\sqrt G + \sqrt{G_2}\big)} .
$$

Here $(y - y_1)/y$ and $(y_2-y)/(1-y)$ both lie in $[0, 1)$. The review subtracted $\Phi/(2y\sqrt{G_1})$ directly. That is the same function, but it loses precision when $y_1$ is tiny.

`rates` gives the smooth remainders. The first terms integrate exactly:

$$
\int_0^\chi \frac{d\chi'}{y} = I(y_1, y_2; \chi),\qquad \int_0^\chi \frac{d\chi'}{1-y} = I(1-y_1, w_2; \chi),
$$

$$
I(p, q; \chi) = \int_0^\chi\frac{d\chi'}{p + (q-p)\sin^2(\chi'/2)} = \frac{2}{\sqrt{pq}}\arctan\Big(\sqrt{q/p}\,\tan\frac\chi2\Big)
$$

on the branch continuous in $\chi$, with $I(p, q; 2\pi) = 2\pi/\sqrt{pq}$. These are multiplied by $\Phi/(2\sqrt{G_1})$ and $\Psi/(2\sqrt{G_2})$ and added at the nodes.

`axis_integral` evaluates $I$ at node $i$ of $N$ from the integer index:

- for $4j \le N$: $\arctan(k\tan(\pi j/N))$;
- otherwise: $\pi/2 - \arctan\big(\tan(\pi(N-2j)/2N)/k\big)$;
- for $i > N/2$: $I(2\pi - \chi) = I(2\pi) - I(\chi)$.

This makes the values at $\chi = 0$, $\pi$ and $2\pi$ exact. The reason is the size of the spike. With $\Phi = 10^{-12}$, $\phi$ swings by about $\pi$ within $\Delta\chi \sim 10^{-12}$ of the pole. Evaluating at the floating-point $2\pi$, which is out by $4\times10^{-16}$, gave a $10^{-4}$ error in $\Delta\phi$.

If $\Phi = 0$ there is neither an axis term nor a remainder in $\phi$; likewise for $\Psi = 0$ and $\psi$.

### Longer paths

**Periodic extension.** `node(table, i)` gives the state at any node, with $m = \lfloor i/N\rfloor$:

$$
(\sigma, \lambda, \phi, \psi)(\chi + 2\pi m) = \text{table}(\chi) + m\,(\Delta\sigma, \Delta\lambda, \Delta\phi, \Delta\psi).
$$

**Frequencies.** The turns per oscillation count an axis passage as half a turn (review §B6), which makes them continuous mod 1 through $\Phi = 0$ and $\Psi = 0$:

$$
\nu_\phi = \frac{\Delta\phi + \pi\,[\Phi = 0,\ y_1 = 0]}{2\pi},\qquad
\nu_\psi = \frac{\Delta\psi + \pi\,[\Psi = 0,\ w_2 = 0]}{2\pi}.
$$

**Between the nodes.** Near an axis, $\phi$ or $\psi$ swings by about $\pi$ within one step, so the drawing needs values between the nodes. `between(orbit, table, i, t)` gives the state at $\chi_i + th$, for $0 \le t \le 1$: one Simpson step of width $th$ from node $i$ for the smooth remainders, plus the closed-form axis terms,

$$
Z(\chi_i + th) = Z_i + \tfrac{th}{6}\big[f(\chi_i) + 4f(\chi_i + \tfrac{th}{2}) + f(\chi_i + th)\big] + c\,\big[I(i + t) - I(i)\big],
$$

where $c\,I$ is the axis term of §6, evaluated by `axis_integral` at the fractional index $i + t$. At $t = 0$ this is node $i$ exactly, and at $t = 1$ it meets node $i + 1$ to rounding. Next to a swing, the smooth remainder varies on the swing's own scale, so one step is good to about $10^{-8}$, like the table (`test_between_follows_the_swing_near_an_axis`).

The drawing (`demo/src/trail.js`) takes 128 steps per oscillation and halves any step in which $\phi$ or $\psi$ turns by more than 10°, as often as needed: a swing $10^{-12}$ wide takes about 40 halvings. Oscillation $m$ is then oscillation 0 turned about the axis by $2\pi m\,\nu_\phi$, with $Z$ changing sign each oscillation for an orbit through the equator, and $\psi$ advanced by $m\,\Delta\psi$.

**What is not integrated.** $x$ is held fixed. Every orbit is radially unstable; see the Hamiltonian cross-check in §10. Coordinate time $t$ is not computed.

## 7. The picture

`embed` gives the paper's projection without the constant twist $\arctan(a/r)$, which at fixed $r$ is a rigid rotation:

$$
X = \sqrt{x+a^2}\,\sin\theta\cos\phi,\qquad Y = \sqrt{x+a^2}\,\sin\theta\sin\phi,\qquad Z = \sqrt{x+b^2}\,\cos\theta,\qquad \text{colour} = \psi .
$$

`latitude` gives $\sin\theta$ and $\cos\theta$ from the absolute $\chi$, which is not reduced mod $2\pi$. In general

$$
\sin\theta = \sqrt{y_1 + (y_2-y_1)\sin^2\tfrac\chi2},\qquad \cos\theta = \sqrt{w_2 + (y_2-y_1)\cos^2\tfrac\chi2}.
$$

Where an orbit passes through an axis, the coordinates are signed, as the paper's figures do with their $\theta \in [0, \pi]$ (plan decision):

- **$\Phi = 0$ and $y_1 = 0$:** $\sin\theta = \sqrt{y_2}\,\sin(\chi/2)$. It changes sign each oscillation, so the path goes straight over the pole.
- **$\Psi = 0$ and $w_2 = 0$:** $\cos\theta = \sqrt{1-y_1}\,\cos(\chi/2)$. $Z$ is negative for $\chi \in (\pi, 3\pi)$, so the path crosses the equator ring onto the lower half and keeps its colour.

**The side panel** plots the same signed coordinates, $\theta = \operatorname{atan2}(\sin\theta, \cos\theta)$, against $\lambda$, with $\phi$ and $\psi$ unwrapped (Mark's decision in milestone 4: the curves are geodesics, so they stay smooth). So $\theta$, $\phi$ and $\psi$ are continuous through the axes:

- through the pole, $\theta$ swings between $-\theta_{\max}$ and $\theta_{\max}$, and $\phi$ advances by $\Delta\phi$ per oscillation, without the half turn that the picture's rotation adds;
- through the equator, $\theta$ runs from $\theta_{\min}$ past $\pi/2$ to $\pi - \theta_{\min}$ and back, over two oscillations;
- through both ($\Phi = \Psi = 0$, so $y_1 = 0$ and $y_2 = 1$), $\theta = \chi/2$, which grows by $\pi$ per oscillation.

`trace` in `demo/src/trail.js` gives these from the drawn oscillation: for oscillation $m$, $\theta \mapsto (-1)^m\theta$ through the pole, $\theta \mapsto \pi - \theta$ for odd $m$ through the equator, and $\theta \mapsto \theta + m\pi$ through both. `test/trail.test.js` checks that the signed $\theta$ and unwrapped $\phi$ put every point where the picture does.

## 8. The allowed region

### Corners and edge

**Corners.** `vertices` gives the paper's circular orbits (3.9)–(3.15), for $a, b \ge 0$, with each square root clamped at 0:

| Corner | Momentum | $x$ |
| --- | --- | --- |
| $\Psi = 0$, counter-rotating | $\Phi = 1 + \sqrt{(1+a)^2 - b^2}$ | $1 + a - b^2 + \sqrt{(1+a)^2-b^2}$ |
| $\Psi = 0$, co-rotating | $\Phi = -1 - \sqrt{(1-a)^2 - b^2}$ | $1 - a - b^2 + \sqrt{(1-a)^2-b^2}$ |
| $\Phi = 0$, counter-rotating | $\Psi = 1 + \sqrt{(1+b)^2 - a^2}$ | $1 + b - a^2 + \sqrt{(1+b)^2-a^2}$ |
| $\Phi = 0$, co-rotating | $\Psi = -1 - \sqrt{(1-b)^2 - a^2}$ | $1 - b - a^2 + \sqrt{(1-b)^2-a^2}$ |

The clamp matters at extremality. For $a = 0.9$, $b = 0.1$, $(1-a)^2 - b^2$ evaluates to $-6.9\times10^{-18}$.

**The edge.** `edge(a, b, dΦ, dΨ)` finds how far the region reaches along a unit direction. It bisects on $s \in [0, 4]$, with 4 doubled while still inside, for 48 steps. The test at each step is `orbit_exists`: the radius, then the band's peak, which is the same decision as `solve_orbit` makes.

- `boundary` uses 96 rays, $\alpha_k = 2\pi k/96$.
- The directions are passed as unit vectors, so the JS reproduces the same edge to the last bit.
- The tests check that the region is star-shaped and convex.
- Away from extremality, the edge is made of circular orbits, as in Fig 4.12.

### Extremal holes

At $a + b = 1$ we have $x_+ = ab$, $x_+ + a^2 = a$ and $x_+ + b^2 = b$, so

$$
\mathcal E_+ = 1 + \Phi + \Psi,
$$

and $\Omega_a = \Omega_b = 1/(a+b) = 1$. On the segment $\Phi + \Psi = -1$ ($\Phi, \Psi \le 0$), co-rotating photons sit exactly at the superradiant bound $E = \Omega_a|\Phi| + \Omega_b|\Psi|$.

Write $\delta = b\Phi + a\Psi$, which is negative there. Then $P_0 = ab\,\mathcal E_+^2$ and $P_1 = -2\delta\,\mathcal E_+ + O(\mathcal E_+^2)$. The orbit is the positive root of $u^3 - P_1u - 2P_0 = 0$, which gives two regimes:

- **Inside ($\mathcal E_+ \to 0^+$):** $x - x_+ \approx \sqrt{2|\delta|\,\mathcal E_+}$ and $K \to P_2 = 1 + ab$. The θ-band stays wide (about 25°–60° at $\Phi = \Psi = -\tfrac12$). **So the orbits close on the horizon while still sweeping a band of latitudes.**
- **Outside ($\mathcal E_+ < 0$):** the root is $x - x_+ \approx ab\,|\mathcal E_+|/|\delta|$, with $K \to 1 + ab - \delta^2/(ab)$, and $\Theta < 0$ everywhere: no orbit.

So **the co-rotating edge of the region at extremality is the straight segment from $(-1, 0)$ to $(0, -1)$**, not a curve of circular orbits. On the segment itself, $V$ is $0/0$ at the horizon, and `radial` returns "at horizon".

The tests check all of this for five extremal holes:

- the edge lies on $\Phi + \Psi = -1$ to $10^{-8}$; that is as extremal as floating-point spins get;
- the $\sqrt{\mathcal E_+}$ approach;
- the limits of $K$;
- that there is no orbit beyond the line.

As $\mathcal E_+ \to 0$, $A_\phi$ and $A_\psi$ grow like $1/(x - x_+)$, so near this edge the orbit is dragged round with the horizon:

| $a$, $b$ | $1 + \Phi + \Psi$ | turns of $\phi$ per oscillation |
| --- | --- | --- |
| 0.5, 0.5 | $10^{-4}$ | about 21 |
| 0.7, 0.3 | $10^{-5}$ | about 68 |

The drawing subdivides by $\phi$ as well as $\psi$ there (§6, "Between the nodes"). The pad (`demo/src/momenta.js`) stops the dot at $1 + \Phi + \Psi = 10^{-4}$, where $\phi$ turns about 21 times per oscillation.

## 9. Guards and tolerances

| Where | What | Why |
| --- | --- | --- |
| Every entry point | $\lvert\Phi\rvert, \lvert\Psi\rvert < 10^{-12}$ become exactly $0$ | The plan's second guard. A ray straight up the $\Psi$ axis gives $\Phi = s\cos(\pi/2) \approx 10^{-16}$. |
| `horizons` | $d \in [-10^{-12}, 0)$ clamped to $0$; no horizon if $s \le 0$ or $d < -10^{-12}$ | Extremal holes after rounding |
| `radial` | $\mathcal E_+ = 0$: no orbit, "at horizon" | The plan's first guard: $V = 0/0$ at the horizon |
| `radial` | Bisection to $10^{-15}$ relative, at most 200 steps | |
| `band` | $\tau = 10^{-12}(1 + \lvert K\rvert + \Phi^2 + \Psi^2)$ separates no band, circular and banded | |
| `band` | Edges by bisection to $10^{-15}$ relative | Tiny $y_1$, $w_2$ keep full precision |
| `vertices`, `horizons` | Square roots clamped at $0$ | |
| `tabulate` | $N = 256$, or $512$ if $\min(y_1, w_2) < 10^{-4}$ | The plan's rule |
| `radial` | A second root above the horizon raises a warning | Never seen |

## 10. How it is checked

| Test | What it shows |
| --- | --- |
| `test_orbit.py::test_the_review_table`, `test_oscillation.py::test_the_review_table` | All 15 rows of the review's table (§B7): $x_+$, $x$, $r$, $K$, $\theta_{\min}$, $\theta_{\max}$, $\Delta\sigma$, $\Delta\lambda$, $\Delta\phi$, $\Delta\psi$, $\nu$, to every printed digit |
| `test_the_circular_orbits_at_45_degrees` | Fig 4.12 and its three companions |
| `test_the_shifted_coefficients_are_p_and_f` | §3's $P(u)$ and $F(u)$ equal the paper's $P$ and $P'\Delta - P\Delta'$ |
| `test_one_orbit_per_pair_of_momenta` | One root of $F$ above the horizon in 3,000 random samples, matching numpy's roots |
| `test_momenta_outside_the_region_have_no_orbit` | The review's out-of-region cases, including its counterexample to §3.1.1 |
| `test_equal_spins_…`, `test_tangherlini_…` | $x$, $K$ depend only on $\Phi + \Psi$ when $a = b$; $x = 2$, $K = 4$, great circles for $a = b = 0$ |
| `test_the_period_in_mino_time_…` | $\Delta\sigma = 2K(m)/\sqrt{G_{\max}}$, $m = 1 - G_{\min}/G_{\max}$, to $10^{-11}$ |
| `test_equal_spins_have_closed_forms` | $\Delta\sigma$, $\Delta\phi$, $\Delta\psi$ in closed form for $a = b$, to $10^{-11}$ |
| `test_the_closed_orbit_…` | $\nu = (1/3, -3/5)$: closes after 15 oscillations |
| `test_the_table_has_converged` | $N = 256/512$ against $N = 4096$, to $10^{-7}$ (the plan's target), including band edges within $10^{-8}$ of an axis |
| `test_between_the_nodes`, `test_between_repeats_with_the_table`, `test_between_follows_the_swing_near_an_axis` | Between the nodes: a table 16 times finer to $10^{-7}$ on 11 orbits, node $i$ exactly at $t = 0$ and node $i + 1$ to $10^{-12}$ at $t = 1$, periodic extension, and scipy's adaptive quadrature of the full rate across the swings near the pole and the equator, to $10^{-7}$ |
| `test_small_momenta_…`, `test_tiny_momenta_count_as_zero` | $\nu$ continuous mod 1 through $\Phi = 0$ and $\Psi = 0$ for momenta from $10^{-12}$ to $10^{-4}$; $10^{-16}$ is exactly $0$ |
| `test_phi_zero_…`, `test_psi_zero_…`, `test_the_path_lies_on_the_spheroid…` | Through the pole, across the equator, on the spheroid, and continuous |
| `test_region.py` | Corners, star shape, convexity, circular-orbit edges, the extremal edge and its limits, and the plan's stress test: 96 directions × 8 fractions on 7 holes, all valid. Up to 0.999 of the way to the edge, every orbit is at least $1.9\times10^{-3}$ outside the horizon; the plan quoted 0.002 from the research solver. |
| `test_hamiltonian.py` | Null geodesics of the corrected (2.1), from an inverse metric built by sympy and Hamilton's equations by DOP853 at rtol $10^{-12}$, from $\theta_{\min}$ over 3 oscillations. The table matches to $2\times10^{-8}$ in $\theta$, $\phi$, $\psi$, and to $1.2\times10^{-7}$ in $\psi$ for Fig 4.10, where $\psi$ swings by about $\pi$ near the equator. $x$ stays within $2\times10^{-8}$ of $x_0$, and $\lvert H\rvert < 10^{-10}$. That is on six orbits, one of them 0.1 outside an extremal horizon. Left free, $x$ leaves the orbit by 0.1% after 5.8–6.4 oscillations. |
| `demo/test/physics.test.js` | The JS port against `tests/fixtures/orbits.json`: 31 orbits, with points at and between the nodes, 7 no-orbit cases, 8 regions of 96 rays. Largest relative difference $3\times10^{-16}$. An orbit and its table take about 0.13 ms; a 96-ray edge about 7 ms. |

## 11. Changes from the research scripts

| Research (`reference_impl.py`, `mp5d.py`) | Here | Why |
| --- | --- | --- |
| Bisection on $F(x)$ from $x_+(1+10^{-12}) + 10^{-14}$ | Bisection on $F(u)$ from $u = 0$, with $P = \Delta W + (x+a^2)(x+b^2)\mathcal E^2$ | Extremal holes (§3, §8) |
| $K = P(x)/\Delta(x)$, Mino constants with $\Delta = (x+a^2)(x+b^2) - x$ | $\Delta = u(u + e)$ | Precision near the horizon |
| $y_2$ by bisection in $y$ | $w_2 = \cos^2\theta_{\max}$ by bisection in $w$, by symmetry (2.2) | Small $\Psi$ (§4) |
| Division of $Q$ by $(1-y)$ from the top | From the constant term | Keeps $-\Phi^2$ exact |
| Spike subtraction $\frac{\Phi}{2y}\big(\frac1{\sqrt G} - \frac1{\sqrt{G_1}}\big)$, only if $y_1 > 0$ | The cancellation-free remainder, whenever $\Phi \neq 0$ | Small $y_1$ (§6) |
| $I(\chi)$ via `atan2` at floating $\chi$ | From the node index, with symmetries | Exact at $0$, $\pi$, $2\pi$ |
| $t$ tabulated | Not computed | The plan's clock is $\lambda$ |
| Momenta used as given | Snapped to $0$ below $10^{-12}$ | The plan's second guard |
