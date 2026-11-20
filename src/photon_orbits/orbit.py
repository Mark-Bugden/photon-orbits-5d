"""The spherical photon orbit for given spins and momenta: the horizon, the radius and the latitude band.

Units μ = 1 and E = 1, so Φ, Ψ and K stand for Φ/E, Ψ/E and K/E². x = r², y = sin²θ and w = 1 − y = cos²θ,
with θ in [0, π/2]. Polynomials are lists of coefficients, highest power first. docs/equations.md gives every
formula here with its derivation and the reasons for the guards.

demo/src/physics/orbit.js is a line-by-line port, so this module uses plain floats, writes squares as
products, and does the same operations in the same order.
"""

import math
import warnings
from dataclasses import dataclass

SNAP = 1e-12            # a momentum smaller than this in size counts as exactly 0
BAND_TOL = 1e-12        # the band test's tolerance, relative to 1 + |K| + Φ² + Ψ²


class SecondRootWarning(UserWarning):
    """F has another root above the horizon: a second spherical orbit with the same momenta, never yet seen."""


@dataclass(frozen=True)
class Orbit:
    """A spherical photon orbit. y2 = 1 − w2; w2 is kept separately so that it is exact near the equator."""

    a: float
    b: float
    Phi: float
    Psi: float
    x_plus: float           # the horizon, x₊
    x: float                # the orbit's radius, x = r²
    K: float                # the Carter-like constant
    y1: float               # sin²θ_min
    y2: float               # sin²θ_max
    w2: float               # cos²θ_max
    delta: float            # Δ(x)
    calE: float             # 𝓔 = 1 + aΦ/(x + a²) + bΨ/(x + b²)
    A_phi: float            # dφ/dσ = Φ/y − A_φ
    A_psi: float            # dψ/dσ = Ψ/(1 − y) − A_ψ
    G1: float               # G(y1)
    G2: float               # G(y2)
    extra_roots: int        # other roots of F above the horizon (always 0 so far)
    ok = True


@dataclass(frozen=True)
class NoOrbit:
    """No spherical orbit: 'no horizon' (|a| + |b| > 1), 'at horizon' (the orbit would sit on the horizon), or
    'no band' (the radius exists, but θ-motion is forbidden: the momenta are outside the allowed region)."""

    reason: str
    x: float | None = None
    K: float | None = None
    ok = False


def snap(value):
    """A momentum, with sizes below SNAP set to exactly 0."""
    return 0.0 if abs(value) < SNAP else float(value)


def sqrt0(value):
    """The square root, with rounding errors below 0 clamped to 0."""
    return math.sqrt(value) if value > 0 else 0.0


def polyval(c, x):
    total = 0.0
    for coefficient in c:
        total = total * x + coefficient
    return total


def quadratic_roots(A, B, C):
    """The real roots of A t² + B t + C in increasing order; a linear equation if A = 0."""
    if A == 0:
        return [] if B == 0 else [-C / B]
    disc = B * B - 4 * A * C
    if disc < 0:
        return []
    root = math.sqrt(disc)
    q = -0.5 * (B - root if B < 0 else B + root)
    if q == 0:
        return [0.0]
    t1, t2 = q / A, C / q
    return [t1, t2] if t1 <= t2 else [t2, t1]


def horizons(a, b):
    """(x₊, e): the outer horizon x₊ = r₊² and e = x₊ − x₋, so Δ(x₊ + u) = u(u + e); None if there is no
    horizon (|a| + |b| > 1, or x₊ = 0). e = 0 for an extremal hole."""
    s = 1 - a * a - b * b
    disc = s * s - 4 * a * a * b * b
    if s <= 0 or disc < -1e-12:     # a disc just below 0 is an extremal hole after rounding
        return None
    e = sqrt0(disc)
    return 0.5 * (s + e), e


def horizon(a, b):
    """x₊ = r₊², the outer horizon; None if there is none."""
    found = horizons(a, b)
    return None if found is None else found[0]


def p_coefficients(a, b, Phi, Psi):
    """P(x), where 𝒳(x) = P(x) − K Δ(x)."""
    A2, B2 = a * a, b * b
    p1 = (A2 - B2) * (Phi * Phi - Psi * Psi) + A2 * B2 + A2 + B2 + 2 * a * Phi + 2 * b * Psi
    m = a * b + b * Phi + a * Psi
    p0 = m * m + (A2 - B2) * (B2 * Phi * Phi - A2 * Psi * Psi)
    return [1.0, A2 + B2, p1, p0]


def shifted_coefficients(a, b, Phi, Psi, x_plus, e):
    """P and F = P′Δ − PΔ′ as polynomials in u = x − x₊, built so that no coefficient loses precision when
    the orbit is near the horizon.

    P = ΔW + (x + a²)(x + b²)𝓔², with W = x + (a² − b²)(Φ²/(x + a²) − Ψ²/(x + b²)), so its Taylor coefficients
    at x₊, where Δ = 0, are P₀ = (x₊ + a²)(x₊ + b²)𝓔₊², P₁ = eW₊ + (2x₊ + a² + b²)𝓔₊² + 2(x₊ + a²)(x₊ + b²)𝓔₊𝓔′₊
    and P₂ = 3x₊ + a² + b². With Δ = u² + eu this gives F = u⁴ + 2eu³ + (eP₂ − P₁)u² − 2P₀u − eP₀.
    """
    A2, B2 = a * a, b * b
    ra, rb = x_plus + A2, x_plus + B2
    calE = 1 + a * Phi / ra + b * Psi / rb
    drag = a * Phi * rb / ra + b * Psi * ra / rb        # −(x + a²)(x + b²)𝓔′ at x₊
    W = x_plus + (A2 - B2) * (Phi * Phi / ra - Psi * Psi / rb)
    p0 = ra * rb * calE * calE
    p1 = e * W + (2 * x_plus + A2 + B2) * calE * calE - 2 * calE * drag
    p2 = 3 * x_plus + A2 + B2
    return [1.0, p2, p1, p0], [1.0, 2 * e, e * p2 - p1, -2 * p0, -e * p0]


def q_coefficients(a, b, Phi, Psi, K):
    """Q(y) = y(1 − y)Θ. With (a, Φ) and (b, Ψ) swapped it is Q(1 − w) as a polynomial in w."""
    A2, B2 = a * a, b * b
    return [A2 - B2, B2 - 2 * A2 - K, Phi * Phi - Psi * Psi + A2 + K, -Phi * Phi]


def extra_roots(F, root):
    """How many more roots the quartic F (in u) has at u > 0, besides root.

    C = F/(u − root) is positive just above u = 0, where F < 0, and at +∞, so its roots there are the sign
    changes between its turning points; they come in pairs, a maximum of V and another minimum.
    """
    cubic = [1.0]
    for coefficient in F[1:-1]:
        cubic.append(coefficient + root * cubic[-1])
    signs = [1.0] + [polyval(cubic, t) for t in quadratic_roots(3.0, 2 * cubic[1], cubic[2]) if t > 0] + [1.0]
    signs = [s for s in signs if s != 0]
    return sum(1 for s, t in zip(signs, signs[1:]) if (s < 0) != (t < 0))


def radial(a, b, Phi, Psi):
    """(x, K, Δ(x), extra_roots) for the orbit's radius x, the root of F above the horizon, with
    K = P(x)/Δ(x); or a NoOrbit."""
    found = horizons(a, b)
    if found is None:
        return NoOrbit("no horizon")
    x_plus, e = found
    P, F = shifted_coefficients(a, b, Phi, Psi, x_plus, e)
    if not P[3] > 0:
        return NoOrbit("at horizon")     # 𝓔(x₊) = 0: P and Δ both vanish there, and V = P/Δ is 0/0
    # P₀ > 0, so F < 0 just above the horizon (V → +∞ there), and F → +∞: bisect on u in (0, hi).
    lo, hi = 0.0, 4.0
    while polyval(F, hi) <= 0:
        hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if polyval(F, mid) > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo <= 1e-15 * hi:
            break
    u = 0.5 * (lo + hi)
    delta = u * (u + e)
    return x_plus + u, polyval(P, u) / delta, delta, extra_roots(F, u)


def reduce(q, at_zero, at_one):
    """Q divided by u if it vanishes at u = 0, and by (1 − u) if it vanishes at u = 1.

    The division by (1 − u) works up from the constant term, (1 − u)(s₀ + s₁u + …) = c₀ + c₁u + …, so the
    constant term −Φ² of the quotient stays exact however small Φ is.
    """
    r = list(q)
    if at_zero:
        r = r[:-1]
    if at_one:
        s = [r[-1]]
        for coefficient in reversed(r[1:-1]):
            s.append(coefficient + s[-1])
        r = s[::-1]
    return r


def band_peak(a, b, Phi, Psi, K):
    """(y_m, R(y_m), R): where R, which has the sign of Θ, is largest on [0, 1]; None if there is nowhere to look.

    R is Q, divided by y if Φ = 0 and by 1 − y if Ψ = 0, so that it stays finite and nonzero on the axes the
    orbit can then reach. Otherwise R(0) = −Φ² < 0 and R(1) = −Ψ² < 0, and the largest value is inside.
    """
    r = reduce(q_coefficients(a, b, Phi, Psi, K), Phi == 0, Psi == 0)
    degree = len(r) - 1
    slope = [0.0] * (3 - degree) + [r[i] * (degree - i) for i in range(degree)]
    candidates = ([0.0] if Phi == 0 else []) + ([1.0] if Psi == 0 else [])
    candidates += [t for t in quadratic_roots(*slope) if 0 < t < 1]
    if not candidates:
        return None
    best, top = candidates[0], polyval(r, candidates[0])
    for t in candidates[1:]:
        value = polyval(r, t)
        if value > top:
            best, top = t, value
    return best, top, r


def lower_root(c, hi):
    """The root of c in (0, hi), where c(0) < 0 < c(hi), by bisection to full relative precision."""
    lo = 0.0
    for _ in range(1100):
        mid = 0.5 * (lo + hi)
        if polyval(c, mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-15 * hi:
            break
    return 0.5 * (lo + hi)


def band(a, b, Phi, Psi, K):
    """(y1, w2) = (sin²θ_min, cos²θ_max), the edges of the band where Θ ≥ 0; None if there is no band.

    The upper edge is found as w2 = cos²θ_max from the same cubic with (a, Φ) and (b, Ψ) swapped, the paper's
    symmetry (2.2), so both edges keep full relative precision however close they come to an axis.
    """
    peak = band_peak(a, b, Phi, Psi, K)
    if peak is None:
        return None
    top_y, top, r = peak
    tol = BAND_TOL * (1 + abs(K) + Phi * Phi + Psi * Psi)
    if top < -tol:
        return None
    if top <= tol:
        return top_y, 1 - top_y                      # a circular orbit at constant θ
    flipped = reduce(q_coefficients(b, a, Psi, Phi, K), Psi == 0, Phi == 0)
    y1 = 0.0 if polyval(r, 0.0) >= 0 else lower_root(r, top_y)
    w2 = 0.0 if polyval(flipped, 0.0) >= 0 else lower_root(flipped, 1 - top_y)
    return y1, w2


def orbit_exists(a, b, Phi, Psi):
    """Whether these momenta have a spherical orbit: the same decision as solve_orbit, without the band's edges."""
    Phi, Psi = snap(Phi), snap(Psi)
    found = radial(a, b, Phi, Psi)
    if isinstance(found, NoOrbit):
        return False
    peak = band_peak(a, b, Phi, Psi, found[1])
    return peak is not None and peak[1] >= -BAND_TOL * (1 + abs(found[1]) + Phi * Phi + Psi * Psi)


def solve_orbit(a, b, Phi, Psi):
    """The spherical photon orbit with momenta (Φ, Ψ) around the hole with spins (a, b), or a NoOrbit."""
    Phi, Psi = snap(Phi), snap(Psi)
    found = radial(a, b, Phi, Psi)
    if isinstance(found, NoOrbit):
        return found
    x, K, delta, extra = found
    edges = band(a, b, Phi, Psi, K)
    if edges is None:
        return NoOrbit("no band", x, K)
    if extra:
        warnings.warn(f"F has {extra} more roots above the horizon for a={a}, b={b}, Phi={Phi}, Psi={Psi}",
                      SecondRootWarning, stacklevel=2)
    y1, w2 = edges
    y2 = 1 - w2
    A2, B2 = a * a, b * b
    calE = 1 + a * Phi / (x + A2) + b * Psi / (x + B2)
    A_phi = a * (x + B2) * calE / delta + (A2 - B2) * Phi / (x + A2)
    A_psi = b * (x + A2) * calE / delta - (A2 - B2) * Psi / (x + B2)
    G1 = (B2 - A2) * (y1 + y1 + y2) - (B2 - 2 * A2 - K)
    G2 = (B2 - A2) * (y2 + y1 + y2) - (B2 - 2 * A2 - K)
    return Orbit(a, b, Phi, Psi, horizon(a, b), x, K, y1, y2, w2, delta, calE, A_phi, A_psi, G1, G2, extra)
