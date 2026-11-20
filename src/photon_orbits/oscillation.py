"""One latitude oscillation of an orbit, tabulated in the angle variable χ, and the path drawn from it.

With Mino time σ (dσ = dλ/ρ²) and y = sin²θ = ȳ − Δy cos χ, every rate depends on χ alone, so one oscillation
(χ from 0 to 2π) is a quadrature and longer paths repeat it. χ = 0 is θ_min and χ = π is θ_max.
docs/equations.md gives the derivations; demo/src/physics/oscillation.js is a line-by-line port.
"""

import math
from dataclasses import dataclass

TWO_PI = 2 * math.pi
AXIS_NEAR = 1e-4        # a band within this of an axis, in sin²θ or cos²θ, gets 512 steps instead of 256


@dataclass(frozen=True)
class Table:
    """σ, λ, φ and ψ at the nodes χ_i = 2πi/n, i = 0 … n, starting from 0 at θ_min."""

    n: int
    sigma: tuple
    lam: tuple
    phi: tuple
    psi: tuple

    @property
    def step(self):
        return TWO_PI / self.n


def steps_for(orbit):
    return 512 if min(orbit.y1, orbit.w2) < AXIS_NEAR else 256


def axis_integral(at0, atpi, i, n):
    """∫₀^χ dχ′ / (at0 + (atpi − at0) sin²(χ′/2)) at χ = 2πi/n (n even, 0 ≤ i ≤ n), for at0, atpi > 0.

    It is (2/√(at0·atpi)) arctan(k tan(χ/2)) with k = √(atpi/at0), on the branch continuous in χ. When the
    band nearly touches an axis the integrand is a spike of width ~1e-12 at χ = 0 or π, so the angle comes
    from the index i, and the symmetries about χ = π/2 and π make the values at 0, π and 2π exact. i need not
    be an integer: `between` uses it between the nodes.
    """
    scale = 2 / math.sqrt(at0 * atpi)
    k = math.sqrt(atpi / at0)
    j = i if 2 * i <= n else n - i                  # fold (π, 2π] onto [0, π)
    if 4 * j <= n:
        angle = math.atan(k * math.tan(math.pi * j / n))
    else:
        angle = 0.5 * math.pi - math.atan(math.tan(math.pi * (n - 2 * j) / (2 * n)) / k)
    return scale * angle if j == i else scale * (math.pi - angle)


def rates(orbit, chi):
    """d(σ, λ, φ, ψ)/dχ at χ, less the axis terms Φ/(2y√G₁) and Ψ/(2(1 − y)√G₂), which tabulate integrates exactly."""
    A2, B2 = orbit.a * orbit.a, orbit.b * orbit.b
    s, c = math.sin(0.5 * chi), math.cos(0.5 * chi)
    width = orbit.y2 - orbit.y1
    rise = width * s * s            # y − y1
    fall = width * c * c            # y2 − y
    y = orbit.y1 + rise
    w = orbit.w2 + fall             # 1 − y, without the rounding of 1 − y near the equator
    root = math.sqrt((B2 - A2) * (y + orbit.y1 + orbit.y2) - (B2 - 2 * A2 - orbit.K))
    dsigma = 0.5 / root
    dlam = (orbit.x + A2 + (B2 - A2) * y) * dsigma
    dphi = -orbit.A_phi * dsigma
    if orbit.Phi != 0:
        root1 = math.sqrt(orbit.G1)
        dphi -= orbit.Phi * (B2 - A2) * rise / (2 * y * root * root1 * (root + root1))
    dpsi = -orbit.A_psi * dsigma
    if orbit.Psi != 0:
        root2 = math.sqrt(orbit.G2)
        dpsi += orbit.Psi * (B2 - A2) * fall / (2 * w * root * root2 * (root + root2))
    return dsigma, dlam, dphi, dpsi


def tabulate(orbit, n=None):
    """One oscillation by Simpson's rule on n steps (256, or 512 near an axis), with the axis terms in closed form."""
    n = n or steps_for(orbit)
    if n % 2:
        raise ValueError("the number of steps must be even")
    h = TWO_PI / n
    f = [rates(orbit, 0.5 * h * j) for j in range(2 * n + 1)]
    columns = [[0.0] for _ in range(4)]
    for i in range(n):
        for k, column in enumerate(columns):
            column.append(column[-1] + h / 6 * (f[2 * i][k] + 4 * f[2 * i + 1][k] + f[2 * i + 2][k]))
    sigma, lam, phi, psi = columns
    if orbit.Phi != 0:
        scale = orbit.Phi / (2 * math.sqrt(orbit.G1))
        for i in range(n + 1):
            phi[i] += scale * axis_integral(orbit.y1, orbit.y2, i, n)
    if orbit.Psi != 0:
        scale = orbit.Psi / (2 * math.sqrt(orbit.G2))
        for i in range(n + 1):
            psi[i] += scale * axis_integral(1 - orbit.y1, orbit.w2, i, n)
    return Table(n, tuple(sigma), tuple(lam), tuple(phi), tuple(psi))


def node(table, i):
    """(χ, σ, λ, φ, ψ) at node i, for any integer i: oscillation m = ⌊i/n⌋ is the table shifted by m increments."""
    m, j = divmod(i, table.n)
    return (j * table.step + m * TWO_PI,
            table.sigma[j] + m * table.sigma[-1],
            table.lam[j] + m * table.lam[-1],
            table.phi[j] + m * table.phi[-1],
            table.psi[j] + m * table.psi[-1])


def between(orbit, table, i, t):
    """(χ, σ, λ, φ, ψ) at χ = (i + t)·2π/n, a fraction 0 ≤ t ≤ 1 of the way from node i to node i + 1.

    One Simpson step of width t·h from node i, with the axis terms in closed form: it is node i at t = 0 and
    meets node i + 1 at t = 1. Near an axis φ or ψ swings by about π between two nodes, and this follows the
    swing, so the drawing can subdivide it.
    """
    n = table.n
    m, j = divmod(i, n)
    h = table.step
    chi = (j + t) * h
    f0, f1, f2 = rates(orbit, j * h), rates(orbit, (j + 0.5 * t) * h), rates(orbit, chi)
    sigma, lam, phi, psi = (column[j] + t * h / 6 * (f0[k] + 4 * f1[k] + f2[k])
                            for k, column in enumerate((table.sigma, table.lam, table.phi, table.psi)))
    if orbit.Phi != 0:
        scale = orbit.Phi / (2 * math.sqrt(orbit.G1))
        phi += scale * (axis_integral(orbit.y1, orbit.y2, j + t, n) - axis_integral(orbit.y1, orbit.y2, j, n))
    if orbit.Psi != 0:
        scale = orbit.Psi / (2 * math.sqrt(orbit.G2))
        psi += scale * (axis_integral(1 - orbit.y1, orbit.w2, j + t, n) - axis_integral(1 - orbit.y1, orbit.w2, j, n))
    return (chi + m * TWO_PI,
            sigma + m * table.sigma[-1],
            lam + m * table.lam[-1],
            phi + m * table.phi[-1],
            psi + m * table.psi[-1])


def latitude(orbit, chi):
    """(sin θ, cos θ) at χ, signed where the orbit passes through an axis, so that it carries straight on."""
    s, c = math.sin(0.5 * chi), math.cos(0.5 * chi)
    width = orbit.y2 - orbit.y1
    through_pole = orbit.Phi == 0 and orbit.y1 == 0
    through_equator = orbit.Psi == 0 and orbit.w2 == 0
    sin_theta = math.sqrt(orbit.y2) * s if through_pole else math.sqrt(orbit.y1 + width * s * s)
    cos_theta = math.sqrt(1 - orbit.y1) * c if through_equator else math.sqrt(orbit.w2 + width * c * c)
    return sin_theta, cos_theta


def embed(orbit, chi, phi):
    """(X, Y, Z) in the paper's picture: the a-plane flat, and the b-plane's radius as height. Colour it by ψ."""
    sin_theta, cos_theta = latitude(orbit, chi)
    across = math.sqrt(orbit.x + orbit.a * orbit.a) * sin_theta
    return across * math.cos(phi), across * math.sin(phi), math.sqrt(orbit.x + orbit.b * orbit.b) * cos_theta


def frequencies(orbit, table):
    """(ν_φ, ν_ψ): the turns of φ and ψ per oscillation, counting a pass through an axis as half a turn."""
    dphi = table.phi[-1] + (math.pi if orbit.Phi == 0 and orbit.y1 == 0 else 0.0)
    dpsi = table.psi[-1] + (math.pi if orbit.Psi == 0 and orbit.w2 == 0 else 0.0)
    return dphi / TWO_PI, dpsi / TWO_PI
