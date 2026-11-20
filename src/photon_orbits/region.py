"""The allowed region 𝒜(a, b): the momenta (Φ, Ψ) that have a spherical orbit.

𝒜 is a convex "curved diamond" around the origin, so it is star-shaped and each ray from the origin leaves it
once; its edge is found by bisection along rays. demo/src/physics/region.js is a line-by-line port.
"""

import math

from .orbit import orbit_exists, sqrt0

REACH = 4.0             # 𝒜 lies inside |Φ|, |Ψ| < 3 for spins up to 1; the search starts from here
EDGE_STEPS = 48         # bisection steps along a ray: 4/2⁴⁸ ≈ 1.4e-14


def vertices(a, b):
    """The corners of 𝒜 on the axes, as (Φ, Ψ, x), for a, b ≥ 0: the equatorial circular orbits (Ψ = 0, eqs
    3.9–3.12), counter- then co-rotating, and the polar ones (Φ = 0, eqs 3.13–3.15)."""
    sp, sm = sqrt0((1 + a) * (1 + a) - b * b), sqrt0((1 - a) * (1 - a) - b * b)
    tp, tm = sqrt0((1 + b) * (1 + b) - a * a), sqrt0((1 - b) * (1 - b) - a * a)
    return [(1 + sp, 0.0, 1 + a - b * b + sp),
            (-1 - sm, 0.0, 1 - a - b * b + sm),
            (0.0, 1 + tp, 1 + b - a * a + tp),
            (0.0, -1 - tm, 1 - b - a * a + tm)]


def edge(a, b, d_phi, d_psi):
    """How far 𝒜 reaches along the unit direction (d_phi, d_psi): the largest s with an orbit at s·(d_phi, d_psi)."""
    lo, hi = 0.0, REACH
    while orbit_exists(a, b, hi * d_phi, hi * d_psi):
        hi *= 2
    for _ in range(EDGE_STEPS):
        mid = 0.5 * (lo + hi)
        if orbit_exists(a, b, mid * d_phi, mid * d_psi):
            lo = mid
        else:
            hi = mid
    return lo


def directions(rays):
    """rays unit vectors (cos α, sin α), α = 2πk/rays, starting along +Φ and turning towards +Ψ."""
    return [(math.cos(2 * math.pi * k / rays), math.sin(2 * math.pi * k / rays)) for k in range(rays)]


def boundary(a, b, rays=96):
    """The edge of 𝒜 as rays points (Φ, Ψ)."""
    points = []
    for d_phi, d_psi in directions(rays):
        s = edge(a, b, d_phi, d_psi)
        points.append((s * d_phi, s * d_psi))
    return points
