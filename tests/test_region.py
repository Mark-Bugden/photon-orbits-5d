"""The allowed region of momenta: its corners, its shape, its edge, and extremal black holes."""

import math
import warnings

import numpy as np
import pytest

from photon_orbits import SecondRootWarning, boundary, directions, edge, orbit_exists, solve_orbit, tabulate, vertices
from photon_orbits.orbit import radial
from reference_values import STRESS_HOLES, VERTICES

HOLES = [(0, 0), (0.5, 0), (0.6, 0.3), (0.4, 0.4), (0.9, 0.09), (0.3, 0.6), (0.99, 0)]


@pytest.mark.parametrize("row", VERTICES, ids=lambda row: f"a={row[0]} b={row[1]}")
def test_the_corners(row):
    a, b = row[:2]
    expected = [(row[2], 0, row[3]), (row[4], 0, row[5]), (0, row[6], row[7]), (0, row[8], row[9])]
    for (Phi, Psi, x), (want_Phi, want_Psi, want_x) in zip(vertices(a, b), expected):
        assert (Phi, Psi, x) == pytest.approx((want_Phi, want_Psi, want_x), abs=1.5e-6)
        # The ray search finds the same corner, and the orbit just inside it is the circular orbit there.
        reach = math.hypot(Phi, Psi)
        assert edge(a, b, Phi / reach, Psi / reach) == pytest.approx(reach, abs=1e-9)
        inside = solve_orbit(a, b, 0.9999999 * Phi, 0.9999999 * Psi)
        assert inside.x == pytest.approx(x, abs=1e-5)


@pytest.mark.parametrize("a, b", HOLES)
def test_the_region_is_star_shaped_and_convex(a, b):
    rays = directions(24)
    reach = [edge(a, b, *d) for d in rays]
    for (d_phi, d_psi), s in zip(rays, reach):
        for f in np.linspace(0, 0.999, 25):
            assert orbit_exists(a, b, f * s * d_phi, f * s * d_psi)
        for f in np.linspace(1.001, 2, 10):
            assert not orbit_exists(a, b, f * s * d_phi, f * s * d_psi)
    points = boundary(a, b, 48)
    rng = np.random.default_rng(5)
    for _ in range(100):
        i, j = rng.integers(0, len(points), 2)
        assert orbit_exists(a, b, *(0.9999 * 0.5 * (np.array(points[i]) + points[j])))


@pytest.mark.parametrize("a, b", [(0.6, 0.3), (0.9, 0.09), (0.4, 0.4), (0.5, 0)])
def test_the_edge_is_circular_orbits(a, b):
    """Off extremality, the edge is made of circular orbits of constant θ, as in Fig 4.12."""
    for d_phi, d_psi in directions(16):
        s = edge(a, b, d_phi, d_psi)
        orbit = solve_orbit(a, b, (1 - 1e-9) * s * d_phi, (1 - 1e-9) * s * d_psi)
        assert math.asin(math.sqrt(1 - orbit.w2)) - math.asin(math.sqrt(orbit.y1)) < 2e-3


@pytest.mark.parametrize("a, b", [(0.9, 0.1), (0.5, 0.5), (0.7, 0.3), (0.99, 0.01), (0.25, 0.75)])
def test_extremal_holes_meet_the_horizon_along_phi_plus_psi_minus_1(a, b):
    """At a + b = 1 the co-rotating edge (Φ, Ψ < 0) is the line Φ + Ψ = −1, where 𝓔(x₊) = 1 + Φ + Ψ = 0.

    Approaching it, the orbit closes on the horizon, x − x₊ ∝ √(1 + Φ + Ψ), with a wide θ band; past it there
    is no θ-motion. Neither the paper nor the review has this: the review's solver stopped about 3e-4 short.
    """
    for k in range(1, 24):
        angle = math.pi + 0.5 * math.pi * k / 24
        d_phi, d_psi = math.cos(angle), math.sin(angle)
        s = edge(a, b, d_phi, d_psi)
        assert s * (d_phi + d_psi) == pytest.approx(-1, abs=1e-8)      # floating-point spins are extremal to ~1e-8
        gaps = []
        for gap in (1e-2, 1e-4, 1e-6):
            orbit = solve_orbit(a, b, (1 - gap) * s * d_phi, (1 - gap) * s * d_psi)
            assert orbit.ok
            assert math.asin(math.sqrt(1 - orbit.w2)) - math.asin(math.sqrt(orbit.y1)) > 0.01
            gaps.append(orbit.x - orbit.x_plus)
        assert gaps[1] / gaps[0] == pytest.approx(0.1, rel=0.1) and gaps[2] / gaps[1] == pytest.approx(0.1, rel=0.1)
        assert not orbit_exists(a, b, (1 + 1e-6) * s * d_phi, (1 + 1e-6) * s * d_psi)


@pytest.mark.parametrize("a, b", [(0.9, 0.1), (0.5, 0.5), (0.7, 0.3)])
def test_near_the_extremal_edge_the_orbit_follows_the_depressed_cubic(a, b):
    """At extremality F = u(u³ − P₁u − 2P₀) with u = x − x₊, P₀ = ab𝓔₊² and P₁ ≈ −2(bΦ + aΨ)𝓔₊, 𝓔₊ = 1 + Φ + Ψ.
    Inside the edge u ≈ √(2|bΦ + aΨ|𝓔₊) and K → 1 + ab; outside u ≈ ab|𝓔₊|/|bΦ + aΨ| and
    K → 1 + ab − (bΦ + aΨ)²/(ab), with no θ-motion (docs/equations.md, "Extremal holes")."""
    for t in (0.0, 0.3, 0.5, 1.0):
        for eps in (1e-6, -1e-6):
            Phi, Psi = -(1 - eps) * t, -(1 - eps) * (1 - t)
            d = abs(b * Phi + a * Psi)
            x, K, _, _ = radial(a, b, Phi, Psi)
            if eps > 0:
                assert x - a * b == pytest.approx(math.sqrt(2 * d * eps), rel=2e-3)
                assert K == pytest.approx(1 + a * b, abs=5e-3)
                assert orbit_exists(a, b, Phi, Psi)
            else:
                assert x - a * b == pytest.approx(a * b * -eps / d, rel=1e-4)
                assert K == pytest.approx(1 + a * b - d * d / (a * b), abs=1e-4)
                assert not orbit_exists(a, b, Phi, Psi)


def test_on_the_line_itself_the_orbit_would_sit_on_the_horizon():
    """𝓔(x₊) = 0 exactly: P and Δ both vanish at the horizon, V = P/Δ is 0/0 there, and there is no orbit."""
    found = solve_orbit(0.5, 0.5, -0.5, -0.5)
    assert not found.ok and found.reason == "at horizon"
    assert not orbit_exists(0.5, 0.5, -0.25, -0.75)


@pytest.mark.parametrize("a, b", STRESS_HOLES)
def test_the_plans_extremal_stress_test(a, b):
    """96 directions × 8 fractions of the way to the edge: every orbit exists, its table is finite, and none
    warns of a second root. Up to 0.999 of the way, the orbit stays at least 0.001 outside the horizon."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", SecondRootWarning)
        for d_phi, d_psi in directions(96):
            s = edge(a, b, d_phi, d_psi)
            for f in (0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999):
                orbit = solve_orbit(a, b, f * s * d_phi, f * s * d_psi)
                assert orbit.ok and orbit.x - orbit.x_plus > 1e-3
                assert orbit.G1 > 0 and orbit.G2 > 0
                table = tabulate(orbit)
                assert all(math.isfinite(v) for v in table.sigma + table.lam + table.phi + table.psi)
