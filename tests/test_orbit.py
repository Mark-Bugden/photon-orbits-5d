"""The horizon, the radius, K and the latitude band, against the paper review's values and independent checks."""

import math
import warnings

import numpy as np
import pytest

from photon_orbits import NoOrbit, SecondRootWarning, horizon, orbit_exists, solve_orbit
from photon_orbits.orbit import horizons, p_coefficients, polyval, q_coefficients, radial, shifted_coefficients
from reference_values import FIG_4_12, TABLE


@pytest.mark.parametrize("row", TABLE, ids=lambda row: f"a={row[0]} b={row[1]} Phi={row[2]} Psi={row[3]}")
def test_the_review_table(row):
    a, b, Phi, Psi, x_plus, x, r, K, theta_min, theta_max = row[:10]
    orbit = solve_orbit(a, b, Phi, Psi)
    assert orbit.ok
    assert orbit.x_plus == pytest.approx(x_plus, abs=1.5e-6)
    assert orbit.x == pytest.approx(x, abs=1.5e-6)
    assert math.sqrt(orbit.x) == pytest.approx(r, abs=1.5e-6)
    assert orbit.K == pytest.approx(K, abs=1.5e-6)
    assert math.degrees(math.asin(math.sqrt(orbit.y1))) == pytest.approx(theta_min, abs=1.5e-4)
    assert math.degrees(math.acos(math.sqrt(orbit.w2))) == pytest.approx(theta_max, abs=1.5e-4)
    assert orbit.extra_roots == 0
    assert orbit.G1 > 0 and orbit.G2 > 0


@pytest.mark.parametrize("Phi, Psi, x, K", FIG_4_12)
def test_the_circular_orbits_at_45_degrees(Phi, Psi, x, K):
    orbit = solve_orbit(0.5, 0.25, Phi, Psi)
    assert orbit.x == pytest.approx(x, abs=1.5e-6)
    assert orbit.K == pytest.approx(K, abs=1.5e-6)
    # The momenta are rounded to 6 digits, so the band is a hair wide around 45°, by up to 0.03°.
    assert math.degrees(math.asin(math.sqrt(orbit.y1))) == pytest.approx(45, abs=0.05)
    assert math.degrees(math.acos(math.sqrt(orbit.w2))) == pytest.approx(45, abs=0.05)


def test_the_horizon():
    assert horizon(0, 0) == 1
    assert horizon(0.5, 0) == pytest.approx(0.75)
    for a, b in [(0.3, 0.2), (0.9, 0.099), (0.5, 0.5), (0.9, 0.1), (0.99, 0.01)]:
        x_plus, e = horizons(a, b)
        assert (x_plus + a * a) * (x_plus + b * b) - x_plus == pytest.approx(0, abs=1e-15)
        assert e == pytest.approx(math.sqrt(max((1 - (a + b) ** 2) * (1 - (a - b) ** 2), 0)), abs=1e-8)
    assert horizons(0.9, 0.1)[1] == 0 and horizon(0.9, 0.1) == pytest.approx(0.09, abs=1e-16)   # extremal: x₊ = ab
    for a, b in [(0.6, 0.6), (1.0, 0.0), (0.0, 1.0), (0.8, 0.3)]:
        assert horizon(a, b) is None
        assert solve_orbit(a, b, 0.1, 0.1) == NoOrbit("no horizon")


def test_the_shifted_coefficients_are_p_and_f():
    """P and F in u = x − x₊ agree with P(x) from the paper's 𝒳 and with F = P′Δ − PΔ′ formed directly."""
    rng = np.random.default_rng(3)
    for _ in range(300):
        a = rng.uniform(0, 0.99)
        b = rng.uniform(0, 1 - a)
        Phi, Psi = rng.uniform(-3, 3, 2)
        x_plus, e = horizons(a, b)
        P_u, F_u = shifted_coefficients(a, b, Phi, Psi, x_plus, e)
        P = np.array(p_coefficients(a, b, Phi, Psi))
        D = np.array([1, a * a + b * b - 1, a * a * b * b])
        F = np.polysub(np.polymul(np.polyder(P), D), np.polymul(P, np.polyder(D)))
        for u in (0.05, 0.3, 1.0, 2.5):
            assert polyval(P_u, u) == pytest.approx(np.polyval(P, x_plus + u), rel=1e-11, abs=1e-12)
            assert polyval(F_u, u) == pytest.approx(np.polyval(F, x_plus + u), rel=1e-10, abs=1e-12)


def test_one_orbit_per_pair_of_momenta():
    """F has exactly one root above the horizon (numpy's roots of the quartic in x), and it is the one found."""
    rng = np.random.default_rng(7)
    for _ in range(3000):
        a = rng.uniform(0, 0.99)
        b = rng.uniform(0, 1 - a) if rng.random() < 0.5 else (1 - a) * (1 - 10 ** rng.uniform(-4, -1))
        Phi, Psi = rng.uniform(-3, 3, 2)
        x_plus = horizon(a, b)
        P = np.array(p_coefficients(a, b, Phi, Psi))
        D = np.array([1, a * a + b * b - 1, a * a * b * b])
        roots = np.roots(np.polysub(np.polymul(np.polyder(P), D), np.polymul(P, np.polyder(D))))
        above = [r.real for r in roots if abs(r.imag) < 1e-9 and r.real > x_plus + 1e-6]
        x, K, delta, extra = radial(a, b, Phi, Psi)
        assert len(above) == 1 and extra == 0, (a, b, Phi, Psi, above)
        assert x == pytest.approx(above[0], rel=1e-8)
        assert K == pytest.approx(np.polyval(P, x) / np.polyval(D, x), rel=1e-9)
        assert delta == pytest.approx(np.polyval(D, x), rel=1e-9)


@pytest.mark.parametrize("a, b, Phi, Psi, x, K", [
    (0.4, 0.4, 1.5, -1.5, 1.659921, 3.664599),     # the review's counterexample to §3.1.1 of the paper
    (0.55, 0.2, -1.3, 0.6, None, None),
    (0, 0, 1.5, 0.6, 2, 4),                         # Tangherlini: |Φ| + |Ψ| > 2
])
def test_momenta_outside_the_region_have_no_orbit(a, b, Phi, Psi, x, K):
    found = solve_orbit(a, b, Phi, Psi)
    assert not found.ok and found.reason == "no band"
    assert not orbit_exists(a, b, Phi, Psi)
    if x is not None:
        assert found.x == pytest.approx(x, abs=1.5e-6) and found.K == pytest.approx(K, abs=1.5e-6)


def test_equal_spins_depend_only_on_the_sum_of_momenta():
    for Phi in (0.35, -0.3, 0.6, 1.0, -0.65):
        orbit = solve_orbit(0.4, 0.4, Phi, 0.7 - Phi)
        assert orbit.x == pytest.approx(1.967354452441, abs=1e-11)
        assert orbit.K == pytest.approx(4.232684973306, abs=1e-11)


def test_tangherlini_orbits_all_have_x_2_and_k_4():
    rng = np.random.default_rng(11)
    for _ in range(200):
        Phi, Psi = rng.uniform(-2, 2, 2)
        if abs(Phi) + abs(Psi) < 1.999:
            orbit = solve_orbit(0, 0, Phi, Psi)
            assert orbit.x == pytest.approx(2, abs=1e-13) and orbit.K == pytest.approx(4, abs=1e-12)


def test_tiny_momenta_count_as_zero():
    """Below 1e-12 a momentum is exactly 0. Without that, Φ = 1e-16 (from a ray straight up the Ψ axis) gave
    the research script a band edge of y₁ ≈ 0 while it still divided by y."""
    for a, b in [(0.6, 0.3), (0.9, 0.1), (0.8, 0), (0, 0)]:
        assert solve_orbit(a, b, 1e-16, 0.5) == solve_orbit(a, b, 0.0, 0.5)
        assert solve_orbit(a, b, 0.5, -1e-16) == solve_orbit(a, b, 0.5, 0.0)
        assert solve_orbit(a, b, math.cos(math.pi / 2), 1.0) == solve_orbit(a, b, 0.0, 1.0)


def test_small_but_real_momenta_keep_their_band_edges_exact():
    """sin²θ_min ≈ Φ²/(a² + K) and cos²θ_max ≈ Ψ²/(b² + K) are far below the spacing of floats near 1, but both
    edges are roots of Q to full relative precision: the upper one is found in w = cos²θ."""
    for a, b in [(0.6, 0.3), (0.9, 0.09)]:
        for value in (1e-12, 1e-10, 1e-8, 1e-6):
            orbit = solve_orbit(a, b, value, -value)
            assert 0 < orbit.y1 < 1e-9 and 0 < orbit.w2 < 1e-9
            assert abs(polyval(q_coefficients(a, b, value, -value, orbit.K), orbit.y1)) <= 1e-12 * value * value
            assert abs(polyval(q_coefficients(b, a, -value, value, orbit.K), orbit.w2)) <= 1e-12 * value * value


def test_no_second_root_warnings_on_extremal_holes():
    with warnings.catch_warnings():
        warnings.simplefilter("error", SecondRootWarning)
        for a, b in [(0.9, 0.1), (0.5, 0.5), (0.7, 0.3), (0.99, 0.01)]:
            for k in range(200):
                angle = math.pi + 0.5 * math.pi * k / 199
                for s in (0.7, 0.9, 0.999999, 1.0, 1.000001, 1.1):
                    solve_orbit(a, b, s * math.cos(angle) / (-math.cos(angle) - math.sin(angle)),
                                s * math.sin(angle) / (-math.cos(angle) - math.sin(angle)))
