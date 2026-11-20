"""One oscillation in χ, its increments and frequencies, periodic extension, and the drawn path."""

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ellipk

from photon_orbits import between, embed, frequencies, latitude, node, solve_orbit, tabulate
from photon_orbits.oscillation import steps_for
from reference_values import TABLE

ids = lambda row: f"a={row[0]} b={row[1]} Phi={row[2]} Psi={row[3]}"  # noqa: E731


@pytest.mark.parametrize("row", TABLE, ids=ids)
def test_the_review_table(row):
    a, b, Phi, Psi = row[:4]
    d_sigma, d_lam, d_phi, d_psi, nu_phi, nu_psi = row[10:]
    orbit = solve_orbit(a, b, Phi, Psi)
    table = tabulate(orbit)
    assert table.sigma[-1] == pytest.approx(d_sigma, abs=1.5e-6)
    assert table.lam[-1] == pytest.approx(d_lam, abs=1.5e-6)
    assert table.phi[-1] == pytest.approx(d_phi, abs=1.5e-6)
    assert table.psi[-1] == pytest.approx(d_psi, abs=1.5e-6)
    assert frequencies(orbit, table) == pytest.approx((nu_phi, nu_psi), abs=1.5e-6)


@pytest.mark.parametrize("row", [row for row in TABLE if row[8] != row[9]], ids=ids)
def test_the_period_in_mino_time_is_a_complete_elliptic_integral(row):
    """Δσ = 2K(m)/√G_max with m = 1 − G_min/G_max (paper review §B3)."""
    orbit = solve_orbit(*row[:4])
    big, small = max(orbit.G1, orbit.G2), min(orbit.G1, orbit.G2)
    assert tabulate(orbit).sigma[-1] == pytest.approx(2 * ellipk(1 - small / big) / math.sqrt(big), abs=1e-11)


@pytest.mark.parametrize("a, Phi, Psi", [(0.4, 1.2, -0.5), (0.3, -0.7, -0.9), (0.45, 0.6, 1.1), (0.2, 0.0, 0.8)])
def test_equal_spins_have_closed_forms(a, Phi, Psi):
    """For a = b, G = a² + K is constant, Δσ = π/√(a² + K) and Δφ = π[sgn Φ − A_φ/√(a² + K)], with A_φ = A_ψ."""
    orbit = solve_orbit(a, a, Phi, Psi)
    table = tabulate(orbit)
    root = math.sqrt(a * a + orbit.K)
    assert orbit.A_phi == pytest.approx(orbit.A_psi, abs=1e-14)
    assert table.sigma[-1] == pytest.approx(math.pi / root, abs=1e-12)
    sign = lambda v: 0.0 if v == 0 else math.copysign(1, v)  # noqa: E731   (Φ = 0: no axis term, as in Δφ_eff − π)
    assert table.phi[-1] == pytest.approx(math.pi * (sign(Phi) - orbit.A_phi / root), abs=1e-11)
    assert table.psi[-1] == pytest.approx(math.pi * (sign(Psi) - orbit.A_psi / root), abs=1e-11)


def test_tangherlini_orbits_are_great_circles():
    for Phi, Psi in [(0.7, 0.9), (-1.2, 0.3), (0.0, 1.5), (1.0, 0.0), (0.0, 0.0), (-0.4, -1.5)]:
        orbit = solve_orbit(0, 0, Phi, Psi)
        table = tabulate(orbit)
        assert table.sigma[-1] == pytest.approx(math.pi / 2, abs=1e-13)
        assert table.lam[-1] == pytest.approx(math.pi, abs=1e-12)
        nu = frequencies(orbit, table)
        assert [v % 1 for v in nu] == pytest.approx([0.5, 0.5], abs=1e-12)


def test_the_closed_orbit_closes_after_15_oscillations():
    """The review's resonance (§B6): ν = (1/3, −3/5) at a = 0.6, b = 0.3."""
    orbit = solve_orbit(0.6, 0.3, 0.1051348934, -0.1354415960)
    nu_phi, nu_psi = frequencies(orbit, tabulate(orbit))
    assert 15 * nu_phi == pytest.approx(5, abs=1e-8)
    assert 15 * nu_psi == pytest.approx(-9, abs=1e-8)


@pytest.mark.parametrize("a, b, Phi, Psi", [
    *[row[:4] for row in TABLE],
    (0.6, 0.3, 1e-3, 0.5), (0.6, 0.3, -0.8, -1e-3),        # the review's hard cases: a band edge near an axis
    (0.6, 0.3, 1e-8, 0.5), (0.9, 0.09, 0.5, 1e-9), (0.9, 0.1, -0.49, -0.495),
])
def test_the_table_has_converged(a, b, Phi, Psi):
    """The default 256 or 512 steps agree with 4096 steps to 1e-7, the accuracy the plan asks for near the axes."""
    orbit = solve_orbit(a, b, Phi, Psi)
    coarse, fine = tabulate(orbit), tabulate(orbit, 4096)
    step = fine.n // coarse.n
    for name in ("sigma", "lam", "phi", "psi"):
        assert getattr(coarse, name) == pytest.approx(getattr(fine, name)[::step], abs=1e-7), name


def test_512_steps_near_an_axis():
    assert steps_for(solve_orbit(0.6, 0.3, 1.0, -1.0)) == 256
    assert steps_for(solve_orbit(0.9, 0.09, 0.1, 0.01)) == 512        # θ_max = 89.69°: cos²θ_max = 3e-5
    assert steps_for(solve_orbit(0.8, 0.0, 0.0, 0.5)) == 512          # through the pole


@pytest.mark.parametrize("a, b, Phi, Psi", [(0.6, 0.3, 0.5, 0.0), (0.6, 0.3, 0.0, 0.5), (0.9, 0.09, 0.0, -0.6),
                                            (0.6, 0.3, 0.0, -0.5), (0.4, 0.4, 0.8, 0.0)])
def test_small_momenta_change_the_path_continuously(a, b, Phi, Psi):
    """ν, counting an axis passage as half a turn, is continuous mod 1 as Φ (or Ψ) passes through 0."""
    at_zero = frequencies(*(lambda o: (o, tabulate(o)))(solve_orbit(a, b, Phi, Psi)))
    move = (lambda v: (v, Psi)) if Phi == 0 else (lambda v: (Phi, v))
    for value in (1e-12, -1e-12, 1e-10, 1e-8, -1e-8, 1e-6, -1e-4):
        orbit = solve_orbit(a, b, *move(value))
        table = tabulate(orbit)
        assert all(math.isfinite(v) for v in table.phi + table.psi + table.lam + table.sigma)
        for got, want in zip(frequencies(orbit, table), at_zero):
            assert abs(math.remainder(got - want, 1)) <= abs(value) + 1e-9, value


def test_periodic_extension():
    orbit = solve_orbit(0.6, 0.3, 1.0, -1.0)
    table = tabulate(orbit)
    n = table.n
    for i in (0, 7, n // 2, n - 1):
        base = node(table, i)
        for m in (-2, 1, 5):
            shifted = node(table, i + m * n)
            assert shifted[0] == pytest.approx(base[0] + 2 * math.pi * m, abs=1e-12)
            for k, increment in enumerate((table.sigma[-1], table.lam[-1], table.phi[-1], table.psi[-1]), start=1):
                assert shifted[k] == pytest.approx(base[k] + m * increment, abs=1e-12)
    assert node(table, n)[1:] == pytest.approx(node(table, 0)[1:] + np.array(
        [table.sigma[-1], table.lam[-1], table.phi[-1], table.psi[-1]]), abs=1e-15)


BETWEEN = [
    (0.6, 0.3, 1.0, -1.0), (0.9, 0.09, 0.1, 0.01), (0.8, 0.0, 0.0, 0.5), (0.6, 0.3, 1.0, 0.0), (0.4, 0.4, 1.2, -0.5),
    (0.6, 0.3, 1e-3, 0.5), (0.6, 0.3, -0.8, -1e-3), (0.6, 0.3, 1e-8, 0.5), (0.9, 0.09, 0.5, 1e-9),   # near an axis
    (0.5, 0.5, -0.3, -0.6999), (0.7, 0.3, 0.0, -0.99999),                     # near the extremal edge: φ winds fast
]


@pytest.mark.parametrize("a, b, Phi, Psi", BETWEEN)
def test_between_the_nodes(a, b, Phi, Psi):
    """Between two nodes, the state agrees with a table 16 times finer, including next to χ = 0, π and 2π."""
    orbit = solve_orbit(a, b, Phi, Psi)
    table = tabulate(orbit)
    n = table.n
    fine = tabulate(orbit, 16 * n)
    for i in sorted({*range(0, n, 7), 1, n // 2 - 1, n // 2, n - 1}):
        assert between(orbit, table, i, 0.0) == node(table, i)
        assert between(orbit, table, i, 1.0)[1:] == pytest.approx(node(table, i + 1)[1:], abs=1e-12)
        for k in range(1, 16):
            got, want = between(orbit, table, i, k / 16), node(fine, 16 * i + k)
            assert got[0] == pytest.approx(want[0], abs=1e-12)
            assert got[1:] == pytest.approx(want[1:], abs=1e-7), (i, k)


def test_between_repeats_with_the_table():
    orbit = solve_orbit(0.9, 0.09, 0.1, 0.01)
    table = tabulate(orbit)
    n = table.n
    increments = np.array([2 * math.pi, table.sigma[-1], table.lam[-1], table.phi[-1], table.psi[-1]])
    for i, t in [(0, 0.3), (n // 2 - 1, 0.9), (n - 1, 0.5)]:
        base = np.array(between(orbit, table, i, t))
        for m in (-3, 1, 40):
            assert between(orbit, table, i + m * n, t) == pytest.approx(base + m * increments, abs=1e-11)


@pytest.mark.parametrize("a, b, Phi, Psi, angle, chi0", [
    (0.6, 0.3, 1e-3, 0.5, "phi", 0.0),          # φ swings by about π within ~1e-3 of χ = 0, near the pole
    (0.6, 0.3, -0.8, -1e-3, "psi", math.pi),    # ψ swings by about π within ~1e-3 of χ = π, near the equator
    (0.9, 0.09, 0.1, 0.01, "psi", math.pi),     # Fig 4.10
])
def test_between_follows_the_swing_near_an_axis(a, b, Phi, Psi, angle, chi0):
    """Across the swing, the change from the node agrees with adaptive quadrature of the full rate, axis term and all.

    The smooth remainder varies on the swing's own scale, so one Simpson step is good to about 1e-8 here, like the table.
    """
    orbit = solve_orbit(a, b, Phi, Psi)
    table = tabulate(orbit)
    A2, B2 = a * a, b * b
    width = orbit.y2 - orbit.y1

    def rate(chi):
        y = orbit.y1 + width * math.sin(0.5 * chi) ** 2
        G = (B2 - A2) * (y + orbit.y1 + orbit.y2) - (B2 - 2 * A2 - orbit.K)
        if angle == "phi":
            return (orbit.Phi / y - orbit.A_phi) / (2 * math.sqrt(G))
        return (orbit.Psi / (1 - y) - orbit.A_psi) / (2 * math.sqrt(G))

    column = 3 if angle == "phi" else 4
    centre = round(chi0 / table.step)
    swept = []
    for i in (centre - 1, centre):
        start = node(table, i)
        for t in (0.02, 0.1, 0.37, 0.5, 0.63, 0.9, 0.98):
            chi = (i + t) * table.step
            exact, _ = quad(rate, i * table.step, chi, limit=400, epsabs=1e-11, epsrel=1e-11)
            got = between(orbit, table, i, t)[column] - start[column]
            assert got == pytest.approx(exact, abs=1e-7), (i, t)      # the table's own accuracy
            swept.append(between(orbit, table, i, t)[column])
    assert abs(swept[-1] - swept[0]) > 1.5          # at least half of the swing happens within these two steps


@pytest.mark.parametrize("row", TABLE, ids=ids)
def test_the_path_lies_on_the_spheroid_and_is_continuous(row):
    a, b = row[:2]
    orbit = solve_orbit(*row[:4])
    table = tabulate(orbit)
    points = [embed(orbit, chi, phi) for chi, _, _, phi, _ in (node(table, i) for i in range(3 * table.n + 1))]
    for X, Y, Z in points:
        assert (X * X + Y * Y) / (orbit.x + a * a) + Z * Z / (orbit.x + b * b) == pytest.approx(1, abs=1e-12)
    jumps = [math.dist(p, q) for p, q in zip(points, points[1:])]
    assert max(jumps) < 0.2          # signed coordinates carry the path straight through the axes


def test_phi_zero_orbits_pass_through_the_pole():
    orbit = solve_orbit(0.8, 0.0, 0.0, 0.5)
    table = tabulate(orbit)
    X, Y, Z = embed(orbit, 0.0, 0.0)
    assert (X, Y) == (0, 0) and Z == pytest.approx(math.sqrt(orbit.x))
    # sin θ changes sign each oscillation: the path goes over the pole rather than bouncing off it.
    assert latitude(orbit, 0.5 * math.pi)[0] > 0 and latitude(orbit, 2.5 * math.pi)[0] < 0
    assert frequencies(orbit, table)[0] == pytest.approx((table.phi[-1] + math.pi) / (2 * math.pi))


def test_psi_zero_orbits_cross_the_equator_and_keep_their_colour():
    orbit = solve_orbit(0.6, 0.3, 1.0, 0.0)
    table = tabulate(orbit)
    assert orbit.w2 == 0
    # Z, the height above the equator ring, is positive, then negative from χ = π to 3π, then positive again.
    assert [embed(orbit, chi, 0.0)[2] > 0 for chi in (0.5 * math.pi, 1.5 * math.pi, 2.5 * math.pi, 3.5 * math.pi)] \
        == [True, False, False, True]
    psi = [node(table, i)[4] for i in range(2 * table.n + 1)]
    assert max(abs(p - q) for p, q in zip(psi, psi[1:])) < 0.05       # ψ, the colour, changes smoothly
