"""The cross-check against the full metric: integrate null geodesics of the corrected eq. (2.1) in Boyer-Lindquist
coordinates, from Hamilton's equations for H = ½ g^{μν} p_μ p_ν, and compare them with the reference model.

The inverse metric and H come from sympy, straight from the line element, so they share nothing with 𝒳, Θ, the
Mino-time rates or anything else in photon_orbits. The integration also holds x free, so it shows the radial
instability: started on the orbit, it stays there for about 3 oscillations and has left by about 7.
"""

import math

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import solve_ivp

from photon_orbits import latitude, node, solve_orbit, tabulate


@pytest.fixture(scope="module")
def hamilton():
    """(H, rates): H and the right-hand sides d(x, θ, p_x, p_θ, φ, ψ)/dλ, with E = 1, p_t = −1, μ = 1."""
    x, th, px, pth, a, b, Phi, Psi = sp.symbols("x theta p_x p_theta a b Phi Psi", real=True)
    s2, c2 = sp.sin(th) ** 2, sp.cos(th) ** 2
    rho2 = x + a**2 * c2 + b**2 * s2
    delta = (x + a**2) * (x + b**2) - x
    # ds² = −dt² + (μ/ρ²)(dt + a sin²θ dφ + b cos²θ dψ)² + ρ²/(4Δ) dx² + ρ² dθ² + (x + a²) sin²θ dφ² + (x + b²) cos²θ dψ²
    k = [1, a * s2, b * c2]
    block = sp.Matrix(3, 3, lambda i, j: k[i] * k[j] / rho2)
    block[0, 0] -= 1
    block[1, 1] += (x + a**2) * s2
    block[2, 2] += (x + b**2) * c2
    p = sp.Matrix([-1, Phi, Psi])
    H = ((p.T * block.inv(method="LU") * p)[0, 0] + 4 * delta / rho2 * px**2 + pth**2 / rho2) / 2
    args = (x, th, px, pth, a, b, Phi, Psi)
    rates = [sp.diff(H, px), sp.diff(H, pth), -sp.diff(H, x), -sp.diff(H, th), sp.diff(H, Phi), sp.diff(H, Psi)]
    return sp.lambdify(args, H, "math", cse=True), sp.lambdify(args, rates, "math", cse=True)


def integrate(hamilton, orbit, lam_end, **options):
    """From θ_min, at rest in x and θ, with φ = ψ = 0: the same start as the table's first node."""
    _, rates = hamilton
    params = (orbit.a, orbit.b, orbit.Phi, orbit.Psi)
    start = [orbit.x, math.asin(math.sqrt(orbit.y1)), 0.0, 0.0, 0.0, 0.0]
    return solve_ivp(lambda lam, Y: rates(*Y[:4], *params), (0, lam_end), start, method="DOP853",
                     rtol=1e-12, atol=1e-13, dense_output=True, **options)


@pytest.mark.parametrize("a, b, Phi, Psi", [
    (0.5, 0.0, 2.0, 0.0),           # Fig 4.3: Ψ = 0, so the path crosses the equator
    (0.6, 0.3, 1.0, -1.0),          # Fig 4.7
    (0.9, 0.1, 0.8, 0.8),           # Fig 4.8: extremal
    (0.9, 0.09, 0.1, 0.01),         # Fig 4.10: close to both axes
    (0.15, 0.7, 0.9, -0.6),         # b > a
    (0.9, 0.1, -0.5, -0.49),        # extremal, 0.01 inside the co-rotating edge: 0.1 outside the horizon
])
def test_the_table_follows_the_geodesic(hamilton, a, b, Phi, Psi):
    H, _ = hamilton
    orbit = solve_orbit(a, b, Phi, Psi)
    table = tabulate(orbit)
    nodes = [node(table, i) for i in range(0, 3 * table.n + 1, 4)]
    path = integrate(hamilton, orbit, nodes[-1][2])
    assert path.success
    x, theta, px, pth, phi, psi = path.sol(np.array([lam for _, _, lam, _, _ in nodes]))
    params = (a, b, Phi, Psi)
    assert max(abs(H(*y, *params)) for y in zip(x, theta, px, pth)) < 1e-10   # a null geodesic throughout
    assert np.max(np.abs(x - orbit.x)) < 1e-7                                # the radius holds for 3 oscillations
    ours = np.array([math.atan2(*latitude(orbit, chi)) for chi, *_ in nodes])
    assert np.max(np.abs(theta - ours)) < 1e-7
    assert np.max(np.abs(phi - [n[3] for n in nodes])) < 1e-6
    assert np.max(np.abs(psi - [n[4] for n in nodes])) < 1e-6


@pytest.mark.parametrize("a, b, Phi, Psi", [(0.6, 0.3, 1.0, -1.0), (0.0, 0.0, 0.7, 0.9)])
def test_every_orbit_is_radially_unstable(hamilton, a, b, Phi, Psi):
    """Rounding errors grow by about e⁴ an oscillation: x leaves the orbit by 0.1% within 3 to 9 oscillations."""
    orbit = solve_orbit(a, b, Phi, Psi)
    per = tabulate(orbit).lam[-1]
    leaves = lambda lam, Y: abs(Y[0] - orbit.x) - 1e-3 * orbit.x  # noqa: E731
    leaves.terminal = True
    path = integrate(hamilton, orbit, 12 * per, events=leaves)
    assert path.status == 1
    assert 3 < path.t[-1] / per < 9
