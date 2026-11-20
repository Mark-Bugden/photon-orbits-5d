"""Writes tests/fixtures/orbits.json: results of the reference model that the demo's JS port must reproduce.

demo/test/physics.test.js replays every case to 1e-9, and tests/test_fixtures.py checks that this file still
matches the model. Keys use the JS names. Run with `uv run python -m photon_orbits.fixtures` after any change
to the model.
"""

import json
from pathlib import Path

from .orbit import horizons, solve_orbit
from .oscillation import embed, frequencies, node, tabulate
from .region import directions, edge, vertices

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "tests" / "fixtures" / "orbits.json"

EVERY = 8           # the table's nodes are kept every EVERY steps (256 and 512 are multiples of it)

ORBITS = [
    # the paper review's unit-test table (§B7)
    (0, 0, 0, 0), (0, 0, 0.7, 0.9), (0.5, 0, 0, 0), (0.5, 0, 2, 0), (0.8, 0, 0, 0.5), (0, 0.99, 0, 0),
    (0.4, 0.4, 1.2, -0.5), (0.6, 0.3, 1, -1), (0.9, 0.1, 0.8, 0.8), (0.9, 0.09, 0.1, 0.01), (0.55, 0.2, -0.6, 0.4),
    (0.15, 0.7, 0.9, -0.6), (0.5, 0.25, -0.929241, 0.903667), (0.6, 0.3, 0, 0), (0.6, -0.3, 1, 1),
    # the other circular orbits at 45°, the closed orbit, and band edges close to an axis
    (0.5, 0.25, -0.764662, -0.733371), (0.5, 0.25, 1.048754, -1.026162), (0.5, 0.25, 1.168941, 1.148716),
    (0.6, 0.3, 0.1051348934, -0.135441596), (0.6, 0.3, 1e-3, 0.5), (0.6, 0.3, -0.8, -1e-3),
    # tiny momenta: snapped to 0 below 1e-12, kept exact above
    (0.6, 0.3, 1e-16, 0.5), (0.6, 0.3, 0.5, -1e-16), (0.6, 0.3, 1e-12, 0.5), (0.6, 0.3, 0.5, 1e-12),
    (0.9, 0.09, 1e-10, -1e-10), (0.6, 0.3, 1.0, 0.0),
    # extremal holes, near the co-rotating edge Φ + Ψ = −1
    (0.9, 0.1, -0.5, -0.49), (0.5, 0.5, -0.3, -0.6999), (0.9, 0.1, -0.999, 0.0), (0.7, 0.3, 0.0, -0.99999),
]

NO_ORBITS = [
    (0.4, 0.4, 1.5, -1.5), (0.55, 0.2, -1.3, 0.6), (0, 0, 1.5, 0.6),       # outside the region: no band
    (0.5, 0.5, -0.5, -0.5), (0.9, 0.1, -0.6, -0.5),                        # on and past the extremal edge
    (0.6, 0.6, 0.1, 0.1), (1.0, 0.0, 0.1, 0.1),                            # no horizon
]

REGIONS = [(0, 0), (0.5, 0), (0.6, 0.3), (0.4, 0.4), (0.9, 0.09), (0.9, 0.1), (0.5, 0.5), (0.3, 0.699)]

HORIZONS = [(0, 0), (0.5, 0), (0.6, 0.3), (0.9, 0.1), (0.5, 0.5), (0.99, 0.01), (0.6, 0.6), (1.0, 0.0), (0.0, 1.0)]


def orbit_case(a, b, Phi, Psi):
    orbit = solve_orbit(a, b, Phi, Psi)
    table = tabulate(orbit)
    n = table.n
    points = []
    for i in (0, n // 4, n // 2, 3 * n // 4, n, n + n // 3, 2 * n + 5, -7, 37 * n + 11):
        chi, sigma, lam, phi, psi = node(table, i)
        X, Y, Z = embed(orbit, chi, phi)
        points.append({"i": i, "chi": chi, "sigma": sigma, "lam": lam, "phi": phi, "psi": psi, "X": X, "Y": Y, "Z": Z})
    return {
        "a": a, "b": b, "Phi": Phi, "Psi": Psi,
        "orbit": {
            "Phi": orbit.Phi, "Psi": orbit.Psi, "xPlus": orbit.x_plus, "x": orbit.x, "K": orbit.K,
            "y1": orbit.y1, "y2": orbit.y2, "w2": orbit.w2, "delta": orbit.delta, "calE": orbit.calE,
            "Aphi": orbit.A_phi, "Apsi": orbit.A_psi, "G1": orbit.G1, "G2": orbit.G2, "extraRoots": orbit.extra_roots,
        },
        "table": {
            "n": n, "every": EVERY,
            "sigma": list(table.sigma[::EVERY]), "lam": list(table.lam[::EVERY]),
            "phi": list(table.phi[::EVERY]), "psi": list(table.psi[::EVERY]),
        },
        "frequencies": list(frequencies(orbit, table)),
        "points": points,
    }


def no_orbit_case(a, b, Phi, Psi):
    found = solve_orbit(a, b, Phi, Psi)
    return {"a": a, "b": b, "Phi": Phi, "Psi": Psi, "reason": found.reason, "x": found.x, "K": found.K}


def region_case(a, b):
    rays = directions(96)
    return {
        "a": a, "b": b,
        "vertices": [list(v) for v in vertices(a, b)],
        "directions": [list(d) for d in rays],
        "edges": [edge(a, b, *d) for d in rays],
    }


def horizon_case(a, b):
    found = horizons(a, b)
    return {"a": a, "b": b, "xPlus": None if found is None else found[0], "e": None if found is None else found[1]}


def build():
    return {
        "about": "Results of the reference model (src/photon_orbits/) for the JS port, demo/src/physics/, to "
                 "reproduce. Tables keep every 'every'-th node. Written by src/photon_orbits/fixtures.py.",
        "horizons": [horizon_case(*case) for case in HORIZONS],
        "orbits": [orbit_case(*case) for case in ORBITS],
        "noOrbits": [no_orbit_case(*case) for case in NO_ORBITS],
        "regions": [region_case(*case) for case in REGIONS],
    }


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(build(), separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"{OUT.relative_to(ROOT)}: {OUT.stat().st_size / 1000:.0f} KB")


if __name__ == "__main__":
    main()
