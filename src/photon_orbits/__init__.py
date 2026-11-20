"""Reference model of spherical photon orbits around a 5D Myers-Perry black hole (Bugden 2018, arXiv:1801.03389).

The demo's physics, demo/src/physics/, is a port of this package, and its tests replay tests/fixtures/orbits.json,
which fixtures.py writes from it. docs/equations.md sets out the equations as implemented here.
"""

from .orbit import NoOrbit, Orbit, SecondRootWarning, horizon, orbit_exists, solve_orbit
from .oscillation import Table, embed, frequencies, latitude, node, tabulate
from .region import boundary, directions, edge, vertices

__all__ = [
    "NoOrbit", "Orbit", "SecondRootWarning", "Table",
    "boundary", "directions", "edge", "embed", "frequencies", "horizon", "latitude", "node", "orbit_exists",
    "solve_orbit", "tabulate", "vertices",
]
