"""
Map the set of (Phi/E, Psi/E) for which a spherical photon orbit exists, and test
uniqueness of (x, K) given (Phi, Psi). mu = 1, E = 1.
Run: uv run --with numpy --with scipy python explore_region.py
"""
import numpy as np
from mp5d import *

cases = [(0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.3, 0.3), (0.5, 0.5), (0.6, 0.3), (0.3, 0.6),
         (0.5, 0.25), (0.8, 0.0), (0.9, 0.09), (0.89, 0.1), (0.6, -0.3), (0.99, 0.0), (0.495, 0.495),
         (0.7, 0.2999), (0.1, 0.85)]
N = 241
L = 4.6
grid = np.linspace(-L, L, N)

print(f"grid {N}x{N} on [-{L},{L}]^2")
print(" a      b      x+      #crit(F>x+) max | #valid max | frac valid | min chi''| min dt/dsig | max |x| range of valid orbits")
for a, b in cases:
    xp = horizon(a, b)
    maxcrit = 0; maxvalid = 0; nvalid = 0; minchi2 = np.inf; mindt = np.inf
    xs = []
    multi = []
    for Phi in grid:
        for Psi in grid:
            orbs = spherical_orbits(a, b, Phi, Psi, require_band=False)
            maxcrit = max(maxcrit, len(orbs))
            val = [o for o in orbs if o['band'] is not None]
            maxvalid = max(maxvalid, len(val))
            if len(val) > 1:
                multi.append((Phi, Psi, [(o['x'], o['K']) for o in val]))
            if val:
                nvalid += 1
                o = val[0]
                xs.append(o['x'])
                minchi2 = min(minchi2, o['chi2'])
                Aphi, Apsi, At, calE = mino_constants(a, b, Phi, Psi, o['x'])
                for yy in o['band']:
                    rho2 = o['x'] + a*a + (b*b - a*a)*yy
                    mindt = min(mindt, rho2 + At)
    print(f"{a:5.3f} {b:6.3f} {xp:7.4f}   {maxcrit:3d}            | {maxvalid:3d}       | {nvalid/N**2:6.3f}   | {minchi2:8.4f} | {mindt:9.4f} | [{min(xs):.4f}, {max(xs):.4f}]")
    if multi:
        print("    multiple valid orbits at e.g.", multi[:3])
