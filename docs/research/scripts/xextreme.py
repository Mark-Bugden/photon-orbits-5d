"""Where are min/max of x over the allowed region? Sample the boundary finely (x is extremal on the boundary)."""
import numpy as np
from region import boundary_s, valid, corners
for a, b in [(0.6, 0.3), (0.3, 0.6), (0.9, 0.09), (0.5, 0.0), (0.4, 0.4)]:
    al = np.linspace(0, 2*np.pi, 721)[:-1]
    best = []
    for A in al:
        s = boundary_s(a, b, A, tol=1e-9)
        o = valid(a, b, 0.9999999*s*np.cos(A), 0.9999999*s*np.sin(A))
        best.append((o['x'], A, s))
    xs = np.array([q[0] for q in best])
    imin, imax = xs.argmin(), xs.argmax()
    c = corners(a, b)
    print(f"a={a} b={b}: x_min={xs[imin]:.6f} at alpha={np.degrees(best[imin][1]):.1f} deg; x_max={xs[imax]:.6f} at alpha={np.degrees(best[imax][1]):.1f} deg;"
          f" vertices x: Phi+ {c['Phi_plus'][1]:.6f}, Phi- {c['Phi_minus'][1]:.6f}, Psi+ {c['Psi_plus'][1]:.6f}, Psi- {c['Psi_minus'][1]:.6f}")
