"""
Shape of the allowed (Phi, Psi) region (mu = E = 1): star-shapedness, convexity,
corners vs the closed-form equatorial/polar circular orbits, x-range, and plots.
Run: uv run --with numpy --with scipy --with matplotlib python region.py
"""
import numpy as np
from mp5d import *


def valid(a, b, Phi, Psi):
    o = spherical_orbits(a, b, Phi, Psi)
    return o[0] if o else None


def boundary_s(a, b, alpha, smax=8.0, tol=1e-11):
    """largest s with a valid orbit at (s cos alpha, s sin alpha), by bisection"""
    lo, hi = 0.0, smax
    assert valid(a, b, 0.0, 0.0) is not None
    while hi - lo > tol:
        m = 0.5*(lo + hi)
        if valid(a, b, m*np.cos(alpha), m*np.sin(alpha)) is not None:
            lo = m
        else:
            hi = m
    return lo


def corners(a, b):
    Sp = np.sqrt((1 + a)**2 - b**2); Sm = np.sqrt((1 - a)**2 - b**2)
    Tp = np.sqrt((1 + b)**2 - a**2); Tm = np.sqrt((1 - b)**2 - a**2)
    return dict(Phi_plus=(1 + Sp, 1 + a - b*b + Sp), Phi_minus=(-1 - Sm, 1 - a - b*b + Sm),
                Psi_plus=(1 + Tp, 1 + b - a*a + Tp), Psi_minus=(-1 - Tm, 1 - b - a*a + Tm))


if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cases = [(0.0, 0.0), (0.5, 0.0), (0.49, 0.49), (0.6, 0.3), (0.9, 0.09), (0.3, 0.6), (0.7, 0.2), (0.4, 0.4)]
    fig, axs = plt.subplots(2, 4, figsize=(18, 9))
    for (a, b), ax in zip(cases, axs.ravel()):
        print(f"=== a={a}, b={b}, x+={horizon(a, b):.6f}")
        c = corners(a, b)
        for k, (val, xr) in c.items():
            alpha = {'Phi_plus': 0.0, 'Psi_plus': np.pi/2, 'Phi_minus': np.pi, 'Psi_minus': 1.5*np.pi}[k]
            s = boundary_s(a, b, alpha)
            o = valid(a, b, 0.999999*s*np.cos(alpha), 0.999999*s*np.sin(alpha))
            print(f"  {k:9s}: closed form {val:+.8f} (x={xr:.6f}) | ray bisection {np.sign(np.cos(alpha)+np.sin(alpha))*s:+.8f} (x={(o['x'] if o else float('nan')):.6f})")
        alphas = np.linspace(0, 2*np.pi, 361)[:-1]
        S = np.array([boundary_s(a, b, al, tol=1e-7) for al in alphas])
        # star-shapedness: everything inside along the ray valid, slightly outside invalid
        bad = 0
        for al, s in zip(alphas[::6], S[::6]):
            for f in np.linspace(0, 0.995, 40):
                if valid(a, b, f*s*np.cos(al), f*s*np.sin(al)) is None:
                    bad += 1
            for f in np.linspace(1.005, 2.0, 40):
                if valid(a, b, f*s*np.cos(al), f*s*np.sin(al)) is not None:
                    bad += 1
        # boundary orbits have degenerate theta band
        degen = []
        for al, s in zip(alphas[::15], S[::15]):
            o = valid(a, b, s*np.cos(al), s*np.sin(al))
            if o is not None:
                degen.append(np.degrees(o['theta_max'] - o['theta_min']))
        # convexity test: midpoints of random boundary point pairs inside?
        P = np.c_[S*np.cos(alphas), S*np.sin(alphas)]
        rng = np.random.default_rng(1)
        nonconv = 0
        for _ in range(300):
            i, j = rng.integers(0, len(P), 2)
            m = 0.5*(P[i] + P[j])*0.9999
            if valid(a, b, *m) is None:
                nonconv += 1
        # x range over a fine interior grid
        xs = []
        for al, s in zip(alphas[::4], S[::4]):
            for f in np.linspace(0, 0.999, 15):
                o = valid(a, b, f*s*np.cos(al), f*s*np.sin(al))
                xs.append(o['x'])
        print(f"  star-shaped w.r.t. origin violations: {bad};  convexity violations (300 random chords): {nonconv}")
        print(f"  theta-band width on boundary (deg): max {max(degen):.2e};  x range over region [{min(xs):.5f}, {max(xs):.5f}]")
        print(f"  extent: Phi in [{P[:,0].min():.5f}, {P[:,0].max():.5f}], Psi in [{P[:,1].min():.5f}, {P[:,1].max():.5f}]")
        # plot: x as colour on a grid
        g = np.linspace(-4.2, 4.2, 161)
        Xg = np.full((len(g), len(g)), np.nan)
        for i, Ps in enumerate(g):
            for j, Ph in enumerate(g):
                o = valid(a, b, Ph, Ps)
                if o is not None:
                    Xg[i, j] = o['x']
        im = ax.pcolormesh(g, g, Xg, shading='auto', cmap='viridis')
        ax.plot(np.r_[P[:, 0], P[0, 0]], np.r_[P[:, 1], P[0, 1]], 'r-', lw=1)
        ax.axhline(0, color='k', lw=0.4); ax.axvline(0, color='k', lw=0.4)
        ax.set_aspect('equal'); ax.set_title(f"a={a}, b={b}  (colour: x=r^2)")
        ax.set_xlabel('Phi/E'); ax.set_ylabel('Psi/E')
        plt.colorbar(im, ax=ax, shrink=0.7)
    plt.tight_layout()
    plt.savefig('region_maps.png', dpi=110)
    print("saved region_maps.png")
