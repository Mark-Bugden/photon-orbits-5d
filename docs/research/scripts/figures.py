"""
Reproduce the paper's Figures 4.2-4.12 with the angle-variable (chi) / Mino-time
formulation, report the orbit data for the preset table, and render the paper's
projection (x1, y1, z1) coloured by psi for visual comparison (scratch PNGs).
Run: uv run --with numpy --with scipy --with matplotlib --with sympy python figures.py
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from mp5d import *


def orbit_chi(a, b, Phi, Psi, o, chi_max, N=6000, mu=1.0, chi0=0.0):
    """Angle-variable integration. Returns dict of arrays sampled uniformly in chi.
    y = ybar - dy cos(chi);  signed sin/cos theta used when Phi = 0 / Psi = 0."""
    x0, K = o['x'], o['K']
    y1, y2 = o['band']
    ybar, dy = 0.5*(y1 + y2), 0.5*(y2 - y1)
    Aphi, Apsi, At, calE = mino_constants(a, b, Phi, Psi, x0, mu)

    def f(chi, Z):
        y = ybar - dy*np.cos(chi)
        dsig = 1.0/(2*np.sqrt(G_of_y(a, b, K, y1, y2, y)))
        rho2 = x0 + a*a + (b*b - a*a)*y
        dphi = (Phi/y if Phi != 0 else 0.0) - Aphi
        dpsi = (Psi/(1 - y) if Psi != 0 else 0.0) - Apsi
        return [dsig, rho2*dsig, (rho2 + At)*dsig, dphi*dsig, dpsi*dsig]
    chis = np.linspace(chi0, chi0 + chi_max, N)
    sol = solve_ivp(f, (chis[0], chis[-1]), [0, 0, 0, 0, 0], t_eval=chis, method='DOP853',
                    rtol=1e-11, atol=1e-12)
    sig, lam, t, ph, ps = sol.y
    y = ybar - dy*np.cos(chis)
    sth = np.sqrt(y2)*np.sin(chis/2) if Phi == 0 else np.sqrt(np.clip(y, 0, 1))
    cth = np.sqrt(1 - y1)*np.cos(chis/2) if Psi == 0 else np.sqrt(np.clip(1 - y, 0, 1))
    return dict(chi=chis, sig=sig, lam=lam, t=t, phi=ph, psi=ps, y=y, sth=sth, cth=cth,
                Aphi=Aphi, Apsi=Apsi, At=At)


def solve_constant_theta(a, b, theta0, guess, mu=1.0):
    """circular orbit at fixed theta0: unknowns (x, K, Phi, Psi)"""
    y0 = np.sin(theta0)**2

    def eqs(v):
        x, K, Phi, Psi = v
        chi = chi_coef(a, b, Phi, Psi, K, mu)
        Q = Q_coef(a, b, Phi, Psi, K)
        return [np.polyval(chi, x), np.polyval(np.polyder(chi), x),
                np.polyval(Q, y0), np.polyval(np.polyder(Q), y0)]
    sol, info, ier, msg = fsolve(eqs, guess, full_output=True, xtol=1e-14)
    return sol, ier, np.max(np.abs(eqs(sol)))


figs = [
    # name, a, b, Phi, Psi, chi-periods to draw, caption
    ("4.2", 0.5, 0.0, 0.0, 0.0, 2, "one latitudinal oscillation, zero momenta"),
    ("4.3", 0.5, 0.0, 2.0, 0.0, 60, "many latitudinal oscillations, Phi=2"),
    ("4.4", 0.8, 0.0, 0.0, 0.5, 2, "Psi=0.5, Phi=0 (short)"),
    ("4.5", 0.8, 0.0, 0.0, 0.5, 60, "Psi=0.5, Phi=0 (many)"),
    ("4.6", 0.0, 0.99, 0.0, 0.0, 2, "football, zero momenta"),
    ("4.7", 0.6, 0.3, 1.0, -1.0, 40, "band, Phi=1, Psi=-1"),
    ("4.8", 0.9, 0.1, 0.8, 0.8, None, "cusps, extremal BH, lambda in [0,5]"),
    ("4.10", 0.9, 0.09, 0.1, 0.01, None, "non-monotonic phi/psi, lambda in [1,5]"),
]

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import hsv_to_rgb
    rows = []
    fig = plt.figure(figsize=(20, 10))
    for k, (name, a, b, Phi, Psi, nper, cap) in enumerate(figs):
        o = spherical_orbits(a, b, Phi, Psi)[0]
        y1, y2 = o['band']
        Lsig = theta_period_sigma(a, b, o['K'], y1, y2)
        one = orbit_chi(a, b, Phi, Psi, o, 2*np.pi, N=4001)
        lam_per, t_per, dphi, dpsi = one['lam'][-1], one['t'][-1], one['phi'][-1], one['psi'][-1]
        # monotonicity of phi, psi along the orbit (Mino-time rates at band edges)
        rphi = [((Phi/yy) if Phi != 0 else 0) - one['Aphi'] for yy in (max(y1, 1e-300), y2)]
        rpsi = [((Psi/(1-yy)) if Psi != 0 else 0) - one['Apsi'] for yy in (y1, min(y2, 1-1e-16))]
        print(f"Fig {name}: a={a} b={b} Phi={Phi} Psi={Psi} | x={o['x']:.6f} r={o['r']:.6f} K={o['K']:.6f} "
              f"| theta in [{np.degrees(o['theta_min']):.3f}, {np.degrees(o['theta_max']):.3f}] deg | x+={horizon(a,b):.5f}")
        print(f"      per chi-period: dlambda={lam_per:.5f} dt={t_per:.5f} dphi={dphi:.5f} dpsi={dpsi:.5f} dsigma={Lsig:.5f}"
              f" | dphi/dsig at (y1,y2)=({rphi[0]:.4g},{rphi[1]:.4g}) dpsi/dsig at (y1,y2)=({rpsi[0]:.4g},{rpsi[1]:.4g}) | chi''={o['chi2']:.3f}")
        rows.append((name, a, b, Phi, Psi, o, lam_per, t_per, dphi, dpsi))
        if nper is None:
            nper = 5.0/lam_per + 1
        orb = orbit_chi(a, b, Phi, Psi, o, 2*np.pi*nper, N=int(3000*nper) if nper < 5 else 30000)
        if name == "4.8":     # fitted phase from fit_figs.py: start just after theta_min
            pass
        X = np.sqrt(o['x'] + a*a)*orb['sth']*np.cos(orb['phi'])
        Y = np.sqrt(o['x'] + a*a)*orb['sth']*np.sin(orb['phi'])
        Z = np.sqrt(o['x'] + b*b)*orb['cth']
        hue = (orb['psi'] % (2*np.pi))/(2*np.pi)
        ax = fig.add_subplot(2, 5, k + 1, projection='3d')
        ax.scatter(X, Y, Z, c=hsv_to_rgb(np.c_[hue, np.ones_like(hue), np.ones_like(hue)]), s=0.3)
        u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
        ax.plot_wireframe(np.sqrt(o['x']+a*a)*np.cos(u)*np.sin(v), np.sqrt(o['x']+a*a)*np.sin(u)*np.sin(v),
                          np.sqrt(o['x']+b*b)*np.cos(v), color='0.85', lw=0.3)
        ax.set_title(f"Fig {name}: {cap}", fontsize=8)
        ax.set_box_aspect((1, 1, np.sqrt(o['x']+b*b)/np.sqrt(o['x']+a*a)))
    # Fig 4.12: constant theta = pi/4
    print("\nFig 4.12: circular orbit at theta = pi/4, a=0.5, b=0.25")
    sols = []
    for gPhi in (-1.0, 1.0):
        for gPsi in (-1.0, 1.0):
            for gx in (0.8, 1.5, 2.5):
                s, ier, res = solve_constant_theta(0.5, 0.25, np.pi/4, [gx, 3.0, gPhi, gPsi])
                if ier == 1 and res < 1e-10 and s[0] > horizon(0.5, 0.25):
                    if not any(np.allclose(s, q, atol=1e-8) for q in sols):
                        sols.append(s)
    for s in sols:
        o = spherical_orbits(0.5, 0.25, s[2], s[3])
        print(f"   x={s[0]:.6f} r={np.sqrt(s[0]):.6f} K={s[1]:.6f} Phi={s[2]:+.6f} Psi={s[3]:+.6f}  | re-solved from (Phi,Psi): "
              f"x={o[0]['x']:.6f}, theta band [{np.degrees(o[0]['theta_min']):.4f},{np.degrees(o[0]['theta_max']):.4f}] deg")
    o = spherical_orbits(0.5, 0.25, -0.9292, 0.9036)
    print("   paper's rounded values Phi=-0.9292, Psi=0.9036 ->", [(round(q['x'], 6), round(q['K'], 6), np.degrees(q['theta_min']), np.degrees(q['theta_max'])) for q in o])
    plt.tight_layout(); plt.savefig('figures_repro.png', dpi=100)
    print("saved figures_repro.png")
