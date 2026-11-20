"""
Independent numerical check: integrate Hamilton's equations built directly from the
inverse of the metric (2.1) (sympy-inverted, NOT the paper's chi/Theta), starting on the
spherical orbit predicted by the polynomial solver, and compare with
 (i)  constancy of x, the predicted theta band, H = 0,
 (ii) the Mino-time / angle-variable formulation recommended for the demo.
Run: uv run --with numpy --with scipy --with sympy python integrate_check.py
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from mp5d import *

# ---------------- Hamiltonian from the metric (sympy)
t_, x_, th_, ph_, ps_ = sp.symbols('t x theta phi psi', real=True)
pt_, px_, pth_, pph_, pps_ = sp.symbols('p_t p_x p_theta p_phi p_psi', real=True)
mu_, a_, b_ = sp.symbols('mu a b', real=True)
s2, c2 = sp.sin(th_)**2, sp.cos(th_)**2
rho2 = x_ + a_**2*c2 + b_**2*s2
Delta = (x_ + a_**2)*(x_ + b_**2) - mu_*x_
kv = [1, a_*s2, b_*c2]
G = sp.Matrix(3, 3, lambda i, j: mu_/rho2*kv[i]*kv[j])
G[0, 0] += -1; G[1, 1] += (x_ + a_**2)*s2; G[2, 2] += (x_ + b_**2)*c2
Gi = G.inv(method='LU')
pv = sp.Matrix([pt_, pph_, pps_])
H = sp.Rational(1, 2)*((pv.T*Gi*pv)[0, 0] + 4*Delta/rho2*px_**2 + pth_**2/rho2)
vars_q = [t_, x_, th_, ph_, ps_]
vars_p = [pt_, px_, pth_, pph_, pps_]
dHdp = [sp.diff(H, p) for p in vars_p]
dHdq = [sp.diff(H, q) for q in (x_, th_)]
args = (x_, th_, pt_, px_, pth_, pph_, pps_, mu_, a_, b_)
fH = sp.lambdify(args, H, 'numpy')
fdHdp = sp.lambdify(args, dHdp, 'numpy')
fdHdq = sp.lambdify(args, dHdq, 'numpy')


def rhs_ham(lam, Y, a, b, mu):
    t, x, th, ph, ps, px, pth = Y
    E, Phi, Psi = 1.0, rhs_ham.Phi, rhs_ham.Psi
    A = (x, th, -E, px, pth, Phi, Psi, mu, a, b)
    qd = fdHdp(*A)
    pd = fdHdq(*A)
    return [qd[0], qd[1], qd[2], qd[3], qd[4], -pd[0], -pd[1]]


def ham_orbit(a, b, Phi, Psi, x0, th0, lam_max, mu=1.0, n=4000):
    rhs_ham.Phi, rhs_ham.Psi = Phi, Psi
    Y0 = [0.0, x0, th0, 0.0, 0.0, 0.0, 0.0]
    sol = solve_ivp(rhs_ham, (0, lam_max), Y0, args=(a, b, mu), method='DOP853',
                    rtol=1e-12, atol=1e-13, dense_output=True)
    lam = np.linspace(0, lam_max, n)
    return lam, sol.sol(lam)


# ---------------- angle-variable (chi) formulation in Mino time
def mino_orbit(a, b, Phi, Psi, o, n_periods, mu=1.0, N=20000, chi0=0.0):
    """Integrate d/dchi of (sigma, lambda, t, phi, psi) with y = ybar - dy cos(chi)."""
    x0, K = o['x'], o['K']
    y1, y2 = o['band']
    ybar, dy = 0.5*(y1 + y2), 0.5*(y2 - y1)
    Aphi, Apsi, At, calE = mino_constants(a, b, Phi, Psi, x0, mu)

    def f(chi, Z):
        y = ybar - dy*np.cos(chi)
        Gy = G_of_y(a, b, K, y1, y2, y)
        dsig = 1.0/(2*np.sqrt(Gy))
        rho2 = x0 + a*a + (b*b - a*a)*y
        dphi = (Phi/y if Phi != 0 else 0.0) - Aphi
        dpsi = (Psi/(1 - y) if Psi != 0 else 0.0) - Apsi
        return [dsig, rho2*dsig, (rho2 + At)*dsig, dphi*dsig, dpsi*dsig]
    chis = np.linspace(chi0, chi0 + 2*np.pi*n_periods, N)
    sol = solve_ivp(f, (chis[0], chis[-1]), [0, 0, 0, 0, 0], t_eval=chis, method='DOP853', rtol=1e-12, atol=1e-13)
    return chis, sol.y


def check_case(name, a, b, Phi, Psi, nper=3, mu=1.0, verbose=True):
    orbs = spherical_orbits(a, b, Phi, Psi, mu)
    assert len(orbs) == 1, orbs
    o = orbs[0]
    y1, y2 = o['band']
    K = o['K']
    Lsig = theta_period_sigma(a, b, K, y1, y2)
    # Mino-chi orbit
    chi0 = np.pi if (y1 == 0.0 and y2 < 1.0) else 0.0   # start at a turning point
    chis, Z = mino_orbit(a, b, Phi, Psi, o, nper, chi0=chi0)
    sig, lam, tt, ph, ps = Z
    lam_per = lam[-1]/nper
    # Hamiltonian orbit from theta = theta_min (turning point) for nper theta-periods
    th0 = np.arcsin(np.sqrt(y1)) if chi0 == 0.0 else np.arcsin(np.sqrt(y2))
    L, Y = ham_orbit(a, b, Phi, Psi, o['x'], th0, lam[-1], mu)
    t_h, x_h, th_h, ph_h, ps_h, px_h, pth_h = Y
    Hs = fH(x_h, th_h, -1.0, px_h, pth_h, Phi, Psi, mu, a, b)
    # compare at equal lambda: interpolate the chi solution
    th_m = np.arcsin(np.sqrt(np.clip(0.5*(y1+y2) - 0.5*(y2-y1)*np.cos(chis), 0, 1)))
    th_m_on_L = np.interp(L, lam, th_m)
    th_h = np.minimum(th_h % np.pi, np.pi - th_h % np.pi)   # fold the paper's [0,pi] range onto [0,pi/2]
    ph_m_on_L = np.interp(L, lam, ph); ps_m_on_L = np.interp(L, lam, ps); t_m_on_L = np.interp(L, lam, tt)
    res = dict(name=name, a=a, b=b, Phi=Phi, Psi=Psi, x=o['x'], r=o['r'], K=K,
               th_min_pred=np.degrees(np.arcsin(np.sqrt(y1))), th_max_pred=np.degrees(np.arcsin(np.sqrt(y2))),
               th_min_num=np.degrees(th_h.min()), th_max_num=np.degrees(th_h.max()),
               max_dx=np.abs(x_h - o['x']).max(), max_H=np.abs(Hs).max(),
               dth=np.abs(th_h - th_m_on_L).max(),
               dph=np.abs(ph_h - ph_m_on_L).max(),
               dps=np.abs(ps_h - ps_m_on_L).max(),
               dt=np.abs(t_h - t_m_on_L).max(),
               Lsig=Lsig, Lsig_num=sig[-1]/nper, lam_per=lam_per,
               dphi_per=ph[-1]/nper, dpsi_per=ps[-1]/nper, dt_per=tt[-1]/nper)
    if verbose:
        print(f"--- {name}: a={a} b={b} Phi={Phi} Psi={Psi}")
        print(f"    x={o['x']:.10f} r={o['r']:.10f} K={K:.10f} chi''={o['chi2']:.4f}")
        print(f"    theta band predicted [{res['th_min_pred']:.6f}, {res['th_max_pred']:.6f}] deg;"
              f" Hamiltonian integration [{res['th_min_num']:.6f}, {res['th_max_num']:.6f}] deg")
        print(f"    max|x - x0| over {nper} theta-periods (lambda in [0,{lam[-1]:.3f}]) = {res['max_dx']:.2e}; max|H| = {res['max_H']:.2e}")
        print(f"    Ham vs Mino-chi: max|dtheta|={res['dth']:.2e} max|dphi|={res['dph']:.2e} max|dpsi|={res['dps']:.2e} max|dt|={res['dt']:.2e}")
        print(f"    Mino period: elliptic K formula {Lsig:.10f}, quadrature {sig[-1]/nper:.10f}; lambda-period {lam_per:.6f}")
        print(f"    per theta-period: dphi={res['dphi_per']:.6f} dpsi={res['dpsi_per']:.6f} dt={res['dt_per']:.6f}")
    return res


if __name__ == "__main__":
    cases = [
        ("Tangherlini", 0.0, 0.0, 0.7, 0.9),
        ("Fig4.3 b=0", 0.5, 0.0, 2.0, 0.0),
        ("Fig4.7", 0.6, 0.3, 1.0, -1.0),
        ("Fig4.8 (extremal)", 0.9, 0.1, 0.8, 0.8),
        ("Fig4.10", 0.9, 0.09, 0.1, 0.01),
        ("generic a>b", 0.55, 0.2, -0.6, 0.4),
        ("generic b>a", 0.15, 0.7, 0.9, -0.6),
        ("a=b", 0.4, 0.4, 1.2, -0.5),
        ("negative b", 0.5, -0.3, 0.8, 1.0),
        ("closed 15-osc orbit", 0.6, 0.3, 0.1051348934, -0.1354415960),
        ("Phi=0 polar-crossing", 0.8, 0.0, 0.0, 0.5),
    ]
    for c in cases:
        check_case(*c)
