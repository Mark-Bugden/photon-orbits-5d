"""
Core numerics for spherical photon orbits of the 5D Myers-Perry black hole.

Conventions (same as Bugden 2018 / Frolov-Stojkovic 2003):
  metric (2.1) with x = r^2, rho^2 = x + a^2 cos^2 th + b^2 sin^2 th,
  Delta = (x+a^2)(x+b^2) - mu x,  p_t = -E, p_phi = Phi, p_psi = Psi.
  E = 1 (so Phi, Psi, K are Phi/E, Psi/E, K/E^2), mu = 1 by default.
  theta in [0, pi/2]; y = sin^2 theta in [0, 1].

The radial function is chi(x) = P(x) - K Delta(x) with P cubic (leading coeff 1).
Spherical orbits = critical points of V(x) = P(x)/Delta(x) on x > x_+, with K = V(x0).
"""
import numpy as np
from numpy.polynomial import polynomial as npoly
from scipy.optimize import brentq
from scipy.special import ellipk

TWO_PI = 2*np.pi


def horizon(a, b, mu=1.0):
    s = mu - a*a - b*b
    disc = s*s - 4*a*a*b*b
    if -1e-12 < disc < 0:
        disc = 0.0                      # extremal (round-off)
    if s <= 0 or disc < 0:
        return None
    return 0.5*(s + np.sqrt(disc))


def delta_coef(a, b, mu=1.0):            # highest power first
    return np.array([1.0, a*a + b*b - mu, a*a*b*b])


def P_coef(a, b, Phi, Psi, mu=1.0):
    A2, B2 = a*a, b*b
    P1 = (Phi**2 - Psi**2)*(A2 - B2) + A2*B2 + mu*(A2 + B2 + 2*a*Phi + 2*b*Psi)
    P0 = mu*(a*b + b*Phi + a*Psi)**2 + (A2 - B2)*(B2*Phi**2 - A2*Psi**2)
    return np.array([1.0, A2 + B2, P1, P0])


def F_coef(a, b, Phi, Psi, mu=1.0):
    """F = P' Delta - P Delta' (quartic, leading coeff 1). Roots = spherical-orbit radii."""
    P = P_coef(a, b, Phi, Psi, mu)
    D = delta_coef(a, b, mu)
    return np.polysub(np.polymul(np.polyder(P), D), np.polymul(P, np.polyder(D)))


def chi_coef(a, b, Phi, Psi, K, mu=1.0):
    return np.polysub(P_coef(a, b, Phi, Psi, mu), K*np.concatenate([[0.0], delta_coef(a, b, mu)]))


def Q_coef(a, b, Phi, Psi, K):
    """Q(y) = y(1-y) Theta, y = sin^2 theta (E = 1); highest power first."""
    A2, B2 = a*a, b*b
    return np.array([A2 - B2, B2 - 2*A2 - K, Phi**2 - Psi**2 + A2 + K, -Phi**2])


def theta_band(a, b, Phi, Psi, K, tol=1e-12):
    """Return (y1, y2), y = sin^2 theta, with Theta >= 0 exactly on [y1, y2], or None.

    Theta = Q(y)/(y(1-y)).  If Phi = 0 (Psi = 0) the factor y (1-y) is divided out
    first, because then Theta at y = 0 (y = 1) is finite and equals R(0) (R(1)).
    Since Q(0) = -Phi^2 <= 0 and Q(1) = -Psi^2 <= 0 the band is unique.
    y1 == y2 means a circular orbit at constant theta."""
    R = Q_coef(a, b, Phi, Psi, K)
    if Phi == 0:
        R, _ = np.polydiv(R, np.array([1.0, 0.0]))       # divide by y
    if Psi == 0:
        R, _ = np.polydiv(R, np.array([-1.0, 1.0]))      # divide by (1 - y)
    R = np.atleast_1d(R)
    # endpoints can only carry Theta >= 0 if the corresponding momentum vanishes exactly
    cands = ([0.0] if Phi == 0 else []) + ([1.0] if Psi == 0 else [])
    if len(R) > 2:
        for r in np.roots(np.polyder(R)):
            if abs(r.imag) < 1e-12 and 0.0 < r.real < 1.0:
                cands.append(r.real)
    if not cands:
        return None
    vals = [np.polyval(R, c) for c in cands]
    i = int(np.argmax(vals))
    ym, Rm = cands[i], vals[i]
    scale = 1 + abs(K) + Phi**2 + Psi**2
    if Rm < -tol*scale:
        return None
    if Rm <= tol*scale:          # double root: circular (constant-theta) orbit
        return (ym, ym)
    f = lambda yy: np.polyval(R, yy)
    y1 = 0.0 if f(0.0) >= 0 else brentq(f, 0.0, ym, xtol=1e-15, rtol=1e-15)
    y2 = 1.0 if f(1.0) >= 0 else brentq(f, ym, 1.0, xtol=1e-15, rtol=1e-15)
    return (y1, y2)


def K_theta_min(a, b, Phi, Psi):
    """Smallest K for which theta-motion is possible: min_y [Phi^2/y + Psi^2/(1-y) - a^2(1-y) - b^2 y]."""
    A2, B2 = a*a, b*b
    if Phi == 0 and Psi == 0:
        return -max(A2, B2)
    return _min_scalar(Phi, Psi, A2, B2)


def _min_scalar(Phi, Psi, A2, B2):
    # derivative is monotone increasing (g convex); bisection on g'
    dg = lambda yy: -Phi**2/yy**2 + Psi**2/(1 - yy)**2 + A2 - B2
    eps = 1e-15
    lo_, hi_ = eps, 1 - eps
    if dg(lo_) >= 0:
        yy = lo_ if Phi != 0 else 0.0
    elif dg(hi_) <= 0:
        yy = hi_ if Psi != 0 else 1.0
    else:
        yy = brentq(dg, lo_, hi_, xtol=1e-16, rtol=1e-15)
    val = (Phi**2/yy if Phi != 0 else 0.0) + (Psi**2/(1 - yy) if Psi != 0 else 0.0) - A2*(1 - yy) - B2*yy
    return val


def spherical_orbits(a, b, Phi, Psi, mu=1.0, require_band=True, return_all=False):
    """All spherical photon orbits with given (Phi, Psi): list of dicts sorted by x."""
    xp = horizon(a, b, mu)
    if xp is None:
        raise ValueError("no horizon")
    F = F_coef(a, b, Phi, Psi, mu)
    P = P_coef(a, b, Phi, Psi, mu)
    D = delta_coef(a, b, mu)
    out = []
    for r in np.roots(F):
        if abs(r.imag) > 1e-9*max(1, abs(r.real)) or r.real <= xp*(1 + 1e-12):
            continue
        x0 = r.real
        # polish with Newton on F
        for _ in range(3):
            fx = np.polyval(F, x0); dfx = np.polyval(np.polyder(F), x0)
            if dfx != 0:
                x0 -= fx/dfx
        Dx = np.polyval(D, x0)
        if not (x0 > xp and Dx > 1e-10):
            continue
        K = np.polyval(P, x0)/Dx
        chi = chi_coef(a, b, Phi, Psi, K, mu)
        chi2 = np.polyval(np.polyder(chi, 2), x0)
        band = theta_band(a, b, Phi, Psi, K)
        d = dict(x=x0, r=np.sqrt(x0), K=K, chi2=chi2, unstable=chi2 > 0, band=band)
        if band is not None:
            d['theta_min'] = np.arcsin(np.sqrt(band[0]))
            d['theta_max'] = np.arcsin(np.sqrt(band[1]))
        if (band is not None) or (not require_band):
            out.append(d)
    out.sort(key=lambda d: d['x'])
    return out


# ----------------------------------------------------------------------------- Mino-time constants
def mino_constants(a, b, Phi, Psi, x, mu=1.0):
    """Constants so that, with d sigma = d lambda / rho^2:
       dphi/dsigma = Phi/y - A_phi, dpsi/dsigma = Psi/(1-y) - A_psi,
       dt/dsigma = rho^2 + A_t,  dlambda/dsigma = rho^2 = x + a^2 + (b^2-a^2) y."""
    A2, B2 = a*a, b*b
    D = (x + A2)*(x + B2) - mu*x
    calE = 1 + a*Phi/(x + A2) + b*Psi/(x + B2)
    A_phi = mu*a*(x + B2)*calE/D + (A2 - B2)*Phi/(x + A2)
    A_psi = mu*b*(x + A2)*calE/D - (A2 - B2)*Psi/(x + B2)
    A_t = mu*(x + A2)*(x + B2)*calE/D
    return A_phi, A_psi, A_t, calE


def G_of_y(a, b, K, y1, y2, y):
    """Q(y) = (y - y1)(y2 - y) G(y), G linear and > 0 on the band."""
    A2, B2 = a*a, b*b
    return (B2 - A2)*(y + y1 + y2) - (B2 - 2*A2 - K)


def theta_period_sigma(a, b, K, y1, y2):
    """Mino-time period of the y = sin^2 theta oscillation (complete elliptic integral)."""
    G1, G2 = G_of_y(a, b, K, y1, y2, y1), G_of_y(a, b, K, y1, y2, y2)
    # Lambda = int_0^{2pi} dchi / (2 sqrt(G(chi))), G = G0 - Gc cos chi, G(0)=G1, G(pi)=G2
    m = (G2 - G1)/G2
    return 2*ellipk(m)/np.sqrt(G2)


# ----------------------------------------------------------------------------- full Hamiltonian (BL)
def ginv_numeric(a, b, mu, x, th):
    s2, c2 = np.sin(th)**2, np.cos(th)**2
    rho2 = x + a*a*c2 + b*b*s2
    D = (x + a*a)*(x + b*b) - mu*x
    gtt = ((a*a - b*b)*s2 - (x + a*a)*(D + mu*(x + b*b))/D)/rho2
    gtp = a*mu*(x + b*b)/(rho2*D)
    gts = b*mu*(x + a*a)/(rho2*D)
    gpp = (1/s2 - ((a*a - b*b)*(x + b*b) + b*b*mu)/D)/rho2
    gss = (1/c2 + ((a*a - b*b)*(x + a*a) - a*a*mu)/D)/rho2
    gps = -a*b*mu/(rho2*D)
    gxx = 4*D/rho2
    gthth = 1/rho2
    return gtt, gtp, gts, gpp, gss, gps, gxx, gthth
