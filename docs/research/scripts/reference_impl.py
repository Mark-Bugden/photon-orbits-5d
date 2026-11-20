"""
JS-portable reference implementation (plain floats, no numpy in the core) of the
recommended algorithm, plus accuracy tests against a high-accuracy DOP853 solution.

  1. solveSpherical(a, b, Phi, Psi): x from bisection on F = P'Delta - P Delta' on (x_+, inf),
     K = P/Delta, theta band from R(y), validity flag.
  2. tabulatePeriod(...): one chi-period (y-oscillation) with fixed-step RK4, using analytic
     subtraction of the Phi/y and Psi/(1-y) axis terms, so the integrand is smooth.
  3. state(chi): periodic extension  (phi, psi, t, lambda)(chi + 2 pi n) = table(chi) + n*increments.

mu = 1, E = 1 throughout.  Run: uv run --with numpy --with scipy python reference_impl.py
"""
import math


# ------------------------------------------------------------------ polynomials (highest power first)
def polyval(c, x):
    s = 0.0
    for ci in c:
        s = s*x + ci
    return s


def polyder(c):
    n = len(c) - 1
    return [c[i]*(n - i) for i in range(n)]


def horizon(a, b, mu=1.0):
    s = mu - a*a - b*b
    disc = s*s - 4*a*a*b*b
    if disc < 0 and disc > -1e-12:
        disc = 0.0
    if s <= 0 or disc < 0:
        return None
    return 0.5*(s + math.sqrt(disc))


def P_coef(a, b, Phi, Psi, mu=1.0):
    A2, B2 = a*a, b*b
    P1 = (Phi*Phi - Psi*Psi)*(A2 - B2) + A2*B2 + mu*(A2 + B2 + 2*a*Phi + 2*b*Psi)
    P0 = mu*(a*b + b*Phi + a*Psi)**2 + (A2 - B2)*(B2*Phi*Phi - A2*Psi*Psi)
    return [1.0, A2 + B2, P1, P0]


def D_coef(a, b, mu=1.0):
    return [1.0, a*a + b*b - mu, a*a*b*b]


def F_eval(a, b, Phi, Psi, x, mu=1.0):
    P, D = P_coef(a, b, Phi, Psi, mu), D_coef(a, b, mu)
    return polyval(polyder(P), x)*polyval(D, x) - polyval(P, x)*polyval(polyder(D), x)


# ------------------------------------------------------------------ step 1: the spherical orbit
def solve_spherical(a, b, Phi, Psi, mu=1.0):
    xp = horizon(a, b, mu)
    if xp is None:
        return dict(ok=False, reason='no horizon')
    lo = xp*(1 + 1e-12) + 1e-14
    if F_eval(a, b, Phi, Psi, lo, mu) > 0:      # only possible if V is finite at x_+ (calE(x_+)=0)
        return dict(ok=False, reason='F>0 at horizon (degenerate)')
    hi = max(4.0, 2*lo)
    while F_eval(a, b, Phi, Psi, hi, mu) <= 0:
        hi *= 2
    for _ in range(200):                         # bisection: F < 0 left of the minimum of V
        mid = 0.5*(lo + hi)
        if F_eval(a, b, Phi, Psi, mid, mu) > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-15*hi:
            break
    x = 0.5*(lo + hi)
    K = polyval(P_coef(a, b, Phi, Psi, mu), x)/polyval(D_coef(a, b, mu), x)
    band = theta_band(a, b, Phi, Psi, K)
    if band is None:
        return dict(ok=False, reason='theta motion forbidden (outside allowed region)', x=x, K=K)
    return dict(ok=True, x=x, K=K, y1=band[0], y2=band[1])


def theta_band(a, b, Phi, Psi, K, tol=1e-12):
    """Theta >= 0 band in y = sin^2 theta. R = Q / (y if Phi==0) / ((1-y) if Psi==0)."""
    A2, B2 = a*a, b*b
    R = [A2 - B2, B2 - 2*A2 - K, Phi*Phi - Psi*Psi + A2 + K, -Phi*Phi]
    if Phi == 0:
        R = R[:-1]                               # divide by y (constant term is 0)
    if Psi == 0:                                 # divide by (1 - y): synthetic division by (y - 1), sign flip
        q = [R[0]]
        for c in R[1:-1]:
            q.append(c + q[-1])
        R = [-c for c in q]
    cands = ([0.0] if Phi == 0 else []) + ([1.0] if Psi == 0 else [])
    dR = polyder(R)
    # interior critical points of R (dR has degree <= 2)
    while len(dR) and dR[0] == 0.0:
        dR = dR[1:]
    if len(dR) == 3:
        A_, B_, C_ = dR
        disc = B_*B_ - 4*A_*C_
        if disc >= 0:
            sq = math.sqrt(disc)
            q = -0.5*(B_ + math.copysign(sq, B_))
            for r in ([q/A_] + ([C_/q] if q != 0 else [])):
                if 0 < r < 1:
                    cands.append(r)
    elif len(dR) == 2:
        r = -dR[1]/dR[0]
        if 0 < r < 1:
            cands.append(r)
    if not cands:
        return None
    ym = max(cands, key=lambda c: polyval(R, c))
    Rm = polyval(R, ym)
    scale = 1 + abs(K) + Phi*Phi + Psi*Psi
    if Rm < -tol*scale:
        return None
    if Rm <= tol*scale:
        return (ym, ym)

    def bisect(lo, hi, f_lo_neg):
        for _ in range(100):
            m = 0.5*(lo + hi)
            if (polyval(R, m) < 0) == f_lo_neg:
                lo = m
            else:
                hi = m
        return 0.5*(lo + hi)
    y1 = 0.0 if polyval(R, 0.0) >= 0 else bisect(0.0, ym, True)
    y2 = 1.0 if polyval(R, 1.0) >= 0 else bisect(ym, 1.0, False)
    return (y1, y2)


# ------------------------------------------------------------------ step 2: one chi-period table
def mino_constants(a, b, Phi, Psi, x, mu=1.0):
    A2, B2 = a*a, b*b
    D = (x + A2)*(x + B2) - mu*x
    calE = 1 + a*Phi/(x + A2) + b*Psi/(x + B2)
    A_phi = mu*a*(x + B2)*calE/D + (A2 - B2)*Phi/(x + A2)
    A_psi = mu*b*(x + A2)*calE/D - (A2 - B2)*Psi/(x + B2)
    A_t = mu*(x + A2)*(x + B2)*calE/D
    return A_phi, A_psi, A_t


def axis_integral(lo_val, hi_val, chi):
    """int_0^chi dchi' / (A - B cos chi'), with A - B = lo_val > 0, A + B = hi_val: continuous in chi."""
    k = math.sqrt(hi_val/lo_val)
    u = 0.5*chi
    return (2.0/math.sqrt(lo_val*hi_val))*(u + math.atan2((k - 1)*math.sin(u)*math.cos(u),
                                                         math.cos(u)**2 + k*math.sin(u)**2))


def tabulate_period(a, b, Phi, Psi, orb, N=256, subtract=True, mu=1.0):
    """Return arrays over chi in [0, 2pi] (N+1 points): sigma, lam, t, phi, psi and the increments."""
    x, K, y1, y2 = orb['x'], orb['K'], orb['y1'], orb['y2']
    A2, B2 = a*a, b*b
    ybar, dy = 0.5*(y1 + y2), 0.5*(y2 - y1)
    Aphi, Apsi, At = mino_constants(a, b, Phi, Psi, x, mu)
    G = lambda y: (B2 - A2)*(y + y1 + y2) - (B2 - 2*A2 - K)
    G1, G2 = G(y1), G(y2)
    use_phi = subtract and Phi != 0 and y1 > 0
    use_psi = subtract and Psi != 0 and y2 < 1

    def rhs(chi):
        y = ybar - dy*math.cos(chi)
        sG = math.sqrt(G(y))
        ds = 1.0/(2*sG)
        rho2 = x + A2 + (B2 - A2)*y
        if Phi == 0:
            fphi = -Aphi*ds
        elif use_phi:      # Phi/(2 y sqrt G) - Phi/(2 y sqrt G1)  (smooth remainder)
            fphi = Phi/(2*y)*(1/sG - 1/math.sqrt(G1)) - Aphi*ds
        else:
            fphi = (Phi/y - Aphi)*ds
        if Psi == 0:
            fpsi = -Apsi*ds
        elif use_psi:      # Psi/(2(1-y)sqrt G) - Psi/(2(1-y) sqrt G2)
            fpsi = Psi/(2*(1 - y))*(1/sG - 1/math.sqrt(G2)) - Apsi*ds
        else:
            fpsi = (Psi/(1 - y) - Apsi)*ds
        return (ds, rho2*ds, (rho2 + At)*ds, fphi, fpsi)

    h = 2*math.pi/N
    chis = [i*h for i in range(N + 1)]
    Z = [(0.0, 0.0, 0.0, 0.0, 0.0)]
    for i in range(N):        # RHS depends on chi only -> RK4 == Simpson's rule
        c = chis[i]
        k1 = rhs(c); k23 = rhs(c + 0.5*h); k4 = rhs(c + h)
        Z.append(tuple(Z[-1][j] + h/6*(k1[j] + 4*k23[j] + k4[j]) for j in range(5)))
    sig = [z[0] for z in Z]; lam = [z[1] for z in Z]; t = [z[2] for z in Z]
    phi = [z[3] for z in Z]; psi = [z[4] for z in Z]
    if use_phi:   # add analytic part:  Phi/(2 sqrt G1) * int dchi / y,  y = ybar - dy cos chi
        for i, c in enumerate(chis):
            phi[i] += Phi/(2*math.sqrt(G1))*axis_integral(y1, y2, c)
    if use_psi:   # 1 - y = (1-ybar) + dy cos chi = A' - B' cos(chi + pi); integrate from chi = 0
        for i, c in enumerate(chis):
            psi[i] += Psi/(2*math.sqrt(G2))*(axis_integral(1 - y2, 1 - y1, c + math.pi)
                                              - axis_integral(1 - y2, 1 - y1, math.pi))
    return dict(chi=chis, sigma=sig, lam=lam, t=t, phi=phi, psi=psi,
                dsigma=sig[-1], dlam=lam[-1], dt=t[-1], dphi=phi[-1], dpsi=psi[-1])


def frequencies(a, b, Phi, Psi, tab, orb):
    """Coordinate-time frequencies. Phi = 0 (Psi = 0) orbits pass through the axis; in the
    embedding this is a phase jump of pi per chi-period (signed sin/cos theta)."""
    dphi = tab['dphi'] + (math.pi if Phi == 0 else 0.0)
    dpsi = tab['dpsi'] + (math.pi if Psi == 0 else 0.0)
    T = tab['dt']
    return dict(Omega_theta=2*math.pi/T, Omega_phi=dphi/T, Omega_psi=dpsi/T,
                nu_phi=dphi/(2*math.pi), nu_psi=dpsi/(2*math.pi))


# ------------------------------------------------------------------ tests
if __name__ == "__main__":
    import numpy as np
    from figures import orbit_chi
    from mp5d import spherical_orbits
    cases = [(0.0, 0.0, 0.7, 0.9), (0.5, 0.0, 2.0, 0.0), (0.6, 0.3, 1.0, -1.0), (0.9, 0.1, 0.8, 0.8),
             (0.9, 0.09, 0.1, 0.01), (0.4, 0.4, 1.2, -0.5), (0.6, 0.3, 1e-3, 0.5), (0.6, 0.3, -0.8, -1e-3),
             (0.8, 0.0, 0.0, 0.5), (0.3, 0.55, 0.0, 0.0)]
    print("Accuracy of fixed-step tabulation of one chi-period vs DOP853 (rtol 1e-12):")
    print(" case                       |  N  | err(dphi) subtract | err(dphi) naive | err(dpsi) subtract | err(dpsi) naive | err(dt)")
    for a, b, Phi, Psi in cases:
        orb = solve_spherical(a, b, Phi, Psi)
        o = spherical_orbits(a, b, Phi, Psi)[0]
        assert abs(orb['x'] - o['x']) < 1e-9*o['x'], (orb, o)
        ref = orbit_chi(a, b, Phi, Psi, dict(x=orb['x'], K=orb['K'], band=(orb['y1'], orb['y2'])), 2*np.pi, N=3)
        for N in (32, 128, 512):
            ts = tabulate_period(a, b, Phi, Psi, orb, N=N, subtract=True)
            tn = tabulate_period(a, b, Phi, Psi, orb, N=N, subtract=False)
            print(f" ({a},{b},{Phi},{Psi})".ljust(28) + f"| {N:3d} | {abs(ts['dphi']-ref['phi'][-1]):.1e}            | {abs(tn['dphi']-ref['phi'][-1]):.1e}         "
                  f"| {abs(ts['dpsi']-ref['psi'][-1]):.1e}            | {abs(tn['dpsi']-ref['psi'][-1]):.1e}         | {abs(ts['dt']-ref['t'][-1]):.1e}")
