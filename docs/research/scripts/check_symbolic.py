"""
Symbolic checks of Bugden, arXiv:1801.03389 (5D Myers-Perry spherical photon orbits).

Conventions: coordinates (t, x, theta, phi, psi), x = r^2, metric as in eq (2.1)
(with the (x+b)^2 typo corrected to (x+b^2)). p_t = -E, p_phi = Phi, p_psi = Psi.
To keep sympy fast, the angular dependence is written with S = sin^2(theta),
so sin^2 -> S, cos^2 -> 1-S, and all checks are rational-function identities
(zero test = cancel(numerator) == 0).
Run:  uv run --with sympy python check_symbolic.py
"""
import sympy as sp

x, S, y = sp.symbols('x S y', real=True)
mu, a, b = sp.symbols('mu a b', real=True)
E, Phi, Psi, K = sp.symbols('E Phi Psi K', real=True)
px, pth = sp.symbols('p_x p_theta', real=True)

s2, c2 = S, 1 - S
rho2 = x + a**2*c2 + b**2*s2
Delta = (x + a**2)*(x + b**2) - mu*x


import random as _rnd
_rnd.seed(12345)


def zero(e, trials=8):
    """Zero test for a rational expression: exact evaluation at random rational points
    (a nonzero rational function vanishes at random points with probability 0)."""
    fs = sorted(e.free_symbols, key=str)
    for _ in range(trials):
        sub = {v: sp.Rational(_rnd.randint(1, 997), _rnd.randint(1, 97)) * (1 if v != S else sp.Rational(1, 1)) for v in fs}
        if S in sub:
            sub[S] = sp.Rational(_rnd.randint(1, 98), 99)   # 0 < S < 1
        val = sp.nsimplify(e.subs(sub)) if False else sp.simplify(e.subs(sub))
        if val != 0:
            return False
    return True


def rz(e):
    return 'OK (exactly 0 at 8 random rational points)' if zero(e) else 'NONZERO'


# ---------------------------------------------------------------- metric
# (t, phi, psi) block of g_{mu nu}
kv = [1, a*s2, b*c2]
G = sp.Matrix(3, 3, lambda i, j: mu/rho2*kv[i]*kv[j])
G[0, 0] += -1
G[1, 1] += (x + a**2)*s2
G[2, 2] += (x + b**2)*c2
Ginv = G.inv(method='LU').applyfunc(lambda e: sp.factor(sp.together(e)))

paper = {
    (0, 0): 1/rho2*((a**2 - b**2)*s2 - (x + a**2)*(Delta + mu*(x + b**2))/Delta),
    (0, 1): a*mu*(x + b**2)/(rho2*Delta),
    (0, 2): b*mu*(x + a**2)/(rho2*Delta),
    (1, 1): 1/rho2*(1/s2 - ((a**2 - b**2)*(x + b**2) + b**2*mu)/Delta),
    (2, 2): 1/rho2*(1/c2 + ((a**2 - b**2)*(x + a**2) - a**2*mu)/Delta),
    (1, 2): -a*b*mu/(rho2*Delta),
}
names = {(0, 0): 'g^tt (3.3c)', (0, 1): 'g^tphi (3.3d)', (0, 2): 'g^tpsi (3.3e)',
         (1, 1): 'g^phiphi (3.3f)', (2, 2): 'g^psipsi (3.3g)', (1, 2): 'g^phipsi (3.3h)'}
print("=== Inverse metric (3.3) ===")
for k_, v in paper.items():
    print(f"{names[k_]}: computed - paper -> {rz(Ginv[k_] - v)}")
print("g^xx = 4 Delta/rho^2 and g^thth = 1/rho^2 follow from g_xx = rho^2/(4Delta), g_thth = rho^2 (diagonal).")

# ---------------------------------------------------------------- HJ separation
calE = E + a*Phi/(x + a**2) + b*Psi/(x + b**2)
calX = Delta*(x*E**2 + (a**2 - b**2)*(Phi**2/(x + a**2) - Psi**2/(x + b**2)) - K) \
    + mu*(x + a**2)*(x + b**2)*calE**2
Theta = E**2*(a**2*c2 + b**2*s2) + K - Phi**2/s2 - Psi**2/c2

pv = sp.Matrix([-E, Phi, Psi])
twoH_ang = (pv.T*Ginv*pv)[0, 0]
twoH = twoH_ang + 4*Delta/rho2*px**2 + pth**2/rho2
print("\n=== Hamilton-Jacobi separation ===")
print("rho^2*2H - [4 Delta p_x^2 + p_th^2 - (X + Delta*Theta)/Delta]:",
      rz(rho2*twoH - (4*Delta*px**2 + pth**2 - (calX + Delta*Theta)/Delta)))
print("f = -(X+Delta Theta)/(2 rho^2 Delta) independent of K:",
      rz(sp.diff(-(calX + Delta*Theta)/(2*rho2*Delta), K)))
H37 = 2*Delta/rho2*px**2 + pth**2/(2*rho2) - (calX + Delta*Theta)/(2*rho2*Delta)
print("(3.7) H - (1/2) g^{mn} p_m p_n:", rz(H37 - twoH/2))

# velocities
qd = Ginv*pv   # (tdot, phidot, psidot)
eq_t = E*rho2 + mu*(x + a**2)*(x + b**2)*calE/Delta
eq_phi = Phi/s2 - mu*a*(x + b**2)*calE/Delta - (a**2 - b**2)*Phi/(x + a**2)
eq_psi = Psi/c2 - mu*b*(x + a**2)*calE/Delta + (a**2 - b**2)*Psi/(x + b**2)
print("\n=== First-order equations (2.3) ===")
print("(2.3a) rho^2 tdot:", rz(rho2*qd[0] - eq_t))
print("(2.3d) rho^2 phidot:", rz(rho2*qd[1] - eq_phi))
print("(2.3e) rho^2 psidot:", rz(rho2*qd[2] - eq_psi))
print("(3.8c) tdot = (1/(2 rho^2 Delta)) d(X+Delta Theta)/dE:",
      rz(qd[0] - sp.diff(calX + Delta*Theta, E)/(2*rho2*Delta)))
print("(3.8d) phidot = -(1/(2 rho^2 Delta)) d(X+Delta Theta)/dPhi:",
      rz(qd[1] + sp.diff(calX + Delta*Theta, Phi)/(2*rho2*Delta)))
print("(3.8e) psidot:", rz(qd[2] + sp.diff(calX + Delta*Theta, Psi)/(2*rho2*Delta)))
print("(2.3b)/(3.5a): rho^2 xdot = 4 Delta p_x, so (rho^2 xdot)^2 = 16 Delta^2 p_x^2 = 4X  <=> p_x^2 = X/(4 Delta^2)  [consistent]")

# ---------------------------------------------------------------- cubic chi
print("\n=== Cubic chi (Sec 3.2), E = 1 ===")
chi = sp.cancel(sp.together(calX.subs(E, 1)))
chi_poly = sp.Poly(chi, x)
print("chi polynomial in x, degree", chi_poly.degree())
a3, a2, a1, a0 = [sp.expand(cf) for cf in chi_poly.all_coeffs()]
pa = [1, a**2 + b**2 - K,
      (Phi**2 - Psi**2)*(a**2 - b**2) + K*(mu - a**2 - b**2) + a**2*b**2 + mu*(a**2 + b**2 + 2*a*Phi + 2*b*Psi),
      Phi**2*b**2*(a**2 - b**2 + mu) + Psi**2*a**2*(b**2 - a**2 + mu) - K*a**2*b**2
      + a*b*mu*(a*b + 2*Phi*b + 2*Psi*a + 2*Phi*Psi)]
for nm, m, pp in zip(['a3', 'a2', 'a1', 'a0'], [a3, a2, a1, a0], pa):
    print(f"{nm}: computed - paper = {sp.expand(m - pp)}")
chiE = sp.Poly(sp.cancel(sp.together(calX)), x)
print("general-E coefficients (homogeneous of degree 2 in (E,Phi,Psi,sqrt K)):")
for cf in chiE.all_coeffs():
    print("    ", sp.expand(cf))

P = sp.expand(chi + K*Delta)
print("P := chi + K*Delta is K-free:", sp.diff(P, K) == 0)
Pp = sp.Poly(P, x)
print("P(x) coefficients:", [sp.expand(cf) for cf in Pp.all_coeffs()])
P0_alt = mu*(a*b + b*Phi + a*Psi)**2 + (a**2 - b**2)*(b**2*Phi**2 - a**2*Psi**2)
print("P0 = mu(ab + b Phi + a Psi)^2 + (a^2-b^2)(b^2 Phi^2 - a^2 Psi^2):",
      sp.expand(Pp.all_coeffs()[-1] - P0_alt) == 0)
# compact form of P
Pform = x*(x + a**2)*(x + b**2) + mu*(x*(1) + (a + b*0)*0)  # placeholder (not used)
T = (x + a**2)*(x + b**2) + a*(x + b**2)*Phi + b*(x + a**2)*Psi   # = (x+a^2)(x+b^2) calE
Pcand = sp.cancel((Delta*(x + (a**2 - b**2)*(Phi**2/(x + a**2) - Psi**2/(x + b**2)))
                   + mu*T**2/((x + a**2)*(x + b**2))))
print("P equals Delta*[x + (a^2-b^2)(Phi^2/(x+a^2) - Psi^2/(x+b^2))] + mu*T^2/((x+a^2)(x+b^2)):",
      zero(Pcand - P))
# try: P = [x (x+a^2)(x+b^2) - ...]?  Express P in terms of alpha=x+a^2, beta=x+b^2
F = sp.expand(sp.diff(P, x)*Delta - P*sp.diff(Delta, x))
Fp = sp.Poly(F, x)
print("F(x) = P' Delta - P Delta' degree:", Fp.degree(), " leading coeff:", Fp.all_coeffs()[0])
print("F coefficients:")
for k_, cf in enumerate(Fp.all_coeffs()):
    print(f"   x^{4-k_}:", sp.factor(cf))

# ---------------------------------------------------------------- horizon / Omega
print("\n=== Horizon & angular velocity sign ===")
xp = (mu - a**2 - b**2 + sp.sqrt((mu - a**2 - b**2)**2 - 4*a**2*b**2))/2
print("Delta(x_+) = 0:", sp.simplify(sp.expand(Delta.subs(x, xp))) == 0)
num = {mu: 1, a: sp.Rational(3, 10), b: sp.Rational(1, 5)}
xpn = sp.nsimplify(xp.subs(num))
for sgn in (+1, -1):
    Oa = sgn*a/(xp + a**2); Ob = sgn*b/(xp + b**2)
    xi = sp.Matrix([1, Oa, Ob])
    nrm = (xi.T*G*xi)[0, 0]
    vals = [sp.N(nrm.subs(num).subs(x, xpn).subs(S, Sv), 12) for Sv in (sp.Rational(1, 10), sp.Rational(1, 2), sp.Rational(9, 10))]
    print(f"xi = d_t {'+' if sgn>0 else '-'} [a/(x_+ + a^2)] d_phi {'+' if sgn>0 else '-'} [b/(x_+ + b^2)] d_psi : g(xi,xi)|_H =", vals)
print("frame dragging omega_phi = -g_tphi/g_phiphi (theta=pi/2, b=0):",
      sp.factor(sp.together((-G[0, 1]/G[1, 1]).subs({S: 1, b: 0}))))

# ---------------------------------------------------------------- flat limit embedding
print("\n=== mu = 0: BL coords are flat spheroidal coordinates of R^4 ===")
th, ph, ps = sp.symbols('theta phi psi', real=True)
X1 = sp.sqrt(x + a**2)*sp.sin(th)*sp.cos(ph); Y1 = sp.sqrt(x + a**2)*sp.sin(th)*sp.sin(ph)
X2 = sp.sqrt(x + b**2)*sp.cos(th)*sp.cos(ps); Y2 = sp.sqrt(x + b**2)*sp.cos(th)*sp.sin(ps)
qs = [x, th, ph, ps]
J = sp.Matrix([[sp.diff(Fn, q) for q in qs] for Fn in (X1, Y1, X2, Y2)])
gfl = (J.T*J).applyfunc(lambda e: sp.trigsimp(sp.simplify(e)))
rho2t = x + a**2*sp.cos(th)**2 + b**2*sp.sin(th)**2
target = sp.diag(rho2t/(4*(x + a**2)*(x + b**2)), rho2t, (x + a**2)*sp.sin(th)**2, (x + b**2)*sp.cos(th)**2)
print("pullback of dX^2 - spatial BL metric at mu=0:", (gfl - target).applyfunc(lambda e: sp.simplify(e)))

# ---------------------------------------------------------------- symmetry
print("\n=== Symmetry (2.2): a<->b, theta->pi/2-theta (S->1-S), phi<->psi (Phi<->Psi) ===")
sw = {a: b, b: a, Phi: Psi, Psi: Phi, S: 1 - S}
print("calX invariant:", rz(calX.subs(sw, simultaneous=True) - calX))
print("Theta invariant:", rz(Theta.subs(sw, simultaneous=True) - Theta))
print("Delta, rho^2 invariant:", rz(Delta.subs(sw, simultaneous=True) - Delta), rz(rho2.subs(sw, simultaneous=True) - rho2))

# ---------------------------------------------------------------- Theta(y)
print("\n=== Theta in y = sin^2 theta (Sec 3.1.2), E = 1 ===")
Th_y = Theta.subs(E, 1).subs(S, y)
Qpaper = (a**2 - b**2)*y**3 + (b**2 - 2*a**2 - K)*y**2 + (Phi**2 - Psi**2 + a**2 + K)*y - Phi**2
print("y(1-y)Theta - Q_paper:", rz(y*(1 - y)*Th_y - Qpaper))
F0 = [Qpaper.subs(y, 0)] + [sp.diff(Qpaper, y, k_).subs(y, 0) for k_ in (1, 2, 3)]
F1 = [Qpaper.subs(y, 1)] + [sp.diff(Qpaper, y, k_).subs(y, 1) for k_ in (1, 2, 3)]
print("F_Q(0) =", [sp.expand(f_) for f_ in F0])
print("F_Q(1) =", [sp.expand(f_) for f_ in F1])
print("dy/dlambda = 2 sin cos thetadot => (rho^2 dy/dlambda)^2 = 4 S(1-S) Theta = 4 Q(y=S):",
      rz(4*S*(1 - S)*Theta.subs(E, 1) - 4*Qpaper.subs(y, S)))
# factorisation Q = (y-y1)(y2-y) G(y), G linear
y1, y2 = sp.symbols('y1 y2')
Gc = (b**2 - a**2)*(y + y1 + y2) - (b**2 - 2*a**2 - K)
print("If y1,y2 are roots of Q: Q - (y-y1)(y2-y)G(y) with G=(b^2-a^2)(y+y1+y2)-(b^2-2a^2-K) is")
rem = sp.expand(Qpaper - (y - y1)*(y2 - y)*Gc)
# reduce using Q(y1)=Q(y2)=0: remainder must be linear combination; check numerically
import random
ok = True
for trial in range(5):
    vals = {a: sp.Rational(random.randint(1, 9), 10), b: sp.Rational(random.randint(0, 9), 10),
            K: sp.Rational(random.randint(20, 60), 10), Phi: sp.Rational(random.randint(-9, 9), 10),
            Psi: sp.Rational(random.randint(-9, 9), 10)}
    Qn = sp.Poly(Qpaper.subs(vals), y)
    rts = [complex(r_) for r_ in Qn.nroots(n=30)]
    r1, r2 = rts[0], rts[1]
    diff = sp.Poly(sp.expand(rem.subs(vals).subs({y1: r1, y2: r2})), y).all_coeffs()
    ok &= all(abs(complex(d_)) < 1e-12 for d_ in diff)
print("   Q = (y-y1)(y2-y) G(y) for any two roots (random numeric tests):", ok)

# a = b case
print("\n=== a = b: max over theta of -Phi^2/sin^2 - Psi^2/cos^2 ===")
P_, Q_ = sp.symbols('P Q', positive=True)
ys = P_/(P_ + Q_)
print("at sin^2 = |Phi|/(|Phi|+|Psi|): value =", sp.factor(sp.simplify(-P_**2/ys - Q_**2/(1 - ys))),
      "  [P=|Phi|, Q=|Psi|]  => Theta_max = a^2 + K - (|Phi|+|Psi|)^2")

# ---------------------------------------------------------------- equatorial orbits Sec 3.4
print("\n=== Equatorial circular orbits (Sec 3.4): mu=1, Psi=0, K=Phi^2-b^2 ===")
chi_eq = sp.expand(chi.subs({mu: 1, Psi: 0, K: Phi**2 - b**2}))
ce = sp.Poly(chi_eq, x).all_coeffs()
paper_c1 = a**2 + Phi**2 + 2*a*Phi + b**4 + 2*a**2*b**2 - 2*Phi**2*b**2
print("x^2 coeff:", sp.expand(ce[1]), "  (paper: a^2 + 2b^2 - Phi^2)")
print("x^1 coeff - paper:", sp.expand(ce[2] - paper_c1))
print("x^0 coeff - paper's printed constant:", sp.factor(sp.expand(ce[3] - paper_c1)))
print("correct x^0 coeff:", sp.factor(ce[3]), " = ", sp.expand(ce[3]))
for sa in (+1, -1):
    for sr in (+1, -1):
        Sq = sp.sqrt((1 + sa*a)**2 - b**2)
        xs = 1 + sa*a - b**2 + sr*Sq
        Phis = sa*(1 + sr*Sq)
        v1 = sp.simplify(sp.expand(chi_eq.subs({x: xs, Phi: Phis})))
        v2 = sp.simplify(sp.expand(sp.diff(chi_eq, x).subs({x: xs, Phi: Phis})))
        print(f"  (3.9)-(3.12) branch (1{'+' if sa>0 else '-'}a), sqrt sign {sr:+d}: chi = {v1}, chi' = {v2}")
# Check with the paper's (incorrect) constant: do (3.9)-(3.10) still solve it?
chi_eq_paper = x**3 + ce[1]*x**2 + paper_c1*x + paper_c1
Sq = sp.sqrt((1 + a)**2 - b**2)
v = sp.simplify(chi_eq_paper.subs({x: 1 + a - b**2 + Sq, Phi: 1 + Sq}).subs({a: sp.Rational(1, 2), b: sp.Rational(1, 4)}))
print("  paper's printed chi at its own solution (a=1/2,b=1/4):", sp.N(v), " (nonzero => printed constant is a typo, solutions are right)")
sols = sp.solve([chi_eq, sp.diff(chi_eq, x)], [x, Phi], dict=True)
print("number of (x,Phi) solutions of chi=chi'=0:", len(sols))
for so in sols:
    print("    x =", sp.simplify(so[x]), "   Phi =", sp.simplify(so[Phi]))
