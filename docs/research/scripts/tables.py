"""
Unit-test table, frequency/resonance checks, Fig 4.12 four solutions, and the
Sec 3.4 inequality claims.  Run:
uv run --with numpy --with scipy --with sympy --with matplotlib python tables.py
"""
import math
import numpy as np
from reference_impl import solve_spherical, tabulate_period, frequencies, horizon, mino_constants
from figures import solve_constant_theta
from mp5d import spherical_orbits, theta_period_sigma

print("=== Unit-test table (mu = E = 1; theta in degrees; increments per chi-period = one y=sin^2(theta) oscillation)")
print("| a | b | Phi | Psi | x_+ | x = r^2 | r | K | theta_min | theta_max | Dsigma | Dlambda | Dt | Dphi | Dpsi | nu_phi | nu_psi |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
cases = [(0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.7, 0.9), (0.5, 0.0, 0.0, 0.0), (0.5, 0.0, 2.0, 0.0),
         (0.8, 0.0, 0.0, 0.5), (0.0, 0.99, 0.0, 0.0), (0.4, 0.4, 1.2, -0.5), (0.6, 0.3, 1.0, -1.0),
         (0.9, 0.1, 0.8, 0.8), (0.9, 0.09, 0.1, 0.01), (0.55, 0.2, -0.6, 0.4), (0.15, 0.7, 0.9, -0.6),
         (0.5, 0.25, -0.929241, 0.903667), (0.6, 0.3, 0.0, 0.0), (0.6, -0.3, 1.0, 1.0)]
rows = []
for a, b, Phi, Psi in cases:
    o = solve_spherical(a, b, Phi, Psi)
    if not o['ok']:
        print("INVALID", a, b, Phi, Psi, o['reason']); continue
    tab = tabulate_period(a, b, Phi, Psi, o, N=2048)
    fr = frequencies(a, b, Phi, Psi, tab, o)
    th1, th2 = math.degrees(math.asin(math.sqrt(o['y1']))), math.degrees(math.asin(math.sqrt(o['y2'])))
    print(f"| {a} | {b} | {Phi} | {Psi} | {horizon(a,b):.6f} | {o['x']:.6f} | {math.sqrt(o['x']):.6f} | {o['K']:.6f} | {th1:.4f} | {th2:.4f} | "
          f"{tab['dsigma']:.6f} | {tab['dlam']:.6f} | {tab['dt']:.6f} | {tab['dphi']:.6f} | {tab['dpsi']:.6f} | {fr['nu_phi']:.6f} | {fr['nu_psi']:.6f} |")

print("\n=== Allowed-region axis intercepts (closed form, Sec 3.4 / eq 3.13, 3.15) for a few BHs")
print("| a | b | Phi_max (Psi=0) | x | Phi_min (Psi=0) | x | Psi_max (Phi=0) | x | Psi_min (Phi=0) | x |")
print("|---|---|---|---|---|---|---|---|---|---|")
for a, b in [(0, 0), (0.5, 0), (0.6, 0.3), (0.4, 0.4), (0.9, 0.09)]:
    Sp, Sm = math.sqrt((1+a)**2-b*b), math.sqrt((1-a)**2-b*b)
    Tp, Tm = math.sqrt((1+b)**2-a*a), math.sqrt((1-b)**2-a*a)
    print(f"| {a} | {b} | {1+Sp:.6f} | {1+a-b*b+Sp:.6f} | {-1-Sm:.6f} | {1-a-b*b+Sm:.6f} | {1+Tp:.6f} | {1+b-a*a+Tp:.6f} | {-1-Tm:.6f} | {1-b-a*a+Tm:.6f} |")

print("\n=== Frequencies / closure")
# Tangherlini: every orbit closes (great circle): nu_phi = nu_psi = 1/2 per chi-period
for Phi, Psi in [(0.7, 0.9), (-1.2, 0.3), (0.0, 1.5), (1.0, 0.0), (0.0, 0.0)]:
    o = solve_spherical(0, 0, Phi, Psi); tab = tabulate_period(0, 0, Phi, Psi, o, N=1024)
    fr = frequencies(0, 0, Phi, Psi, tab, o)
    print(f"  Tangherlini Phi={Phi} Psi={Psi}: nu_phi={fr['nu_phi']:.10f} nu_psi={fr['nu_psi']:.10f} Dt={tab['dt']:.10f} (2pi={2*math.pi:.10f})")
# continuity of nu_phi (mod 1) across Phi = 0 with the +pi convention
for a, b, Psi in [(0.6, 0.3, 0.5), (0.8, 0.0, 0.5)]:
    vals = []
    for Phi in (-1e-4, 0.0, 1e-4):
        o = solve_spherical(a, b, Phi, Psi); tab = tabulate_period(a, b, Phi, Psi, o, N=4096)
        vals.append(frequencies(a, b, Phi, Psi, tab, o)['nu_phi'])
    print(f"  a={a} b={b} Psi={Psi}: nu_phi at Phi=-1e-4, 0, +1e-4: {vals[0]:.6f}, {vals[1]:.6f}, {vals[2]:.6f}   (equal mod 1)")
# a = b closed forms: Dphi = pi[sgn(Phi) - A/sqrt(a^2+K)], Dpsi = pi[sgn(Psi) - A/sqrt(a^2+K)], Dsigma = pi/sqrt(a^2+K)
for a, Phi, Psi in [(0.4, 1.2, -0.5), (0.3, -0.7, -0.9), (0.45, 0.6, 1.1)]:
    o = solve_spherical(a, a, Phi, Psi); tab = tabulate_period(a, a, Phi, Psi, o, N=4096)
    Aphi, Apsi, At = mino_constants(a, a, Phi, Psi, o['x'])
    sG = math.sqrt(a*a + o['K'])
    print(f"  a=b={a} Phi={Phi} Psi={Psi}: Dphi num {tab['dphi']:.10f} closed {math.pi*(math.copysign(1,Phi) - Aphi/sG):.10f};"
          f" Dpsi num {tab['dpsi']:.10f} closed {math.pi*(math.copysign(1,Psi) - Apsi/sG):.10f}; Dsigma num {tab['dsigma']:.10f} closed {math.pi/sG:.10f}")
# a = b: x, K depend only on Phi + Psi
xs = [solve_spherical(0.4, 0.4, L/2 + d, L/2 - d) for L in (0.7,) for d in (-0.3, 0.0, 0.25)]
print("  a=b=0.4, Phi+Psi=0.7, varying Phi-Psi: x =", [round(q['x'], 12) for q in xs], " K =", [round(q['K'], 12) for q in xs])

print("\n=== Fig 4.12: circular orbits at theta = pi/4 for a=0.5, b=0.25 (all four quadrants)")
sols = []
for gPhi in (-1.5, -1.0, -0.5, 0.5, 1.0, 1.5):
    for gPsi in (-1.5, -1.0, -0.5, 0.5, 1.0, 1.5):
        for gx in (0.8, 1.2, 1.6, 2.0, 2.5, 3.0):
            for gK in (2.0, 3.0, 4.5, 6.0):
                s, ier, res = solve_constant_theta(0.5, 0.25, np.pi/4, [gx, gK, gPhi, gPsi])
                if ier == 1 and res < 1e-10 and s[0] > horizon(0.5, 0.25) and spherical_orbits(0.5, 0.25, s[2], s[3]):
                    if not any(np.allclose(s, q, atol=1e-7) for q in sols):
                        sols.append(s)
for s in sorted(sols, key=lambda q: (np.sign(q[2]), np.sign(q[3]))):
    print(f"   Phi={s[2]:+.6f} Psi={s[3]:+.6f}  x={s[0]:.6f} r={math.sqrt(s[0]):.6f} K={s[1]:.6f}")

print("\n=== Sec 3.4 claim: '+' roots >= x_+ and '-' roots <= x_+ for (a+b)^2 < 1 (grid, a,b >= 0)")
viol = 0; n = 0
for a in np.linspace(0, 0.999, 120):
    for b in np.linspace(0, 0.999, 120):
        if (a + b) >= 1:
            continue
        xp = horizon(a, b)
        for s in (+1, -1):
            S = math.sqrt((1 + s*a)**2 - b*b)
            n += 1
            if 1 + s*a - b*b + S < xp - 1e-12 or 1 + s*a - b*b - S > xp + 1e-12:
                viol += 1
print(f"   checked {n} (a,b,branch) combos, violations: {viol}")
# equality cases
for a, b in [(0.5, 0.5), (0.9, 0.1), (1.0, 0.0)]:
    xp = horizon(a, b)
    print(f"   extremal a={a} b={b}: x_+={xp}, co-rotating root 1-a-b^2+sqrt((1-a)^2-b^2) = {1-a-b*b+math.sqrt(max((1-a)**2-b*b,0))}")
