"""Frequency map (Phi,Psi) -> (nu_phi, nu_psi) per theta-oscillation; Jacobian rank; a sample
resonance search (closed orbit) by 2D Newton."""
import math
import numpy as np
from reference_impl import solve_spherical, tabulate_period, frequencies

def nus(a, b, Phi, Psi, N=1024):
    o = solve_spherical(a, b, Phi, Psi)
    if not o['ok']:
        return None
    t = tabulate_period(a, b, Phi, Psi, o, N=N)
    f = frequencies(a, b, Phi, Psi, t, o)
    return np.array([f['nu_phi'], f['nu_psi']])

def jac(a, b, Phi, Psi, h=1e-5):
    J = np.zeros((2, 2))
    for k, (dP, dS) in enumerate([(h, 0), (0, h)]):
        J[:, k] = (nus(a, b, Phi+dP, Psi+dS) - nus(a, b, Phi-dP, Psi-dS))/(2*h)
    return J

for a, b, Phi, Psi in [(0.6, 0.3, 1.0, -1.0), (0.55, 0.2, -0.6, 0.4), (0.4, 0.4, 1.2, -0.5), (0.0, 0.0, 0.7, 0.9), (0.9, 0.1, 0.8, 0.8)]:
    J = jac(a, b, Phi, Psi)
    print(f"a={a} b={b} Phi={Phi} Psi={Psi}: nu={nus(a,b,Phi,Psi).round(6)}, det J = {np.linalg.det(J):.3e}, sv = {np.linalg.svd(J)[1].round(5)}")

# resonance hunt example: a=0.6,b=0.3, target nu = (1/3, -3/5) near Fig 4.7 (nu=(0.3363,-0.6080))
a, b = 0.6, 0.3
target = np.array([1/3, -3/5])
# (try a few targets; report the first that converges inside the region)
p = np.array([1.0, -1.0])
for it in range(20):
    r = nus(a, b, *p, N=2048) - target
    if np.max(np.abs(r)) < 1e-12:
        break
    step = np.linalg.solve(jac(a, b, *p), r)
    lam = 1.0
    while nus(a, b, *(p - lam*step)) is None or np.max(np.abs(nus(a, b, *(p - lam*step)) - target)) > np.max(np.abs(r)):
        lam *= 0.5
        if lam < 1e-6: break
    p = p - lam*step
o = solve_spherical(a, b, *p)
print(f"closed orbit a=0.6 b=0.3 with nu=(1/3,-3/5): Phi={p[0]:.10f} Psi={p[1]:.10f} x={o['x']:.8f} K={o['K']:.8f} (closes after 15 theta-oscillations), residual {np.max(np.abs(r)):.1e}")
# verify closure in R^4 by direct tabulation: after n=15 periods, phi, psi advance by multiples of 2pi
t = tabulate_period(a, b, *p, o, N=4096)
n = 15
print("   n*dphi/2pi =", n*t['dphi']/(2*math.pi), " n*dpsi/2pi =", n*t['dpsi']/(2*math.pi))
