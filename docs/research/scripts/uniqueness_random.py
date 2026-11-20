"""Random test: number of critical points of V(x)=P(x)/Delta(x) outside the horizon
(roots of the quartic F with x > x_+), for random BH spins (incl. near-extremal and
negative a, b) and random (Phi, Psi) in [-6,6]^2 (inside AND outside the allowed region)."""
import numpy as np
from mp5d import *
rng = np.random.default_rng(2024)
counts = {}
worst = []
n = 0
for i in range(200000):
    # sample (a,b) uniformly in the allowed diamond |a|+|b| < 1, biased toward extremal
    u = rng.uniform(0, 1)**0.3 * 0.9999
    w = rng.uniform(0, 1)
    a, b = u*w*rng.choice([-1, 1]), u*(1-w)*rng.choice([-1, 1])
    Phi, Psi = rng.uniform(-6, 6, 2)
    xp = horizon(a, b)
    if xp is None:
        continue
    F = F_coef(a, b, Phi, Psi)
    D = delta_coef(a, b)
    r = np.roots(F)
    good = [z.real for z in r if abs(z.imag) < 1e-7*max(1, abs(z.real)) and z.real > xp and np.polyval(D, z.real) > 1e-10]
    k = len(good)
    counts[k] = counts.get(k, 0) + 1
    n += 1
    if k != 1:
        worst.append((a, b, Phi, Psi, xp, good))
print("samples:", n, " histogram of #critical points of V outside horizon:", counts)
for w_ in worst[:10]:
    print("  ", w_)
