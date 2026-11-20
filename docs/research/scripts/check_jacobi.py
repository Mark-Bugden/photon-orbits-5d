"""Check the closed form y(sigma) = y_e +/- (y2-y1) sn^2(sqrt(G_max) sigma | m), m = 1 - G_min/G_max."""
import numpy as np
from scipy.special import ellipj, ellipk
from mp5d import *
from figures import orbit_chi
for a, b, Phi, Psi in [(0.6, 0.3, 1.0, -1.0), (0.15, 0.7, 0.9, -0.6), (0.9, 0.09, 0.1, 0.01), (0.4, 0.4, 1.2, -0.5)]:
    o = spherical_orbits(a, b, Phi, Psi)[0]
    y1, y2 = o['band']
    G = lambda y: G_of_y(a, b, o['K'], y1, y2, y)
    orb = orbit_chi(a, b, Phi, Psi, o, 4*np.pi, N=4001)
    sig, y = orb['sig'], orb['y']                     # sigma = 0 at chi = 0, i.e. at y = y1 (theta_min)
    G1, G2 = G(y1), G(y2)
    if G1 >= G2:   # expand about y1 (sigma = 0 at theta_min)
        m = 1 - G2/G1
        sn = ellipj(np.sqrt(G1)*sig, m)[0]
        yf = y1 + (y2 - y1)*sn**2
        L = 2*ellipk(m)/np.sqrt(G1)
    else:          # expand about y2: shift sigma origin to theta_max (half period)
        m = 1 - G1/G2
        L = 2*ellipk(m)/np.sqrt(G2)
        sn = ellipj(np.sqrt(G2)*(sig - L/2), m)[0]
        yf = y2 - (y2 - y1)*sn**2
    print(f"a={a} b={b}: m={m:.5f}  max|y_closed - y_numeric| = {np.abs(yf - y).max():.2e};  period closed {L:.10f} vs numeric {sig[-1]/2:.10f}")
