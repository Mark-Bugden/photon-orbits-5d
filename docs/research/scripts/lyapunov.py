"""How long does a full (x, theta, p_x, p_theta) Hamiltonian integration (the paper's
method, Sec 3.3) stay on the unstable spherical orbit?  Stop when |x-x0| > 1e-3 x0."""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from mp5d import *
import integrate_check as ic
for (a,b,Phi,Psi) in [(0.5,0.0,2.0,0.0),(0.6,0.3,1.0,-1.0),(0.0,0.0,0.7,0.9),(0.9,0.09,0.1,0.01)]:
    o=spherical_orbits(a,b,Phi,Psi)[0]
    y1,y2=o['band']; Lsig=theta_period_sigma(a,b,o['K'],y1,y2)
    gsig=np.sqrt(2*o['chi2'])   # radial e-folding rate in Mino time (near-orbit: d(dx)/dsigma = sqrt(2 chi'') dx)
    th0=np.arcsin(np.sqrt(y1)) if y1>0 else 1e-9
    ic.rhs_ham.Phi, ic.rhs_ham.Psi = Phi, Psi
    for rtol in (1e-8, 1e-12):
        ev=lambda lam,Y,*args: abs(Y[1]-o['x'])-1e-3*o['x']; ev.terminal=True
        sol=solve_ivp(ic.rhs_ham,(0,400),[0,o['x'],th0,0,0,0,0],args=(a,b,1.0),method='DOP853',rtol=rtol,atol=rtol*1e-1,events=ev)
        lam_end=sol.t[-1]
        # number of theta periods (lambda-period from Mino table)
        from figures import orbit_chi
        one=orbit_chi(a,b,Phi,Psi,o,2*np.pi,N=201)
        print(f"a={a} b={b} Phi={Phi} Psi={Psi}: e-folds per theta-period ~{gsig*Lsig:.2f}; rtol={rtol:g}: |dx|>1e-3 x0 after lambda={lam_end:.1f} = {lam_end/one['lam'][-1]:.1f} theta-periods", flush=True)
