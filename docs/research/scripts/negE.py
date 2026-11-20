import numpy as np
from mp5d import *
for a,b in [(0.6,0.3),(0.9,0.09),(0.5,0.25),(0.1,0.85),(0.7,0.2999)]:
    xp=horizon(a,b)
    grid=np.linspace(-4.6,4.6,185)
    bad=[]
    for Phi in grid:
        for Psi in grid:
            orbs=[o for o in spherical_orbits(a,b,Phi,Psi) ]
            for o in orbs:
                Aphi,Apsi,At,calE=mino_constants(a,b,Phi,Psi,o['x'])
                d=[o['x']+a*a+(b*b-a*a)*yy+At for yy in o['band']]
                # direct check using inverse metric: tdot = -g^tt*1 + g^tphi Phi + g^tpsi Psi at mid-band
                ym=0.5*(o['band'][0]+o['band'][1]); th=np.arcsin(np.sqrt(ym))
                gtt,gtp,gts,gpp,gss,gps,gxx,gthth=ginv_numeric(a,b,1.0,o['x'],th)
                tdot=-gtt+gtp*Phi+gts*Psi
                if min(d)<0 or tdot<0:
                    bad.append((Phi,Psi,o['x'],calE,d,tdot))
    print(a,b,'x+=',xp,'count negative-tdot orbits:',len(bad))
    if bad:
        bad.sort(key=lambda t: t[0]**2+t[1]**2)
        for row in bad[:3]: print('   Phi=%.3f Psi=%.3f x=%.4f calE=%.4f dt/dsig ends=(%.3f,%.3f) tdot_mid=%.4f'%(row[0],row[1],row[2],row[3],row[4][0],row[4][1],row[5]))
        P=np.array([[r[0],r[1]] for r in bad]); print('   Phi range',P[:,0].min(),P[:,0].max(),' Psi range',P[:,1].min(),P[:,1].max())
        xs=np.array([r[2] for r in bad]); print('   x range',xs.min(),xs.max(), ' ergo x_stat range', 1-max(a,b)**2, 1-min(a,b)**2)
