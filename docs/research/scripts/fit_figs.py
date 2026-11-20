"""Digitise Figs 4.9 / 4.11 (theta, phi, psi, 'x' vs lambda) and fit the initial conditions."""
import numpy as np
from PIL import Image
from scipy.optimize import minimize_scalar
from mp5d import *
from integrate_check import mino_orbit

cols = {'x':(0x5E,0x81,0xB5),'theta':(0xE1,0x9C,0x24),'phi':(0x8F,0xB0,0x32),'psi':(0xEB,0x62,0x35)}
figs = {
  'Fig4.9 (cusp, a=0.9,b=0.1,Phi=Psi=0.8)': dict(fn='lambda', a=0.9, b=0.1, Phi=0.8, Psi=0.8, zero_row=609.0, ppu_y=90.0, col0=22.75, ppu_x=183.25, lam0=0.0),
  'Fig4.11 (confused, a=0.9,b=0.09,Phi=0.1,Psi=0.01)': dict(fn='lambda2', a=0.9, b=0.09, Phi=0.1, Psi=0.01, zero_row=512.3, ppu_y=73.75, col0=40-224.67, ppu_x=224.67, lam0=1.0),
}
for name, f in figs.items():
    im = np.asarray(Image.open(f"figs/{f['fn']}_flipped.png").convert('RGB')).astype(int)
    data = {}
    for k, c in cols.items():
        mask = np.abs(im - np.array(c)).sum(axis=2) < 40
        L = []; V = []
        for col in range(45 if f['fn']=='lambda2' else 25, 935):
            rows = np.where(mask[:, col])[0]
            if len(rows) and rows.max() - rows.min() < 12:      # unambiguous single crossing
                L.append((col - f['col0'])/f['ppu_x'] if f['fn']=='lambda' else (col - 40)/f['ppu_x'] + 1.0)
                V.append((f['zero_row'] - rows.mean())/f['ppu_y'])
        data[k] = (np.array(L), np.array(V))
    o = spherical_orbits(f['a'], f['b'], f['Phi'], f['Psi'])[0]
    print(f"== {name}\n   theory: x = {o['x']:.5f}, r = sqrt(x) = {o['r']:.5f}, K = {o['K']:.5f}, theta in [{o['theta_min']:.4f}, {o['theta_max']:.4f}] rad")
    print(f"   digitised 'x' curve: mean {data['x'][1].mean():.4f} (std {data['x'][1].std():.4f})  -> the curve labelled x is r = sqrt(x)")
    chis, Z = mino_orbit(f['a'], f['b'], f['Phi'], f['Psi'], o, 6, N=60000)
    sig, lam, tt, ph, ps = Z
    y1, y2 = o['band']
    th = np.arcsin(np.sqrt(np.clip(0.5*(y1+y2) - 0.5*(y2-y1)*np.cos(chis), 0, 1)))
    Lt, Vt = data['theta']
    def model_theta(i0):
        lam_rel = lam - lam[i0] + f['lam0']
        return np.interp(Lt, lam_rel, th)
    best = None
    per = len(chis)//6
    for i0 in range(0, 2*per, 20):
        err = np.mean((model_theta(i0) - Vt)**2)
        if best is None or err < best[0]: best = (err, i0)
    i0 = best[1]
    for i in range(max(0, i0-40), i0+40):
        err = np.mean((model_theta(i) - Vt)**2)
        if err < best[0]: best = (err, i)
    i0 = best[1]
    lam_rel = lam - lam[i0] + f['lam0']
    th0 = th[i0]; chi0 = chis[i0] % (2*np.pi)
    dirn = 'increasing' if (chi0 % (2*np.pi)) < np.pi else 'decreasing'
    Lp, Vp = data['phi']; Ls, Vs = data['psi']
    phi_model = np.interp(Lp, lam_rel, ph - ph[i0]); psi_model = np.interp(Ls, lam_rel, ps - ps[i0])
    phi0 = np.median(Vp - phi_model); psi0 = np.median(Vs - psi_model)
    rms_phi = np.sqrt(np.mean((Vp - phi_model - phi0)**2)); rms_psi = np.sqrt(np.mean((Vs - psi_model - psi0)**2))
    print(f"   fitted start (lambda = {f['lam0']}): theta0 = {th0:.4f} rad ({dirn}), phi0 = {phi0:.4f}, psi0 = {psi0:.4f}")
    print(f"   rms residuals (rad): theta {np.sqrt(best[0]):.4f}, phi {rms_phi:.4f}, psi {rms_psi:.4f}")
    # report values at a few lambdas
    for Lq in [f['lam0'], 2.0, 3.0, 4.0, 5.0]:
        print(f"     lambda={Lq:.1f}: theta={np.interp(Lq, lam_rel, th):.3f} phi={np.interp(Lq, lam_rel, ph-ph[i0])+phi0:.3f} psi={np.interp(Lq, lam_rel, ps-ps[i0])+psi0:.3f}")

# ---- extra: values at lambda = 0 for Fig 4.11 and 3D renderings with fitted initial data
if __name__ == "__main__":
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.colors import hsv_to_rgb
    fig = plt.figure(figsize=(12, 6))
    specs = [('4.8 (lambda 0..16)', 0.9, 0.1, 0.8, 0.8, 0.0, 0.3924, +1, 1.7476, 0.4490, 0.0, 16.0),
             ('4.10 (lambda 1..5)', 0.9, 0.09, 0.1, 0.01, 1.0, 1.0214, -1, -0.9799, 2.9571, 1.0, 5.0)]
    for k, (nm, a, b, Phi, Psi, lam_s, th_s, sgn, ph_s, ps_s, l0, l1) in enumerate(specs):
        o = spherical_orbits(a, b, Phi, Psi)[0]
        y1, y2 = o['band']
        ys = np.sin(th_s)**2
        c = np.arccos(np.clip((0.5*(y1+y2) - ys)/(0.5*(y2-y1)), -1, 1))   # chi with y(chi)=ys, chi in [0,pi]
        chi0 = c if sgn > 0 else 2*np.pi - c
        from figures import orbit_chi
        orb = orbit_chi(a, b, Phi, Psi, o, 8*np.pi, N=40000, chi0=chi0 - 4*np.pi)
        i0 = np.argmin(np.abs(orb['chi'] - chi0))
        lam = orb['lam'] - orb['lam'][i0] + lam_s
        ph = orb['phi'] - orb['phi'][i0] + ph_s; ps = orb['psi'] - orb['psi'][i0] + ps_s
        th = np.arcsin(np.sqrt(np.clip(orb['y'], 0, 1)))
        j = np.argmin(np.abs(lam - 0.0))
        print(f"Fig {nm}: extrapolated values at lambda=0: theta={th[j]:.4f}, phi={ph[j]:.4f}, psi={ps[j]:.4f}")
        sel = (lam >= l0) & (lam <= l1)
        X = np.sqrt(o['x']+a*a)*np.sqrt(orb['y'][sel])*np.cos(ph[sel]); Y = np.sqrt(o['x']+a*a)*np.sqrt(orb['y'][sel])*np.sin(ph[sel])
        Z = np.sqrt(o['x']+b*b)*np.sqrt(1-orb['y'][sel])
        hue = (ps[sel] % (2*np.pi))/(2*np.pi)
        ax = fig.add_subplot(1, 2, k+1, projection='3d')
        ax.scatter(X, Y, Z, c=hsv_to_rgb(np.c_[hue, np.ones_like(hue), np.ones_like(hue)]), s=1)
        u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
        ax.plot_wireframe(np.sqrt(o['x']+a*a)*np.cos(u)*np.sin(v), np.sqrt(o['x']+a*a)*np.sin(u)*np.sin(v), np.sqrt(o['x']+b*b)*np.cos(v), color='0.85', lw=0.3)
        ax.view_init(elev=25, azim=60 if k == 0 else 30)
        ax.set_title('Fig ' + nm + ' with fitted initial data')
    plt.tight_layout(); plt.savefig('figs_fitted_3d.png', dpi=100)
