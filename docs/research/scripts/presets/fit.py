"""Fit each paper figure's drawn segment, viewpoint and colour offset, with the reference model's path.

The picture is our path (no arctan(a/r) twist), seen from azimuth α and elevation ε with z up, by an orthographic
camera scaled so the spheroid's outline matches the figure's. Turning the path about z is the same as turning the
camera, so α absorbs the start's φ (and the paper's twist). Unknowns: the segment [χ_t, χ_t + 2πN], α, ε, ψ₀.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from measure import NAMES, hsv, load
from photon_orbits import between, embed, node, solve_orbit, tabulate

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * math.pi
PARAMS = {
    "4.2": (0.5, 0.0, 0.0, 0.0), "4.3": (0.5, 0.0, 2.0, 0.0), "4.4": (0.8, 0.0, 0.0, 0.5), "4.5": (0.8, 0.0, 0.0, 0.5),
    "4.6": (0.0, 0.99, 0.0, 0.0), "4.7": (0.6, 0.3, 1.0, -1.0), "4.8": (0.9, 0.1, 0.8, 0.8),
    "4.10": (0.9, 0.09, 0.1, 0.01), "4.12": (0.5, 0.25, -0.929241, 0.903667),
}
PER = 720  # samples per oscillation


class Path3:
    def __init__(self, fig, osc):
        self.orbit = solve_orbit(*PARAMS[fig])
        assert getattr(self.orbit, "x", None), fig
        self.table = tabulate(self.orbit)
        n = self.table.n
        rows = []
        for k in range(int(osc * PER) + 1):
            u = k * n / PER
            i = math.floor(u)
            t = u - i
            rows.append(node(self.table, i) if t == 0 else between(self.orbit, self.table, i, t))
        rows = np.array(rows)
        self.chi, self.lam, self.phi, self.psi = rows[:, 0], rows[:, 2], rows[:, 3], rows[:, 4]
        self.xyz = np.array([embed(self.orbit, c, p) for c, p in zip(self.chi, self.phi)])
        o = self.orbit
        self.A = math.sqrt(o.x + o.a * o.a)
        self.B = math.sqrt(o.x + o.b * o.b)

    def segment(self, chi_t, n_osc):
        k0 = int(round(chi_t / TWO_PI * PER))
        k1 = int(round((chi_t / TWO_PI + n_osc) * PER))
        return slice(k0, k1 + 1)


def camera(az, el):
    a, e = math.radians(az), math.radians(el)
    R = np.array([-math.sin(a), math.cos(a), 0])
    U = np.array([-math.sin(e) * math.cos(a), -math.sin(e) * math.sin(a), math.cos(e)])
    V = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
    return R, U, V


VIEW_DISTANCE = 3.38  # Mathematica's default |ViewPoint|, in lengths of the box's longest side


def project(path, sl, az, el, m):
    """Perspective, as Mathematica draws: the camera VIEW_DISTANCE box lengths from the centre, z up."""
    R, U, V = camera(az, el)
    P = path.xyz[sl]
    D = VIEW_DISTANCE * 2 * max(path.A, path.B)
    C = D * V
    normal = P / np.array([path.A ** 2, path.A ** 2, path.B ** 2])
    front = np.einsum("ij,ij->i", normal, C - P) > 0
    depth = D - P @ V
    f = m["rx"] / math.tan(math.asin(path.A / D))  # the outline's half-width, about, for the equatorial radius
    px = m["cx"] + f * (P @ R) / depth
    py = m["cy"] - f * (P @ U) / depth
    return px, py, front


class Target:
    def __init__(self, m):
        self.m = m
        w, h = m["size"]
        p = m["path"]
        self.pts = np.c_[p["x"], p["y"]].astype(float)
        self.hue = np.array(p["hue"])
        mask = np.zeros((h, w), bool)
        mask[p["y"], p["x"]] = True
        self.dt = ndimage.distance_transform_edt(~mask)
        self.w, self.h = w, h
        rng = np.random.default_rng(0)
        pick = rng.choice(len(self.pts), min(2000, len(self.pts)), replace=False)
        self.sub = self.pts[pick]
        self.sub_hue = self.hue[pick] * TWO_PI

    def score(self, px, py, front, psi=None):
        """Mean distance both ways, in pixels, capped at 12; plus the spread of the colour offset, if psi given."""
        x, y = px[front], py[front]
        if len(x) < 10:
            return 99.0
        xi = np.clip(np.round(x).astype(int), 0, self.w - 1)
        yi = np.clip(np.round(y).astype(int), 0, self.h - 1)
        d1 = np.minimum(self.dt[yi, xi], 12).mean()
        step = max(1, len(x) // 4000)
        tree = cKDTree(np.c_[x[::step], y[::step]])
        dist, idx = tree.query(self.sub)
        d2 = np.minimum(dist, 12).mean()
        if psi is None:
            return d1 + d2
        near = dist < 4
        if near.sum() < 20:
            return d1 + d2 + 4
        delta = self.sub_hue[near] - psi[front][::step][idx[near]]
        length = math.hypot(np.sin(delta).mean(), np.cos(delta).mean())
        return d1 + d2 + 2 * math.sqrt(-2 * math.log(max(length, 1e-9)))


def fit_psi0(target, path, sl, az, el):
    px, py, front = project(path, sl, az, el, target.m)
    tree = cKDTree(np.c_[px[front], py[front]])
    dist, idx = tree.query(target.pts)
    near = dist < 3
    psi = path.psi[sl][front][idx[near]]
    delta = target.hue[near] * TWO_PI - psi
    mean = math.atan2(np.sin(delta).mean(), np.cos(delta).mean())
    spread = math.sqrt(-2 * math.log(max(1e-12, math.hypot(np.sin(delta).mean(), np.cos(delta).mean()))))
    return mean % TWO_PI, spread, int(near.sum())


def search(fig, target, path, chis, oscs, els, azs=np.arange(0, 360, 4.0)):
    best = (99, None)
    for chi_t in chis:
        for n in oscs:
            sl = path.segment(chi_t, n)
            for el in els:
                for az in azs:
                    s = target.score(*project(path, sl, az, el, target.m), path.psi[sl])
                    if s < best[0]:
                        best = (s, (chi_t, n, el, az))
    return best


def refine(target, path, start, free):
    names = ["chi_t", "n", "el", "az"]
    x0 = np.array(start, float)

    def f(v):
        full = x0.copy()
        full[free] = v
        chi_t, n, el, az = full
        if n <= 0.05 or chi_t < 0:
            return 99
        sl = path.segment(chi_t, n)
        return target.score(*project(path, sl, az, el, target.m), path.psi[sl])

    res = minimize(f, x0[free], method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 600})
    full = x0.copy()
    full[free] = res.x
    return res.fun, tuple(full)


def render(fig, target, path, chi_t, n, el, az, psi0, out):
    from matplotlib.colors import hsv_to_rgb
    paper = Image.open(HERE.parent / "figs" / f"{NAMES[fig]}.png").convert("RGB").transpose(Image.FLIP_TOP_BOTTOM)
    paper = paper.resize(tuple(target.m["size"]), Image.LANCZOS)
    w, h = paper.size
    model = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(model)
    m = target.m
    # the spheroid, as a dense cloud of its front points, and the equator's front arc
    class Shape:
        pass
    shape = Shape()
    shape.A, shape.B = path.A, path.B
    u, v = np.meshgrid(np.linspace(0, TWO_PI, 360), np.linspace(0, math.pi, 180))
    shape.xyz = np.c_[(path.A * np.sin(v) * np.cos(u)).ravel(), (path.A * np.sin(v) * np.sin(u)).ravel(),
                      (path.B * np.cos(v)).ravel()]
    sx, sy, sf = project(shape, slice(None), az, el, m)
    for x, y in zip(sx[sf], sy[sf]):
        d.rectangle([x - 1.5, y - 1.5, x + 1.5, y + 1.5], fill=(12, 12, 12))
    t = np.linspace(0, TWO_PI, 720)
    shape.xyz = np.c_[path.A * np.cos(t), path.A * np.sin(t), np.zeros_like(t)]
    rx_, ry_, rf = project(shape, slice(None), az, el, m)
    for k in range(len(t) - 1):
        if rf[k] and rf[k + 1]:
            d.line([(rx_[k], ry_[k]), (rx_[k + 1], ry_[k + 1])], fill=(230, 230, 230), width=2)
    sl = path.segment(chi_t, n)
    px, py, front = project(path, sl, az, el, m)
    hue = ((path.psi[sl] + psi0) % TWO_PI) / TWO_PI
    rgb = (hsv_to_rgb(np.c_[hue, np.ones_like(hue), np.ones_like(hue)]) * 255).astype(int)
    for k in range(len(px) - 1):
        if front[k] and front[k + 1]:
            d.line([(px[k], py[k]), (px[k + 1], py[k + 1])], fill=tuple(rgb[k]), width=3)
    over = Image.blend(paper.convert("L").convert("RGB"), model, 0.5)
    sheet = Image.new("RGB", (3 * w, h), "white")
    sheet.paste(paper, (0, 0))
    sheet.paste(model, (w, 0))
    sheet.paste(over, (2 * w, 0))
    sheet.save(out)


if __name__ == "__main__":
    measured = json.loads((HERE / "measured.json").read_text())
    plan = json.loads(sys.argv[1])  # {fig: {chis: [...], oscs: [...], els: [...], free: [...], max_osc}}
    results = {}
    for fig, spec in plan.items():
        target = Target(measured[fig])
        path = Path3(fig, spec["max_osc"])
        score, best = search(fig, target, path, spec["chis"], spec["oscs"], spec["els"])
        score2, best2 = refine(target, path, best, spec["free"])
        chi_t, n, el, az = best2
        psi0, spread, used = fit_psi0(target, path, path.segment(chi_t, n), az, el)
        sl = path.segment(chi_t, n)
        lam_t, lam_h = path.lam[sl][0], path.lam[sl][-1]
        results[fig] = dict(score=score2, chi_t=chi_t, n=n, el=el, az=az % 360, psi0=psi0, psi_spread=spread,
                            chi_h=chi_t + TWO_PI * n, lam_t=lam_t, lam_h=lam_h, trail=lam_h - lam_t)
        print(f"Fig {fig}: coarse {score:.2f} {tuple(round(v, 3) for v in best)} -> {score2:.2f}: chi_t {chi_t:.3f}, "
              f"osc {n:.3f}, el {el:.1f}, az {az % 360:.1f}, psi0 {psi0:.3f} (spread {spread:.2f}, {used} px), "
              f"lambda {lam_t:.3f}..{lam_h:.3f} (trail {lam_h - lam_t:.3f})", flush=True)
        render(fig, target, path, chi_t, n, el, az, psi0, HERE / f"fit-{fig}.png")
    for fig, result in results.items():
        (HERE / f"fit-{fig}.json").write_text(json.dumps(result, indent=1))
