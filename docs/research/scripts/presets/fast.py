"""A fast first search for figures that show many oscillations: oscillation k is oscillation 0 turned by k·turn
(and flipped in z, for orbits through the equator), so its projection at each azimuth is precomputed once and the
score of any segment is a sum of cumulative sums. Then the best candidates are rescored exactly and refined."""

import json
import math
import sys

import numpy as np
from scipy import ndimage

from fit import HERE, TWO_PI, Path3, Target, fit_psi0, project, refine, render
from photon_orbits import frequencies

J = 720  # samples per oscillation: Path3.PER


def main(fig, spec):
    measured = json.loads((HERE / "measured.json").read_text())
    measured[fig]["rx"] *= spec.get("rx_scale", 1)  # the drawn sphere's outline, as a fraction of the orbit's
    target = Target(measured[fig])
    m = target.m
    max_osc = spec["max_osc"]
    path = Path3(fig, max_osc + 1)
    orbit, table = path.orbit, path.table
    turn = math.degrees(TWO_PI * frequencies(orbit, table)[0])
    dpsi = table.psi[-1]
    flip = orbit.Psi == 0 and orbit.w2 == 0
    # the hue of the nearest paper path pixel, everywhere
    mask = np.zeros((target.h, target.w), bool)
    mask[m["path"]["y"], m["path"]["x"]] = True
    hue_img = np.zeros_like(target.dt)
    hue_img[m["path"]["y"], m["path"]["x"]] = target.hue * TWO_PI
    _, (iy, ix) = ndimage.distance_transform_edt(~mask, return_indices=True)
    nearest_hue = hue_img[iy, ix]

    base = slice(0, J)  # oscillation 0, samples j = 0 … J − 1
    psi0 = path.psi[base]
    best_all = []
    for el in spec["els"]:
        tables = {}
        for parity in ([0, 1] if flip else [0]):
            xyz = path.xyz[base].copy()
            if parity:
                xyz[:, 2] *= -1

            class One:
                pass
            one = One()
            one.xyz, one.A, one.B = xyz, path.A, path.B
            d = np.zeros((360, J))
            w = np.zeros((360, J))
            c = np.zeros((360, J), complex)
            for az in range(360):
                px, py, front = project(one, slice(None), az, el, m)
                xi = np.clip(np.round(px).astype(int), 0, target.w - 1)
                yi = np.clip(np.round(py).astype(int), 0, target.h - 1)
                f = front.astype(float)
                w[az] = f
                d[az] = f * np.minimum(target.dt[yi, xi], 12)
                c[az] = f * np.exp(1j * (nearest_hue[yi, xi] - psi0))
            cum = lambda a: np.concatenate([np.zeros((360, 1), a.dtype), np.cumsum(a, axis=1)], axis=1)
            tables[parity] = (cum(d), cum(w), cum(c))
        az0 = np.arange(360)
        step_t = J // spec.get("chi_steps", 48)
        step_n = J // spec.get("n_steps", 16)
        for s0 in range(0, J, step_t):
            D = np.zeros(360)
            W = np.zeros(360)
            C = np.zeros(360, complex)
            g = s0
            end = s0 + int(spec["min_osc"] * J)
            last = s0 + int(max_osc * J)
            while g < last:
                k, j = divmod(g, J)
                cd, cw, cc = tables[k % 2 if flip else 0]
                rows = np.round(az0 - k * turn).astype(int) % 360
                stop = min(J, j + (last - g))
                # evaluate every segment end that falls in this oscillation
                for e in range(((g // step_n) + 1) * step_n, k * J + stop + 1, step_n):
                    if e < end:
                        continue
                    je = e - k * J
                    Dt = D + cd[rows, je] - cd[rows, j]
                    Wt = W + cw[rows, je] - cw[rows, j]
                    Ct = C + (cc[rows, je] - cc[rows, j]) * np.exp(-1j * k * dpsi)
                    Wt = np.maximum(Wt, 1)
                    score = Dt / Wt + spec.get("colour", 4) * (1 - np.abs(Ct) / Wt)
                    i = int(np.argmin(score))
                    best_all.append((float(score[i]), s0 * TWO_PI / J, (e - s0) / J, el, float(i)))
                D += cd[rows, stop] - cd[rows, j]
                W += cw[rows, stop] - cw[rows, j]
                C += (cc[rows, stop] - cc[rows, j]) * np.exp(-1j * k * dpsi)
                g = k * J + stop
    best_all.sort()
    # the quick score has no term for covering the paper's path, so it favours short segments: keep the best few
    # for each length, rescore them exactly, with coverage both ways, then refine the best of those
    per_n = {}
    for entry in best_all:
        per_n.setdefault(round(entry[2], 4), []).append(entry)
    shortlist = [e for entries in per_n.values() for e in entries[:spec.get("keep", 3)]]
    rescored = []
    for score, chi_t, n, el, az in shortlist:
        sl = path.segment(chi_t, n)
        rescored.append((target.score(*project(path, sl, az, el, m), path.psi[sl]), (chi_t, n, el, az), score))
    rescored.sort(key=lambda r: r[0])
    for r in rescored[:8]:
        print(f"  {r[0]:.2f} (fast {r[2]:.2f}): chi_t {r[1][0]:.3f}, osc {r[1][1]:.3f}, el {r[1][2]}, az {r[1][3]:.0f}")
    score2, (chi_t, n, el, az) = refine(target, path, rescored[0][1], [0, 1, 2, 3])
    psi_0, spread, used = fit_psi0(target, path, path.segment(chi_t, n), az, el)
    sl = path.segment(chi_t, n)
    lam_t, lam_h = path.lam[sl][0], path.lam[sl][-1]
    result = dict(rx_scale=spec.get("rx_scale", 1), score=score2, chi_t=chi_t, n=n, el=el, az=az % 360, psi0=psi_0, psi_spread=spread,
                  chi_h=chi_t + TWO_PI * n, lam_t=lam_t, lam_h=lam_h, trail=lam_h - lam_t)
    print(f"Fig {fig}: {score2:.2f}: chi_t {chi_t:.3f}, osc {n:.3f}, el {el:.1f}, az {az % 360:.1f}, psi0 {psi_0:.3f} "
          f"(spread {spread:.2f}, {used} px), lambda {lam_t:.3f}..{lam_h:.3f} (trail {lam_h - lam_t:.3f})", flush=True)
    render(fig, target, path, chi_t, n, el, az, psi_0, HERE / f"fit-{fig}.png")
    (HERE / f"fit-{fig}.json").write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    for fig, spec in json.loads(sys.argv[1]).items():
        main(fig, spec)
