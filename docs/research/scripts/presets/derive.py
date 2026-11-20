"""The presets in the demo's terms, from the fits, with the λ-plot fits of Figs 4.9 and 4.11 as a cross-check."""

import json
import math

from fit import HERE, PARAMS, TWO_PI, Path3
from photon_orbits import frequencies, solve_orbit, tabulate


def chi_at_theta(orbit, theta, rising):
    """χ in [0, 2π) where θ has this value, rising (towards θ_max, χ < π) or falling."""
    y = math.sin(theta) ** 2
    s = math.sqrt(max(0.0, min(1.0, (y - orbit.y1) / (orbit.y2 - orbit.y1))))
    half = math.asin(s)
    return 2 * half if rising else TWO_PI - 2 * half


def state(path, chi):
    """λ, φ, ψ at χ, by interpolating the dense samples."""
    k = chi / TWO_PI * 720
    i = int(k)
    t = k - i
    pick = lambda column: column[i] + t * (column[i + 1] - column[i])
    return pick(path.lam), pick(path.phi), pick(path.psi)


fits = {fig: json.loads((HERE / f"fit-{fig}.json").read_text()) for fig in PARAMS}
for fig, f in fits.items():
    path = Path3(fig, f["chi_h"] / TWO_PI + 1)
    o = path.orbit
    nu = frequencies(o, path.table)
    print(f"Fig {fig}: osc {f['n']:.3f}, chi {f['chi_t']:.3f} -> {f['chi_h']:.3f}, trail {f['trail']:.3f}, "
          f"view az {f['az']:.1f} el {f['el']:.1f}, psi0 {f['psi0']:.3f}; nu_phi {nu[0]:.4f}, "
          f"turns of phi shown {f['n'] * nu[0]:.3f}")

print()
for fig, theta, rising, lam_paper, phi_paper, psi_paper, lam_end in [
        ("4.8", 0.3924, True, 0.0, 1.7476, 0.4490, None), ("4.10", 1.0214, False, 1.0, -0.9799, 2.9571, 5.0)]:
    f = fits[fig]
    path = Path3(fig, f["chi_h"] / TWO_PI + 2)
    o = path.orbit
    chi_s = chi_at_theta(o, theta, rising)
    lam_s, phi_s, psi_s = state(path, chi_s)
    twist = math.atan(o.a / math.sqrt(o.x))
    print(f"Fig {fig}: the λ-plot fit's start (paper λ = {lam_paper}) is our χ = {chi_s:.3f} (3D fit's tail: "
          f"{f['chi_t']:.3f}), λ = {lam_s:.3f}")
    print(f"  colour: the λ-plot's ψ there is {psi_paper:.3f}; the 3D fit gives {(psi_s + f['psi0']) % TWO_PI:.3f}")
    rotate = (phi_paper - twist) - phi_s
    rotate_plain = phi_paper - phi_s
    print(f"  the paper's picture is ours turned by {math.degrees(rotate) % 360:.1f}° with the twist "
          f"({math.degrees(twist):.1f}°), {math.degrees(rotate_plain) % 360:.1f}° without; so the paper's camera "
          f"azimuth is {(f['az'] + math.degrees(rotate)) % 360:.1f}° with it, {(f['az'] + math.degrees(rotate_plain)) % 360:.1f}° without")
    if lam_end is not None:
        print(f"  the λ-plot's range ends {lam_end - lam_paper:.3f} later: head χ for that: ", end="")
        target = lam_s + (lam_end - lam_paper)
        k = next(i for i in range(len(path.lam)) if path.lam[i] >= target)
        print(f"{path.chi[k]:.3f} (3D fit's head {f['chi_h']:.3f})")
