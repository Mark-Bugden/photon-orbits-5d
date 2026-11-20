// The allowed region 𝒜(a, b): the momenta (Φ, Ψ) that have a spherical orbit.
//
// 𝒜 is a convex "curved diamond" around the origin, so each ray from the origin leaves it once, and its edge
// is found by bisection along rays. A line-by-line port of src/photon_orbits/region.py.

import { orbitExists, sqrt0 } from './orbit.js';

const REACH = 4; // 𝒜 lies inside |Φ|, |Ψ| < 3 for spins up to 1; the search starts from here
const EDGE_STEPS = 48; // bisection steps along a ray: 4/2⁴⁸ ≈ 1.4e-14

/**
 * The corners of 𝒜 on the axes, as [Φ, Ψ, x], for a, b ≥ 0: the equatorial circular orbits (Ψ = 0, eqs
 * 3.9–3.12), counter- then co-rotating, and the polar ones (Φ = 0, eqs 3.13–3.15).
 */
export function vertices(a, b) {
  const sp = sqrt0((1 + a) * (1 + a) - b * b);
  const sm = sqrt0((1 - a) * (1 - a) - b * b);
  const tp = sqrt0((1 + b) * (1 + b) - a * a);
  const tm = sqrt0((1 - b) * (1 - b) - a * a);
  return [
    [1 + sp, 0, 1 + a - b * b + sp],
    [-1 - sm, 0, 1 - a - b * b + sm],
    [0, 1 + tp, 1 + b - a * a + tp],
    [0, -1 - tm, 1 - b - a * a + tm],
  ];
}

/** How far 𝒜 reaches along the unit direction (dPhi, dPsi): the largest s with an orbit at s·(dPhi, dPsi). */
export function edge(a, b, dPhi, dPsi) {
  let lo = 0;
  let hi = REACH;
  while (orbitExists(a, b, hi * dPhi, hi * dPsi)) hi *= 2;
  for (let step = 0; step < EDGE_STEPS; step++) {
    const mid = 0.5 * (lo + hi);
    if (orbitExists(a, b, mid * dPhi, mid * dPsi)) lo = mid;
    else hi = mid;
  }
  return lo;
}

/** rays unit vectors [cos α, sin α], α = 2πk/rays, starting along +Φ and turning towards +Ψ. */
export function directions(rays) {
  return Array.from({ length: rays }, (_, k) => [Math.cos((2 * Math.PI * k) / rays), Math.sin((2 * Math.PI * k) / rays)]);
}

/** The edge of 𝒜 as rays points [Φ, Ψ]. */
export function boundary(a, b, rays = 96) {
  return directions(rays).map(([dPhi, dPsi]) => {
    const s = edge(a, b, dPhi, dPsi);
    return [s * dPhi, s * dPsi];
  });
}
