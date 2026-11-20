// The momentum pad's geometry: how far the allowed region 𝒜(a, b) reaches along a ray or a segment, and where the
// dot goes when the spins change. Pure functions; pad.js draws them and drives them from the pointer and keys.
//
// At an extremal hole (a + b = 1) the co-rotating edge of 𝒜 is the segment Φ + Ψ = −1, where the orbits close on
// the horizon and φ winds without bound (docs/equations.md §8). There the dot stops a hair short of the edge,
// at 1 + Φ + Ψ = GAP, where φ turns about 21 times per oscillation for a = b = ½.

import { edge, orbitExists } from './physics/index.js';

export const GAP = 1e-4;
const STEPS = 48; // bisection steps along a segment, as along a ray in physics/region.js

/** Whether the spins make an extremal hole: a + b = 1, to the sliders' precision. */
export const extremal = (a, b) => a + b > 1 - 1e-9;

/** Whether the pad allows these momenta: they have an orbit, and at an extremal hole they keep GAP from its edge. */
export function allowed(a, b, Phi, Psi) {
  return orbitExists(a, b, Phi, Psi) && !(extremal(a, b) && 1 + Phi + Psi < GAP * (1 - 1e-9));
}

/** How far the pad allows the dot along the unit direction (dPhi, dPsi): the edge of 𝒜, or short of it by GAP. */
export function reach(a, b, dPhi, dPsi) {
  const s = edge(a, b, dPhi, dPsi);
  const toward = dPhi + dPsi;
  return extremal(a, b) && toward < 0 ? Math.min(s, (1 - GAP) / -toward) : s;
}

/** The last allowed point on the segment from `from`, which is allowed, to `to`: `to` itself if it is allowed. */
export function clampAlong(a, b, [Phi0, Psi0], [Phi1, Psi1]) {
  if (allowed(a, b, Phi1, Psi1)) return [Phi1, Psi1];
  let lo = 0;
  let hi = 1;
  for (let step = 0; step < STEPS; step++) {
    const mid = 0.5 * (lo + hi);
    if (allowed(a, b, Phi0 + mid * (Phi1 - Phi0), Psi0 + mid * (Psi1 - Psi0))) lo = mid;
    else hi = mid;
  }
  return [Phi0 + lo * (Phi1 - Phi0), Psi0 + lo * (Psi1 - Psi0)];
}

/** The allowed point nearest (Phi, Psi) along its ray from the origin: itself if it is allowed, else the edge. */
export function toRegion(a, b, Phi, Psi) {
  if (allowed(a, b, Phi, Psi)) return [Phi, Psi];
  const s = Math.hypot(Phi, Psi);
  const dPhi = Phi / s;
  const dPsi = Psi / s;
  const r = reach(a, b, dPhi, dPsi);
  return [r * dPhi, r * dPsi];
}

/**
 * Where (Phi, Psi) sits in the region, as { direction, fraction }: its direction from the origin, a unit vector
 * (null at the origin), and its fraction of the way to the edge along it, from 0 to 1. A point on an axis has a
 * direction exactly along it, so it stays on the axis when the spins change.
 */
export function locate(a, b, Phi, Psi) {
  const s = Math.hypot(Phi, Psi);
  if (s === 0) return { direction: null, fraction: 0 };
  const direction = [Phi / s, Psi / s];
  return { direction, fraction: Math.min(1, s / reach(a, b, ...direction)) };
}

/** The point at this fraction of the way to the edge, along this direction: the inverse of locate. */
export function place(a, b, { direction, fraction }) {
  if (direction === null || fraction === 0) return [0, 0];
  const s = fraction * reach(a, b, ...direction);
  return [s * direction[0], s * direction[1]];
}
