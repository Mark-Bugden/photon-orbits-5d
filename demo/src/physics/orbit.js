// The spherical photon orbit for given spins and momenta: the horizon, the radius and the latitude band.
//
// A line-by-line port of src/photon_orbits/orbit.py, the reference, so it does the same operations in the same
// order; docs/equations.md explains every formula and guard. test/physics.test.js replays the reference's
// results from tests/fixtures/orbits.json. Units μ = 1 and E = 1, so Φ, Ψ and K stand for Φ/E, Ψ/E and K/E².
// x = r², y = sin²θ and w = 1 − y = cos²θ. Polynomials are arrays of coefficients, highest power first.

export const SNAP = 1e-12; // a momentum smaller than this in size counts as exactly 0
export const BAND_TOL = 1e-12; // the band test's tolerance, relative to 1 + |K| + Φ² + Ψ²

/** A momentum, with sizes below SNAP set to exactly 0. */
export const snap = (value) => (Math.abs(value) < SNAP ? 0 : value);

/** The square root, with rounding errors below 0 clamped to 0. */
export const sqrt0 = (value) => (value > 0 ? Math.sqrt(value) : 0);

export function polyval(c, x) {
  let total = 0;
  for (const coefficient of c) total = total * x + coefficient;
  return total;
}

/** The real roots of A t² + B t + C in increasing order; a linear equation if A = 0. */
export function quadraticRoots(A, B, C) {
  if (A === 0) return B === 0 ? [] : [-C / B];
  const disc = B * B - 4 * A * C;
  if (disc < 0) return [];
  const root = Math.sqrt(disc);
  const q = -0.5 * (B < 0 ? B - root : B + root);
  if (q === 0) return [0];
  const t1 = q / A;
  const t2 = C / q;
  return t1 <= t2 ? [t1, t2] : [t2, t1];
}

/**
 * { xPlus, e }: the outer horizon x₊ = r₊² and e = x₊ − x₋, so Δ(x₊ + u) = u(u + e); null if there is no
 * horizon (|a| + |b| > 1, or x₊ = 0). e = 0 for an extremal hole.
 */
export function horizons(a, b) {
  const s = 1 - a * a - b * b;
  const disc = s * s - 4 * a * a * b * b;
  if (s <= 0 || disc < -1e-12) return null; // a disc just below 0 is an extremal hole after rounding
  const e = sqrt0(disc);
  return { xPlus: 0.5 * (s + e), e };
}

/** x₊ = r₊², the outer horizon; null if there is none. */
export function horizon(a, b) {
  const found = horizons(a, b);
  return found === null ? null : found.xPlus;
}

/** P(x), where 𝒳(x) = P(x) − K Δ(x). */
export function pCoefficients(a, b, Phi, Psi) {
  const A2 = a * a;
  const B2 = b * b;
  const p1 = (A2 - B2) * (Phi * Phi - Psi * Psi) + A2 * B2 + A2 + B2 + 2 * a * Phi + 2 * b * Psi;
  const m = a * b + b * Phi + a * Psi;
  const p0 = m * m + (A2 - B2) * (B2 * Phi * Phi - A2 * Psi * Psi);
  return [1, A2 + B2, p1, p0];
}

/**
 * { P, F }: P and F = P′Δ − PΔ′ as polynomials in u = x − x₊, built so that no coefficient loses precision
 * when the orbit is near the horizon (docs/equations.md, "The radius").
 */
export function shiftedCoefficients(a, b, Phi, Psi, xPlus, e) {
  const A2 = a * a;
  const B2 = b * b;
  const ra = xPlus + A2;
  const rb = xPlus + B2;
  const calE = 1 + (a * Phi) / ra + (b * Psi) / rb;
  const drag = (a * Phi * rb) / ra + (b * Psi * ra) / rb; // −(x + a²)(x + b²)𝓔′ at x₊
  const W = xPlus + (A2 - B2) * ((Phi * Phi) / ra - (Psi * Psi) / rb);
  const p0 = ra * rb * calE * calE;
  const p1 = e * W + (2 * xPlus + A2 + B2) * calE * calE - 2 * calE * drag;
  const p2 = 3 * xPlus + A2 + B2;
  return { P: [1, p2, p1, p0], F: [1, 2 * e, e * p2 - p1, -2 * p0, -e * p0] };
}

/** Q(y) = y(1 − y)Θ. With (a, Φ) and (b, Ψ) swapped it is Q(1 − w) as a polynomial in w. */
export function qCoefficients(a, b, Phi, Psi, K) {
  const A2 = a * a;
  const B2 = b * b;
  return [A2 - B2, B2 - 2 * A2 - K, Phi * Phi - Psi * Psi + A2 + K, -Phi * Phi];
}

// How many more roots the quartic F (in u) has at u > 0, besides root: the sign changes of F/(u − root),
// which is positive just above u = 0 and at +∞, between its turning points.
function extraRoots(F, root) {
  const cubic = [1];
  for (const coefficient of F.slice(1, -1)) cubic.push(coefficient + root * cubic[cubic.length - 1]);
  const turning = quadraticRoots(3, 2 * cubic[1], cubic[2]).filter((t) => t > 0);
  const signs = [1, ...turning.map((t) => polyval(cubic, t)), 1].filter((s) => s !== 0);
  let changes = 0;
  for (let i = 1; i < signs.length; i++) if (signs[i - 1] < 0 !== signs[i] < 0) changes++;
  return changes;
}

const noOrbit = (reason, x = null, K = null) => ({ ok: false, reason, x, K });

// { x, K, delta, extra } for the orbit's radius x, the root of F above the horizon, with K = P(x)/Δ(x);
// or a no-orbit result.
function radial(a, b, Phi, Psi) {
  const found = horizons(a, b);
  if (found === null) return noOrbit('no horizon');
  const { xPlus, e } = found;
  const { P, F } = shiftedCoefficients(a, b, Phi, Psi, xPlus, e);
  if (!(P[3] > 0)) return noOrbit('at horizon'); // 𝓔(x₊) = 0: P and Δ both vanish there, and V = P/Δ is 0/0
  // P₀ > 0, so F < 0 just above the horizon (V → +∞ there), and F → +∞: bisect on u in (0, hi).
  let lo = 0;
  let hi = 4;
  while (polyval(F, hi) <= 0) hi *= 2;
  for (let step = 0; step < 200; step++) {
    const mid = 0.5 * (lo + hi);
    if (polyval(F, mid) > 0) hi = mid;
    else lo = mid;
    if (hi - lo <= 1e-15 * hi) break;
  }
  const u = 0.5 * (lo + hi);
  const delta = u * (u + e);
  return { x: xPlus + u, K: polyval(P, u) / delta, delta, extra: extraRoots(F, u) };
}

// Q divided by u if it vanishes at u = 0, and by (1 − u) if it vanishes at u = 1, working up from the
// constant term so that it stays exact.
function reduce(q, atZero, atOne) {
  let r = q.slice();
  if (atZero) r = r.slice(0, -1);
  if (atOne) {
    const s = [r[r.length - 1]];
    for (const coefficient of r.slice(1, -1).reverse()) s.push(coefficient + s[s.length - 1]);
    r = s.reverse();
  }
  return r;
}

// { top: y_m, value: R(y_m), r: R }: where R, which has the sign of Θ, is largest on [0, 1]; null if there is
// nowhere to look. R is Q, divided by y if Φ = 0 and by 1 − y if Ψ = 0.
function bandPeak(a, b, Phi, Psi, K) {
  const r = reduce(qCoefficients(a, b, Phi, Psi, K), Phi === 0, Psi === 0);
  const degree = r.length - 1;
  const slope = [...new Array(3 - degree).fill(0), ...r.slice(0, degree).map((c, i) => c * (degree - i))];
  const candidates = [...(Phi === 0 ? [0] : []), ...(Psi === 0 ? [1] : [])];
  for (const t of quadraticRoots(...slope)) if (t > 0 && t < 1) candidates.push(t);
  if (candidates.length === 0) return null;
  let top = candidates[0];
  let value = polyval(r, top);
  for (const t of candidates.slice(1)) {
    const v = polyval(r, t);
    if (v > value) {
      top = t;
      value = v;
    }
  }
  return { top, value, r };
}

// The root of c in (0, hi), where c(0) < 0 < c(hi), by bisection to full relative precision.
function lowerRoot(c, hi) {
  let lo = 0;
  for (let step = 0; step < 1100; step++) {
    const mid = 0.5 * (lo + hi);
    if (polyval(c, mid) < 0) lo = mid;
    else hi = mid;
    if (hi - lo <= 1e-15 * hi) break;
  }
  return 0.5 * (lo + hi);
}

/**
 * { y1, w2 } = (sin²θ_min, cos²θ_max), the edges of the band where Θ ≥ 0; null if there is no band. The upper
 * edge comes from the same cubic with (a, Φ) and (b, Ψ) swapped, so both edges keep full relative precision.
 */
export function band(a, b, Phi, Psi, K) {
  const peak = bandPeak(a, b, Phi, Psi, K);
  if (peak === null) return null;
  const { top, value, r } = peak;
  const tol = BAND_TOL * (1 + Math.abs(K) + Phi * Phi + Psi * Psi);
  if (value < -tol) return null;
  if (value <= tol) return { y1: top, w2: 1 - top }; // a circular orbit at constant θ
  const flipped = reduce(qCoefficients(b, a, Psi, Phi, K), Psi === 0, Phi === 0);
  const y1 = polyval(r, 0) >= 0 ? 0 : lowerRoot(r, top);
  const w2 = polyval(flipped, 0) >= 0 ? 0 : lowerRoot(flipped, 1 - top);
  return { y1, w2 };
}

/** Whether these momenta have a spherical orbit: the same decision as solveOrbit, without the band's edges. */
export function orbitExists(a, b, Phi, Psi) {
  Phi = snap(Phi);
  Psi = snap(Psi);
  const found = radial(a, b, Phi, Psi);
  if (found.ok === false) return false;
  const peak = bandPeak(a, b, Phi, Psi, found.K);
  return peak !== null && peak.value >= -BAND_TOL * (1 + Math.abs(found.K) + Phi * Phi + Psi * Psi);
}

/**
 * The spherical photon orbit with momenta (Φ, Ψ) around the hole with spins (a, b): { ok: true, … } with the
 * fields of the reference's Orbit, or { ok: false, reason, x, K }, where reason is 'no horizon', 'at horizon'
 * (the orbit would sit on the horizon) or 'no band' (θ-motion is forbidden: outside the allowed region).
 */
export function solveOrbit(a, b, Phi, Psi) {
  Phi = snap(Phi);
  Psi = snap(Psi);
  const found = radial(a, b, Phi, Psi);
  if (found.ok === false) return found;
  const { x, K, delta, extra } = found;
  const edges = band(a, b, Phi, Psi, K);
  if (edges === null) return noOrbit('no band', x, K);
  if (extra) console.warn(`F has ${extra} more roots above the horizon for a=${a}, b=${b}, Phi=${Phi}, Psi=${Psi}`);
  const { y1, w2 } = edges;
  const y2 = 1 - w2;
  const A2 = a * a;
  const B2 = b * b;
  const calE = 1 + (a * Phi) / (x + A2) + (b * Psi) / (x + B2);
  const Aphi = (a * (x + B2) * calE) / delta + ((A2 - B2) * Phi) / (x + A2);
  const Apsi = (b * (x + A2) * calE) / delta - ((A2 - B2) * Psi) / (x + B2);
  const G1 = (B2 - A2) * (y1 + y1 + y2) - (B2 - 2 * A2 - K);
  const G2 = (B2 - A2) * (y2 + y1 + y2) - (B2 - 2 * A2 - K);
  return { ok: true, a, b, Phi, Psi, xPlus: horizon(a, b), x, K, y1, y2, w2, delta, calE, Aphi, Apsi, G1, G2, extraRoots: extra };
}
