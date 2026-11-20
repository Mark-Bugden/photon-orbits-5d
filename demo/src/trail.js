// The photon's path, ready to draw: one oscillation refined so that no drawn step turns φ or ψ by more than 10°,
// and the stretch of it the trail shows, repeated oscillation by oscillation. Pure functions with no three.js;
// scene/path.js puts what fill() writes on screen, and panel/ what trace() visits.

import { between, embed, frequencies, latitude } from './physics/index.js';

const PI = Math.PI;
const TWO_PI = 2 * PI;
export const BASE_STEPS = 128; // drawn steps per oscillation before refining: plenty for θ, which is smooth in χ
export const MAX_TURN = Math.PI / 18; // 10°: the most φ or ψ may turn in one drawn step
const MAX_DEPTH = 48; // halvings of a step: enough to follow a swing near an axis 1e-12 wide

// The state at u = i + t along the table, for 0 ≤ u ≤ n.
function at(orbit, table, u) {
  const i = Math.floor(u);
  return between(orbit, table, i, u - i);
}

// Pushes the states strictly between u0 and u1 that keep each step's turn of φ and ψ within MAX_TURN, in order.
function split(orbit, table, u0, s0, u1, s1, depth, states) {
  if (depth === MAX_DEPTH) return;
  if (Math.abs(s1.phi - s0.phi) <= MAX_TURN && Math.abs(s1.psi - s0.psi) <= MAX_TURN) return;
  const u = 0.5 * (u0 + u1);
  const middle = at(orbit, table, u);
  split(orbit, table, u0, s0, u, middle, depth + 1, states);
  states.push(middle);
  split(orbit, table, u, middle, u1, s1, depth + 1, states);
}

/**
 * One oscillation, from θ_min (χ = 0) to the next θ_min (χ = 2π), refined: { count, chi, lam, theta, phi, psi, x,
 * y, z, period, phiAdvance, psiAdvance, turn, pole, equator, band }. Point k, for 0 ≤ k ≤ count, is at χ = chi[k]
 * and λ = lam[k], with angles theta[k], phi[k] and psi[k] and position (x, y, z)[k]. Oscillation m is this one moved
 * on by m oscillations: λ by m·period, φ and ψ by m·phiAdvance and m·psiAdvance. In the picture it is turned about
 * the axis by m·turn and, if the orbit passes through the equator, has z changing sign when m is odd.
 *
 * θ is signed where the orbit passes through an axis, as the picture draws it (physics/ latitude): through the
 * pole (pole) it runs from 0 at χ = 0 to θ_max and back, and through the equator (equator) from θ_min past π/2 to
 * π − θ_min. band is [θ_min, θ_max], the edges of the band of latitudes.
 */
export function refine(orbit, table) {
  const stride = table.n / BASE_STEPS;
  let previous = at(orbit, table, 0);
  const states = [previous];
  for (let s = 1; s <= BASE_STEPS; s++) {
    const next = at(orbit, table, s * stride);
    split(orbit, table, (s - 1) * stride, previous, s * stride, next, 0, states);
    states.push(next);
    previous = next;
  }
  const count = states.length - 1;
  const columns = ['chi', 'lam', 'theta', 'phi', 'psi', 'x', 'y', 'z'];
  const path = Object.fromEntries(columns.map((name) => [name, new Float64Array(count + 1)]));
  states.forEach((state, k) => {
    const [X, Y, Z] = embed(orbit, state.chi, state.phi);
    path.chi[k] = state.chi;
    path.lam[k] = state.lam;
    path.theta[k] = Math.atan2(...latitude(orbit, state.chi));
    path.phi[k] = state.phi;
    path.psi[k] = state.psi;
    path.x[k] = X;
    path.y[k] = Y;
    path.z[k] = Z;
  });
  return Object.assign(path, {
    count,
    period: table.lam[table.n],
    phiAdvance: table.phi[table.n],
    psiAdvance: table.psi[table.n],
    // Δφ, plus half a turn if the orbit passes over the pole, where sin θ changes sign (latitude in physics/).
    turn: TWO_PI * frequencies(orbit, table)[0],
    // An orbit through the pole passes over it each oscillation, at χ = 0: sin θ changes sign.
    pole: orbit.Phi === 0 && orbit.y1 === 0,
    // An orbit through the equator crosses it each oscillation, at χ = π: cos θ changes sign.
    equator: orbit.Psi === 0 && orbit.w2 === 0,
    band: [Math.asin(Math.sqrt(orbit.y1)), Math.acos(Math.sqrt(orbit.w2))],
  });
}

// The index k of the last entry of values (increasing, from values[0] = 0) at or below v, at most count − 1.
function search(values, count, v) {
  let lo = 0;
  let hi = count;
  while (hi - lo > 1) {
    const middle = (lo + hi) >> 1;
    if (values[middle] <= v) lo = middle;
    else hi = middle;
  }
  return lo;
}

/** Where λ falls on the path, written into place: { m, k, t }, a fraction t of the way from point k to k + 1 of oscillation m. */
export function locate(path, lam, place = {}) {
  const m = Math.floor(lam / path.period);
  const within = lam - m * path.period;
  const k = search(path.lam, path.count, within);
  const span = path.lam[k + 1] - path.lam[k];
  place.m = m;
  place.k = k;
  place.t = span > 0 ? Math.min(1, Math.max(0, (within - path.lam[k]) / span)) : 0;
  return place;
}

/** χ at λ: where the photon is in its oscillation, which stays put when the orbit changes. */
export function chiAt(path, lam) {
  const { m, k, t } = locate(path, lam);
  return m * TWO_PI + path.chi[k] + t * (path.chi[k + 1] - path.chi[k]);
}

/** λ at χ, the inverse of chiAt. */
export function lamAt(path, chi) {
  const m = Math.floor(chi / TWO_PI);
  const within = chi - m * TWO_PI;
  const k = search(path.chi, path.count, within);
  const span = path.chi[k + 1] - path.chi[k];
  const t = span > 0 ? Math.min(1, Math.max(0, (within - path.chi[k]) / span)) : 0;
  return m * path.period + path.lam[k] + t * (path.lam[k + 1] - path.lam[k]);
}

// Scratch space, so that fill and trace allocate nothing per frame.
const tailPlace = {};
const headPlace = {};
const reach = { first: 0, firstT: 0, last: 0, lastT: 0 };
const turned = { m: NaN, cos: 1, sin: 0, z: 1 };
const point = new Float64Array(5); // x, y, z, λ, ψ
const moved = { m: NaN, lam: 0, phi: 0, psi: 0, sign: 1, offset: 0 };

// The stretch of the endless path from λ = head − length to head, written into reach: from a fraction firstT of
// the way along step first to a fraction lastT along step last, where step g runs from point g to point g + 1 and
// g = m·count + k. At most capacity steps: if it needs more, its oldest part is left out.
function span(path, head, length, capacity) {
  locate(path, head - length, tailPlace);
  locate(path, head, headPlace);
  reach.first = tailPlace.m * path.count + tailPlace.k;
  reach.firstT = tailPlace.t;
  reach.last = headPlace.m * path.count + headPlace.k;
  reach.lastT = headPlace.t;
  if (reach.last - reach.first + 1 > capacity) {
    reach.first = reach.last - capacity + 1;
    reach.firstT = 0;
  }
  return reach;
}

// Writes into `point` the position, λ and ψ a fraction t of the way from point g to point g + 1 of the endless path,
// where g = m·count + k.
function pointAt(path, g, t) {
  const m = Math.floor(g / path.count);
  const k = g - m * path.count;
  if (turned.m !== m) {
    turned.m = m;
    turned.cos = Math.cos(m * path.turn);
    turned.sin = Math.sin(m * path.turn);
    turned.z = path.equator && m % 2 !== 0 ? -1 : 1;
  }
  const x = path.x[k] + t * (path.x[k + 1] - path.x[k]);
  const y = path.y[k] + t * (path.y[k + 1] - path.y[k]);
  point[0] = turned.cos * x - turned.sin * y;
  point[1] = turned.sin * x + turned.cos * y;
  point[2] = turned.z * (path.z[k] + t * (path.z[k + 1] - path.z[k]));
  point[3] = m * path.period + path.lam[k] + t * (path.lam[k + 1] - path.lam[k]);
  point[4] = m * path.psiAdvance + path.psi[k] + t * (path.psi[k + 1] - path.psi[k]);
}

/**
 * Writes the trail, the path from λ = head − length up to the photon at λ = head, as line segments: positions and
 * colours each get 6 floats a segment (start, then end). paint(ψ, fraction, colours, offset) writes a point's
 * colour, where fraction runs from 0 at the tail to 1 at the photon. If the trail needs more segments than the
 * arrays hold, its oldest part is left out. Returns the number of segments, and leaves the photon's position and
 * ψ in photon = [x, y, z, ψ].
 */
export function fill(path, head, length, positions, colours, paint, photon) {
  const capacity = Math.floor(positions.length / 6);
  const tail = head - length;
  turned.m = NaN; // the path may have changed since the last call
  const { first, firstT, last, lastT } = span(path, head, length, capacity);
  const fraction = () => (length > 0 ? Math.min(1, Math.max(0, (point[3] - tail) / length)) : 1);

  let segments = 0;
  pointAt(path, first, firstT);
  positions.set(point.subarray(0, 3), 0);
  paint(point[4], fraction(), colours, 0);
  const finish = () => {
    const end = 6 * segments + 3;
    positions[end] = point[0];
    positions[end + 1] = point[1];
    positions[end + 2] = point[2];
    paint(point[4], fraction(), colours, end);
    segments++;
    if (segments < capacity) {
      positions.copyWithin(6 * segments, end, end + 3);
      colours.copyWithin(6 * segments, end, end + 3);
    }
  };
  for (let g = first + 1; g <= last; g++) {
    pointAt(path, g, 0);
    finish();
  }
  pointAt(path, last, lastT);
  finish();
  photon[0] = point[0];
  photon[1] = point[1];
  photon[2] = point[2];
  photon[3] = point[4];
  return segments;
}

// Calls visit(λ, θ, φ, ψ) a fraction t of the way from point g to point g + 1 of the endless path. Oscillation m
// is oscillation 0 moved on: λ, φ and ψ by m times their advances, and θ carried on through the axes. Through the
// pole, θ changes sign each oscillation; through the equator, θ ↦ π − θ; through both, θ = χ/2 and grows by π.
function visitAt(path, g, t, visit) {
  const m = Math.floor(g / path.count);
  const k = g - m * path.count;
  if (moved.m !== m) {
    const odd = m % 2 !== 0;
    moved.m = m;
    moved.lam = m * path.period;
    moved.phi = m * path.phiAdvance;
    moved.psi = m * path.psiAdvance;
    moved.sign = path.pole !== path.equator && odd ? -1 : 1;
    moved.offset = path.pole && path.equator ? m * PI : path.equator && odd ? PI : 0;
  }
  visit(
    moved.lam + path.lam[k] + t * (path.lam[k + 1] - path.lam[k]),
    moved.sign * (path.theta[k] + t * (path.theta[k + 1] - path.theta[k])) + moved.offset,
    moved.phi + path.phi[k] + t * (path.phi[k + 1] - path.phi[k]),
    moved.psi + path.psi[k] + t * (path.psi[k + 1] - path.psi[k]),
  );
}

/**
 * Visits the trail from λ = head − length to the photon at λ = head, oldest first: visit(λ, θ, φ, ψ) at its tail,
 * at the points of the refined path between, and at the photon. φ and ψ are unwrapped, and θ is signed through
 * the axes as the picture draws it, so all three are continuous. A trail of more than `most` points is visited at
 * every second, third, … point, so as to visit about `most` of them: near the extremal edge, a long trail has
 * hundreds of thousands.
 */
export function trace(path, head, length, visit, most = Infinity) {
  moved.m = NaN; // the path may have changed since the last call
  const { first, firstT, last, lastT } = span(path, head, length, Infinity);
  const stride = Math.max(1, Math.ceil((last - first) / most));
  visitAt(path, first, firstT, visit);
  for (let g = first + stride; g <= last; g += stride) visitAt(path, g, 0, visit);
  visitAt(path, last, lastT, visit);
}
