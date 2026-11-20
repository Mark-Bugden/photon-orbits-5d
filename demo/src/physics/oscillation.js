// One latitude oscillation of an orbit, tabulated in the angle variable χ, and the path drawn from it.
//
// With Mino time σ (dσ = dλ/ρ²) and y = sin²θ = ȳ − Δy cos χ, every rate depends on χ alone, so one oscillation
// (χ from 0 to 2π) is a quadrature and longer paths repeat it. χ = 0 is θ_min and χ = π is θ_max. A line-by-line
// port of src/photon_orbits/oscillation.py; docs/equations.md gives the derivations.

const TWO_PI = 2 * Math.PI;
export const AXIS_NEAR = 1e-4; // a band within this of an axis, in sin²θ or cos²θ, gets 512 steps instead of 256

export const stepsFor = (orbit) => (Math.min(orbit.y1, orbit.w2) < AXIS_NEAR ? 512 : 256);

/**
 * ∫₀^χ dχ′ / (at0 + (atpi − at0) sin²(χ′/2)) at χ = 2πi/n (n even, 0 ≤ i ≤ n), for at0, atpi > 0: the closed form
 * (2/√(at0·atpi)) arctan(k tan(χ/2)), k = √(atpi/at0), evaluated from the index i so that the values at
 * 0, π and 2π are exact even when the integrand is a spike of width ~1e-12 there. i need not be an integer:
 * between uses it between the nodes.
 */
export function axisIntegral(at0, atpi, i, n) {
  const scale = 2 / Math.sqrt(at0 * atpi);
  const k = Math.sqrt(atpi / at0);
  const j = 2 * i <= n ? i : n - i; // fold (π, 2π] onto [0, π)
  const angle =
    4 * j <= n
      ? Math.atan(k * Math.tan((Math.PI * j) / n))
      : 0.5 * Math.PI - Math.atan(Math.tan((Math.PI * (n - 2 * j)) / (2 * n)) / k);
  return j === i ? scale * angle : scale * (Math.PI - angle);
}

/** d(σ, λ, φ, ψ)/dχ at χ, less the axis terms Φ/(2y√G₁) and Ψ/(2(1 − y)√G₂), which tabulate integrates exactly. */
export function rates(orbit, chi) {
  const A2 = orbit.a * orbit.a;
  const B2 = orbit.b * orbit.b;
  const s = Math.sin(0.5 * chi);
  const c = Math.cos(0.5 * chi);
  const width = orbit.y2 - orbit.y1;
  const rise = width * s * s; // y − y1
  const fall = width * c * c; // y2 − y
  const y = orbit.y1 + rise;
  const w = orbit.w2 + fall; // 1 − y, without the rounding of 1 − y near the equator
  const root = Math.sqrt((B2 - A2) * (y + orbit.y1 + orbit.y2) - (B2 - 2 * A2 - orbit.K));
  const dsigma = 0.5 / root;
  const dlam = (orbit.x + A2 + (B2 - A2) * y) * dsigma;
  let dphi = -orbit.Aphi * dsigma;
  if (orbit.Phi !== 0) {
    const root1 = Math.sqrt(orbit.G1);
    dphi -= (orbit.Phi * (B2 - A2) * rise) / (2 * y * root * root1 * (root + root1));
  }
  let dpsi = -orbit.Apsi * dsigma;
  if (orbit.Psi !== 0) {
    const root2 = Math.sqrt(orbit.G2);
    dpsi += (orbit.Psi * (B2 - A2) * fall) / (2 * w * root * root2 * (root + root2));
  }
  return [dsigma, dlam, dphi, dpsi];
}

/**
 * One oscillation by Simpson's rule on n steps (256, or 512 near an axis), with the axis terms in closed form:
 * { n, sigma, lam, phi, psi }, each with the values at the n + 1 nodes χ_i = 2πi/n, from 0 at θ_min.
 */
export function tabulate(orbit, n = stepsFor(orbit)) {
  if (n % 2) throw new RangeError('the number of steps must be even');
  const h = TWO_PI / n;
  const f = [];
  for (let j = 0; j <= 2 * n; j++) f.push(rates(orbit, 0.5 * h * j));
  const columns = [0, 1, 2, 3].map(() => new Float64Array(n + 1));
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < 4; k++) {
      columns[k][i + 1] = columns[k][i] + (h / 6) * (f[2 * i][k] + 4 * f[2 * i + 1][k] + f[2 * i + 2][k]);
    }
  }
  const [sigma, lam, phi, psi] = columns;
  if (orbit.Phi !== 0) {
    const scale = orbit.Phi / (2 * Math.sqrt(orbit.G1));
    for (let i = 0; i <= n; i++) phi[i] += scale * axisIntegral(orbit.y1, orbit.y2, i, n);
  }
  if (orbit.Psi !== 0) {
    const scale = orbit.Psi / (2 * Math.sqrt(orbit.G2));
    for (let i = 0; i <= n; i++) psi[i] += scale * axisIntegral(1 - orbit.y1, orbit.w2, i, n);
  }
  return { n, sigma, lam, phi, psi };
}

/** { chi, sigma, lam, phi, psi } at node i, for any integer i: oscillation m = ⌊i/n⌋ is the table shifted by m increments. */
export function node(table, i) {
  const { n } = table;
  const m = Math.floor(i / n);
  const j = i - m * n;
  return {
    chi: j * (TWO_PI / n) + m * TWO_PI,
    sigma: table.sigma[j] + m * table.sigma[n],
    lam: table.lam[j] + m * table.lam[n],
    phi: table.phi[j] + m * table.phi[n],
    psi: table.psi[j] + m * table.psi[n],
  };
}

/**
 * { chi, sigma, lam, phi, psi } at χ = (i + t)·2π/n, a fraction 0 ≤ t ≤ 1 of the way from node i to node i + 1:
 * one Simpson step of width t·h from node i, with the axis terms in closed form. It is node i at t = 0 and meets
 * node i + 1 at t = 1. Near an axis φ or ψ swings by about π between two nodes, and this follows the swing, so
 * the drawing can subdivide it.
 */
export function between(orbit, table, i, t) {
  const { n } = table;
  const m = Math.floor(i / n);
  const j = i - m * n;
  const h = TWO_PI / n;
  const chi = (j + t) * h;
  const f0 = rates(orbit, j * h);
  const f1 = rates(orbit, (j + 0.5 * t) * h);
  const f2 = rates(orbit, chi);
  const [sigma, lam, phi0, psi0] = [table.sigma, table.lam, table.phi, table.psi].map(
    (column, k) => column[j] + ((t * h) / 6) * (f0[k] + 4 * f1[k] + f2[k]),
  );
  let phi = phi0;
  let psi = psi0;
  if (orbit.Phi !== 0) {
    const scale = orbit.Phi / (2 * Math.sqrt(orbit.G1));
    phi += scale * (axisIntegral(orbit.y1, orbit.y2, j + t, n) - axisIntegral(orbit.y1, orbit.y2, j, n));
  }
  if (orbit.Psi !== 0) {
    const scale = orbit.Psi / (2 * Math.sqrt(orbit.G2));
    psi += scale * (axisIntegral(1 - orbit.y1, orbit.w2, j + t, n) - axisIntegral(1 - orbit.y1, orbit.w2, j, n));
  }
  return {
    chi: chi + m * TWO_PI,
    sigma: sigma + m * table.sigma[n],
    lam: lam + m * table.lam[n],
    phi: phi + m * table.phi[n],
    psi: psi + m * table.psi[n],
  };
}

/** [sin θ, cos θ] at χ, signed where the orbit passes through an axis, so that it carries straight on. */
export function latitude(orbit, chi) {
  const s = Math.sin(0.5 * chi);
  const c = Math.cos(0.5 * chi);
  const width = orbit.y2 - orbit.y1;
  const throughPole = orbit.Phi === 0 && orbit.y1 === 0;
  const throughEquator = orbit.Psi === 0 && orbit.w2 === 0;
  const sinTheta = throughPole ? Math.sqrt(orbit.y2) * s : Math.sqrt(orbit.y1 + width * s * s);
  const cosTheta = throughEquator ? Math.sqrt(1 - orbit.y1) * c : Math.sqrt(orbit.w2 + width * c * c);
  return [sinTheta, cosTheta];
}

/** [X, Y, Z] in the paper's picture: the a-plane flat, and the b-plane's radius as height. Colour it by ψ. */
export function embed(orbit, chi, phi) {
  const [sinTheta, cosTheta] = latitude(orbit, chi);
  const across = Math.sqrt(orbit.x + orbit.a * orbit.a) * sinTheta;
  return [across * Math.cos(phi), across * Math.sin(phi), Math.sqrt(orbit.x + orbit.b * orbit.b) * cosTheta];
}

/** [ν_φ, ν_ψ]: the turns of φ and ψ per oscillation, counting a pass through an axis as half a turn. */
export function frequencies(orbit, table) {
  const dphi = table.phi[table.n] + (orbit.Phi === 0 && orbit.y1 === 0 ? Math.PI : 0);
  const dpsi = table.psi[table.n] + (orbit.Psi === 0 && orbit.w2 === 0 ? Math.PI : 0);
  return [dphi / TWO_PI, dpsi / TWO_PI];
}
