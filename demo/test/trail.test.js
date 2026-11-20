import assert from 'node:assert/strict';
import { test } from 'node:test';

import { WHEEL, WHEEL_SIZE, hue, wheelIndex } from '../src/colour.js';
import { embed, solveOrbit, tabulate } from '../src/physics/index.js';
import { BASE_STEPS, MAX_TURN, chiAt, fill, lamAt, locate, refine, trace } from '../src/trail.js';

const TWO_PI = 2 * Math.PI;

const CASES = {
  'Fig 4.7': [0.6, 0.3, 1, -1],
  'Fig 4.10, near both axes': [0.9, 0.09, 0.1, 0.01],
  'through the pole (Fig 4.4)': [0.8, 0, 0, 0.5],
  'through the equator': [0.6, 0.3, 1, 0],
  'through both (Fig 4.2)': [0.5, 0, 0, 0],
  'a swing 1e-8 from the pole': [0.6, 0.3, 1e-8, 0.5],
  'a swing 1e-9 from the equator': [0.9, 0.09, 0.5, 1e-9],
  'near the extremal edge, φ winding fast': [0.5, 0.5, -0.3, -0.6999],
};

const prepare = (args) => {
  const orbit = solveOrbit(...args);
  const table = tabulate(orbit);
  return { orbit, table, path: refine(orbit, table) };
};

for (const [name, args] of Object.entries(CASES)) {
  test(`refines one oscillation so that no step turns φ or ψ by more than 10°: ${name}`, () => {
    const { orbit, table, path } = prepare(args);
    assert.ok(path.count >= BASE_STEPS);
    assert.equal(path.chi[0], 0);
    assert.equal(path.chi[path.count], TWO_PI);
    assert.equal(path.lam[path.count], table.lam[table.n]);
    for (let k = 0; k < path.count; k++) {
      assert.ok(path.chi[k + 1] > path.chi[k] && path.lam[k + 1] > path.lam[k], `χ and λ increase at ${k}`);
      assert.ok(Math.abs(path.phi[k + 1] - path.phi[k]) <= MAX_TURN, `φ at ${k}`);
      assert.ok(Math.abs(path.psi[k + 1] - path.psi[k]) <= MAX_TURN, `ψ at ${k}`);
    }
    const a2 = orbit.a * orbit.a;
    const b2 = orbit.b * orbit.b;
    for (let k = 0; k <= path.count; k++) {
      const onSpheroid = (path.x[k] ** 2 + path.y[k] ** 2) / (orbit.x + a2) + path.z[k] ** 2 / (orbit.x + b2);
      assert.ok(Math.abs(onSpheroid - 1) < 1e-12, `on the spheroid at ${k}`);
    }
  });

  test(`repeats the oscillation exactly as the physics does: ${name}`, () => {
    const { orbit, path } = prepare(args);
    // A trail of several oscillations either side of λ = 0, read back point by point against embed().
    const oscillations = 5;
    const positions = new Float32Array(6 * 20 * path.count);
    const colours = new Float32Array(positions.length);
    const psis = [];
    const paint = (psi, fraction, out, offset) => {
      psis[offset / 3] = psi;
      out[offset] = fraction;
    };
    const head = 2.5 * path.period;
    const segments = fill(path, head, oscillations * path.period, positions, colours, paint, [0, 0, 0, 0]);
    assert.ok(Math.abs(segments - (oscillations * path.count + 1)) <= 1, `${segments} segments`);
    // Segment s starts at point g = first + s, from oscillation m = −3 on (the tail is at λ = −2.5 periods).
    const first = locate(path, head - oscillations * path.period);
    for (let s = 1; s < segments; s++) {
      const g = first.m * path.count + first.k + s;
      const m = Math.floor(g / path.count);
      const k = g - m * path.count;
      const expected = embed(orbit, path.chi[k] + m * TWO_PI, path.phi[k] + m * path.phi[path.count]);
      for (let c = 0; c < 3; c++) {
        assert.ok(Math.abs(positions[6 * s + c] - expected[c]) < 1e-6, `${'xyz'[c]} at point ${g}: ${positions[6 * s + c]} vs ${expected[c]}`);
      }
      assert.ok(Math.abs(psis[2 * s - 1] - (path.psi[k] + m * path.psiAdvance)) < 1e-9, `ψ at point ${g}`); // painted as the end of segment s − 1
    }
  });

  test(`the panel's angles are the picture's, and run on smoothly through the axes: ${name}`, () => {
    const { orbit, path } = prepare(args);
    // The same stretch, as fill() draws it and as trace() visits it: the same points, in the same order.
    const head = 2.3 * path.period;
    const length = 5 * path.period;
    const positions = new Float32Array(6 * 20 * path.count);
    const photon = [0, 0, 0, 0];
    const segments = fill(path, head, length, positions, new Float32Array(positions.length), () => {}, photon);
    const visits = [];
    trace(path, head, length, (lam, theta, phi, psi) => visits.push([lam, theta, phi, psi]));
    assert.equal(visits.length, segments + 1);
    assert.ok(Math.abs(visits[0][0] - (head - length)) < 1e-9, 'from the tail');
    assert.ok(Math.abs(visits.at(-1)[0] - head) < 1e-9, 'to the photon');
    assert.ok(Math.abs(visits.at(-1)[3] - photon[3]) < 1e-9, 'with the photon\'s ψ');
    const across = Math.sqrt(orbit.x + orbit.a ** 2);
    const up = Math.sqrt(orbit.x + orbit.b ** 2);
    visits.forEach(([, theta, phi], s) => {
      // Signed θ and unwrapped φ put each point where the picture does, turned through the axes and all. The two
      // ends fall between points, where fill() interpolates the position and trace() the angles: a step's chord
      // and its arc, at most 10° apart.
      const at = 6 * Math.min(s, segments - 1) + (s === segments ? 3 : 0);
      const tolerance = s === 0 || s === segments ? 0.01 : 1e-6;
      const expected = [across * Math.sin(theta) * Math.cos(phi), across * Math.sin(theta) * Math.sin(phi), up * Math.cos(theta)];
      expected.forEach((value, c) => assert.ok(Math.abs(positions[at + c] - value) < tolerance, `${'xyz'[c]} at visit ${s}: ${positions[at + c]} vs ${value}`));
    });
    for (let s = 1; s < visits.length; s++) {
      const [lam0, theta0, phi0, psi0] = visits[s - 1];
      const [lam1, theta1, phi1, psi1] = visits[s];
      assert.ok(lam1 >= lam0, `λ increases at ${s}`);
      assert.ok(Math.abs(theta1 - theta0) < 0.1, `θ is continuous at ${s}: ${theta0} to ${theta1}`);
      assert.ok(Math.abs(phi1 - phi0) <= MAX_TURN + 1e-12 && Math.abs(psi1 - psi0) <= MAX_TURN + 1e-12, `φ and ψ at ${s}`);
    }
    const thetas = visits.map(([, theta]) => theta);
    const [low, high] = path.band;
    if (path.pole && path.equator) assert.ok(Math.max(...thetas) - Math.min(...thetas) > 4 * Math.PI, 'θ = χ/2 grows');
    else if (path.pole) assert.ok(Math.abs(Math.min(...thetas) + high) < 1e-3 && Math.abs(Math.max(...thetas) - high) < 1e-3, 'θ swings from −θ_max to θ_max');
    else if (path.equator) assert.ok(Math.abs(Math.min(...thetas) - low) < 1e-3 && Math.abs(Math.max(...thetas) - (Math.PI - low)) < 1e-3, 'θ swings from θ_min to π − θ_min');
    else assert.ok(Math.min(...thetas) >= low - 1e-9 && Math.max(...thetas) <= high + 1e-9, 'θ stays in its band');
  });
}

test('a trail with too many points is visited evenly, from its tail to the photon', () => {
  const { path } = prepare(CASES['near the extremal edge, φ winding fast']);
  const head = 1000;
  const length = 300;
  const all = [];
  trace(path, head, length, (lam) => all.push(lam));
  const some = [];
  trace(path, head, length, (lam) => some.push(lam), 20000);
  assert.ok(all.length > 100000, `${all.length} points`);
  assert.ok(some.length <= 20002 && some.length > 15000, `${some.length} visited`);
  assert.equal(some[0], all[0]);
  assert.equal(some.at(-1), all.at(-1));
  for (let k = 1; k < some.length; k++) assert.ok(some[k] >= some[k - 1]);
  // A trail within the budget is visited in full.
  const fig = prepare(CASES['Fig 4.10, near both axes']).path;
  const count = (most) => {
    let n = 0;
    trace(fig, head, length, () => n++, most);
    return n;
  };
  assert.equal(count(20000), count(Infinity));
});

test('the band is the orbit\'s', () => {
  const { orbit, path } = prepare(CASES['Fig 4.7']);
  assert.ok(Math.abs(Math.sin(path.band[0]) ** 2 - orbit.y1) < 1e-12);
  assert.ok(Math.abs(Math.cos(path.band[1]) ** 2 - orbit.w2) < 1e-12);
  assert.deepEqual(prepare(CASES['through the pole (Fig 4.4)']).path.band[0], 0);
  assert.ok(Math.abs(prepare(CASES['through the equator']).path.band[1] - Math.PI / 2) < 1e-12);
});

test('the trail runs from its tail to the photon, and fades from 0 to 1 along it', () => {
  const { orbit, path } = prepare(CASES['Fig 4.7']);
  const positions = new Float32Array(6 * 4096);
  const colours = new Float32Array(positions.length);
  const photon = [0, 0, 0, 0];
  const paint = (psi, fraction, out, offset) => out.set([fraction, psi, 0], offset);
  const head = 7.3;
  const length = 4.1;
  const segments = fill(path, head, length, positions, colours, paint, photon);
  const at = (lam) => {
    const chi = chiAt(path, lam);
    const m = Math.floor(chi / TWO_PI);
    const place = locate(path, lam);
    const phi = m * path.phi[path.count] + path.phi[place.k] + place.t * (path.phi[place.k + 1] - path.phi[place.k]);
    return embed(orbit, chi, phi);
  };
  const end = 6 * (segments - 1) + 3;
  at(head).forEach((value, c) => assert.ok(Math.abs(positions[end + c] - value) < 1e-3, 'the trail ends at the photon'));
  at(head - length).forEach((value, c) => assert.ok(Math.abs(positions[c] - value) < 1e-3, 'and starts at its tail'));
  assert.deepEqual([...photon.slice(0, 3)].map((v) => v.toFixed(5)), [...positions.slice(end, end + 3)].map((v) => v.toFixed(5)));
  assert.equal(colours[0], 0);
  assert.equal(colours[end], 1);
  for (let s = 0; s < segments; s++) {
    assert.ok(colours[6 * s] <= colours[6 * s + 3], 'the fade rises towards the photon');
    for (let c = 0; c < 3; c++) assert.equal(positions[6 * s + c], s === 0 ? positions[c] : positions[6 * s - 3 + c], 'segments join up');
  }
});

test('a trail too long for the buffer keeps its newest part', () => {
  const { path } = prepare(CASES['Fig 4.7']);
  const capacity = 1000;
  const positions = new Float32Array(6 * capacity);
  const colours = new Float32Array(positions.length);
  const photon = [0, 0, 0, 0];
  const segments = fill(path, 50, 300, positions, colours, () => {}, photon);
  assert.equal(segments, capacity);
  const end = 6 * (capacity - 1) + 3;
  for (let c = 0; c < 3; c++) assert.ok(Math.abs(positions[end + c] - photon[c]) < 1e-6);
});

test('χ and λ convert both ways, so the photon keeps its place when the orbit changes', () => {
  const { path } = prepare(CASES['Fig 4.10, near both axes']);
  for (const lam of [0, 0.37, 2.7, -5.2, 123.456]) assert.ok(Math.abs(lamAt(path, chiAt(path, lam)) - lam) < 1e-9, `λ = ${lam}`);
  assert.equal(chiAt(path, 0), 0);
  assert.ok(Math.abs(chiAt(path, path.period) - TWO_PI) < 1e-12);
  const other = prepare([0.6, 0.3, 1, -1]).path;
  const chi = chiAt(path, 10.4);
  assert.ok(Math.abs(chiAt(other, lamAt(other, chi)) - chi) < 1e-9);
});

test('the paper\'s colour wheel', () => {
  const close = (actual, expected) => actual.every((v, i) => Math.abs(v - expected[i]) < 1e-12);
  assert.ok(close(hue(0), [1, 0, 0]), 'red at 0');
  assert.ok(close(hue(1 / 6), [1, 1, 0]), 'yellow at π/3');
  assert.ok(close(hue(2 / 6), [0, 1, 0]), 'green at 2π/3');
  assert.ok(close(hue(3 / 6), [0, 1, 1]), 'cyan at π');
  assert.ok(close(hue(4 / 6), [0, 0, 1]), 'blue at 4π/3');
  assert.ok(close(hue(5 / 6), [1, 0, 1]), 'magenta at 5π/3');
  assert.equal(wheelIndex(0), 0);
  assert.equal(wheelIndex(TWO_PI), 0);
  assert.equal(wheelIndex(-TWO_PI / 4), (3 * WHEEL_SIZE) / 4);
  assert.equal(wheelIndex(TWO_PI * (1 - 1e-9)), 0);
  assert.equal(WHEEL.length, 3 * WHEEL_SIZE);
  assert.ok(Math.abs(WHEEL[3 * (WHEEL_SIZE / 2) + 1] - 1) < 1e-6, 'cyan is full green and blue, in linear light too');
});

test('fills a long trail quickly enough for every frame', (t) => {
  const { path } = prepare(CASES['Fig 4.10, near both axes']);
  const positions = new Float32Array(6 * 65536);
  const colours = new Float32Array(positions.length);
  const photon = [0, 0, 0, 0];
  const paint = (psi, fraction, out, offset) => {
    const k = 3 * wheelIndex(psi);
    out[offset] = WHEEL[k] * fraction;
    out[offset + 1] = WHEEL[k + 1] * fraction;
    out[offset + 2] = WHEEL[k + 2] * fraction;
  };
  let segments = 0;
  const repeats = 20;
  const began = performance.now();
  for (let r = 0; r < repeats; r++) segments = fill(path, 1000 + r, 300, positions, colours, paint, photon);
  const each = (performance.now() - began) / repeats;
  t.diagnostic(`${segments} segments (${path.count} per oscillation): ${each.toFixed(2)} ms a frame`);
  assert.ok(each < 16);

  // And the panel's walk along the same trail.
  let sum = 0;
  const visit = (lam, theta, phi, psi) => {
    sum += theta + phi + psi;
  };
  const start = performance.now();
  for (let r = 0; r < repeats; r++) trace(path, 1000 + r, 300, visit);
  const traced = (performance.now() - start) / repeats;
  t.diagnostic(`trace: ${traced.toFixed(2)} ms a frame (${Number.isFinite(sum)})`);
  assert.ok(traced < 16);
});
