import assert from 'node:assert/strict';
import { test } from 'node:test';

import { WHEEL, WHEEL_SIZE, hue, wheelIndex } from '../src/colour.js';
import { embed, solveOrbit, tabulate } from '../src/physics/index.js';
import { BASE_STEPS, MAX_TURN, chiAt, fill, lamAt, locate, refine } from '../src/trail.js';

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
      assert.ok(Math.abs(psis[2 * s - 1] - (path.psi[k] + m * path.advance)) < 1e-9, `ψ at point ${g}`); // painted as the end of segment s − 1
    }
  });
}

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
});
