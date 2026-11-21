import assert from 'node:assert/strict';
import { test } from 'node:test';

import { GAP, allowed, clampAlong, extremal, locate, place, reach, toRegion } from '../src/momenta.js';
import { directions, edge, orbitExists, solveOrbit } from '../src/physics/index.js';
import { OPENING, PRESETS } from '../src/presets.js';

const HOLES = [[0, 0], [0.5, 0], [0, 0.99], [0.6, 0.3], [0.5, 0.25], [0.9, 0.09], [0.4, 0.4]];
const EXTREMAL = [[0.5, 0.5], [0.7, 0.3], [0.9, 0.1], [0.99, 0.01]];
const RAYS = directions(24);

test('away from extremality, the pad reaches the edge of the region', () => {
  for (const [a, b] of HOLES) {
    assert.equal(extremal(a, b), false);
    for (const [dPhi, dPsi] of RAYS) {
      const s = reach(a, b, dPhi, dPsi);
      assert.equal(s, edge(a, b, dPhi, dPsi));
      assert.ok(allowed(a, b, s * dPhi, s * dPsi));
      assert.ok(!orbitExists(a, b, (s + 1e-9) * dPhi, (s + 1e-9) * dPsi), `just past the edge at ${a}, ${b}`);
    }
  }
});

test('at an extremal hole, the dot stops GAP short of the edge Φ + Ψ = −1', () => {
  for (const [a, b] of EXTREMAL) {
    assert.ok(extremal(a, b));
    for (const [dPhi, dPsi] of RAYS) {
      const s = reach(a, b, dPhi, dPsi);
      const [Phi, Psi] = [s * dPhi, s * dPsi];
      assert.ok(allowed(a, b, Phi, Psi), `allowed at the stop: ${a}, ${b}, ${dPhi}, ${dPsi}`);
      if (dPhi < 0 && dPsi < 0) {
        assert.ok(Math.abs(1 + Phi + Psi - GAP) < 1e-12, `GAP from the edge: ${1 + Phi + Psi}`);
        const orbit = solveOrbit(a, b, Phi, Psi);
        assert.ok(orbit.ok && orbit.x - orbit.xPlus > 5e-3, 'the orbit stays clear of the horizon');
      }
    }
  }
});

test('a slider stops at the edge, and keeps the other momentum exactly', () => {
  for (const [a, b] of [...HOLES, ...EXTREMAL]) {
    for (const from of [[0, 0], [0.3, -0.2], [0, 0.4], [-0.2, 0]]) {
      if (!allowed(a, b, ...from)) continue;
      for (const to of [[3, from[1]], [-3, from[1]], [from[0], 3], [from[0], -3]]) {
        const [Phi, Psi] = clampAlong(a, b, from, to);
        assert.ok(allowed(a, b, Phi, Psi));
        if (to[1] === from[1]) assert.equal(Psi, from[1]);
        else assert.equal(Phi, from[0]);
        const further = [Phi + 1e-9 * Math.sign(to[0] - from[0]), Psi + 1e-9 * Math.sign(to[1] - from[1])];
        assert.ok(!allowed(a, b, ...further), `stops at the edge: ${a}, ${b}, ${from} → ${to}`);
      }
      assert.deepEqual(clampAlong(a, b, from, [0.01, 0.02]), [0.01, 0.02], 'inside, it goes where it is told');
    }
  }
});

test('a point dragged outside goes to the edge along its ray', () => {
  for (const [a, b] of [...HOLES, ...EXTREMAL]) {
    for (const [dPhi, dPsi] of RAYS) {
      const [Phi, Psi] = toRegion(a, b, 3.5 * dPhi, 3.5 * dPsi);
      assert.ok(allowed(a, b, Phi, Psi));
      assert.ok(Math.abs(Phi * dPsi - Psi * dPhi) < 1e-12, 'on the ray');
    }
  }
});

test('the dot keeps its direction and its fraction of the way to the edge when the spins change', () => {
  const start = [0.6, 0.3];
  for (const point of [[0.4, -0.3], [-0.5, -0.4], [0, 0.7], [-0.9, 0], [0, 0], [0.01, 0.002]]) {
    const where = locate(...start, ...point);
    for (const [a, b] of [...HOLES, ...EXTREMAL]) {
      const moved = place(a, b, where);
      assert.ok(allowed(a, b, ...moved), `allowed at ${a}, ${b}`);
      const again = locate(a, b, ...moved);
      assert.ok(Math.abs(again.fraction - where.fraction) < 1e-12);
      if (point[0] === 0) assert.equal(moved[0], 0, 'Φ = 0 stays exactly 0');
      if (point[1] === 0) assert.equal(moved[1], 0, 'Ψ = 0 stays exactly 0');
    }
  }
});

test('every fraction up to the edge is allowed, so the region is star-shaped as the pad assumes', () => {
  for (const [a, b] of [...HOLES, ...EXTREMAL]) {
    for (const direction of RAYS) {
      for (const fraction of [0.25, 0.5, 0.9, 0.999, 1]) {
        assert.ok(allowed(a, b, ...place(a, b, { direction, fraction })), `${a}, ${b}, ${direction}, ${fraction}`);
      }
    }
  }
});

test('every preset has an orbit, no preset names its figure, and the opening state is zero momenta', () => {
  assert.equal(OPENING.name, 'Zero momenta');
  assert.deepEqual([OPENING.Phi, OPENING.Psi], [0, 0]);
  assert.ok(PRESETS.includes(OPENING));
  for (const preset of PRESETS) {
    assert.doesNotMatch(preset.name, /Fig/, preset.name);
    const { a, b, Phi, Psi, span, view } = preset;
    assert.ok(a + b <= 1 + 1e-12, preset.name);
    assert.ok(allowed(a, b, Phi, Psi), preset.name);
    assert.ok(span[0] >= 0 && span[1] > span[0], preset.name);
    assert.equal(view.length, 2);
  }
  const circular = solveOrbit(0.5, 0.25, -0.929241, 0.903667);
  assert.ok(Math.abs(Math.asin(Math.sqrt(circular.y1)) - Math.PI / 4) < 1e-3, 'Fig 4.12 is at θ = 45°');
  assert.ok(Math.abs(Math.acos(Math.sqrt(circular.w2)) - Math.PI / 4) < 1e-3);
});
