import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { test } from 'node:test';

import { boundary, edge, embed, frequencies, horizons, node, solveOrbit, tabulate, vertices } from '../src/physics/index.js';

// Results of the reference model, src/photon_orbits/, written by src/photon_orbits/fixtures.py.
const fixtures = JSON.parse(await readFile(new URL('../../tests/fixtures/orbits.json', import.meta.url), 'utf8'));

const TOLERANCE = 1e-9;
let worst = 0;
function close(actual, expected, label) {
  const difference = Math.abs(actual - expected) / Math.max(1, Math.abs(expected));
  worst = Math.max(worst, difference);
  assert.ok(difference <= TOLERANCE, `${label}: ${actual} vs ${expected}`);
}

test('reproduces the horizons', () => {
  for (const { a, b, xPlus, e } of fixtures.horizons) {
    const found = horizons(a, b);
    if (xPlus === null) assert.equal(found, null, `a=${a} b=${b}`);
    else {
      close(found.xPlus, xPlus, `x₊ at a=${a} b=${b}`);
      close(found.e, e, `e at a=${a} b=${b}`);
    }
  }
});

for (const fixture of fixtures.orbits) {
  const { a, b, Phi, Psi } = fixture;
  test(`reproduces the orbit at a=${a} b=${b} Phi=${Phi} Psi=${Psi}`, () => {
    const orbit = solveOrbit(a, b, Phi, Psi);
    assert.equal(orbit.ok, true);
    for (const [key, value] of Object.entries(fixture.orbit)) close(orbit[key], value, key);

    const table = tabulate(orbit);
    assert.equal(table.n, fixture.table.n);
    for (const key of ['sigma', 'lam', 'phi', 'psi']) {
      fixture.table[key].forEach((value, k) => close(table[key][k * fixture.table.every], value, `${key}[${k * fixture.table.every}]`));
    }
    frequencies(orbit, table).forEach((value, k) => close(value, fixture.frequencies[k], `ν[${k}]`));

    for (const point of fixture.points) {
      const at = node(table, point.i);
      for (const key of ['chi', 'sigma', 'lam', 'phi', 'psi']) close(at[key], point[key], `${key} at node ${point.i}`);
      embed(orbit, at.chi, at.phi).forEach((value, k) => close(value, point['XYZ'[k]], `${'XYZ'[k]} at node ${point.i}`));
    }
  });
}

test('reproduces the momenta with no orbit, and why', () => {
  for (const { a, b, Phi, Psi, reason, x, K } of fixtures.noOrbits) {
    const found = solveOrbit(a, b, Phi, Psi);
    assert.equal(found.ok, false);
    assert.equal(found.reason, reason, `a=${a} b=${b} Phi=${Phi} Psi=${Psi}`);
    if (x === null) assert.equal(found.x, null);
    else {
      close(found.x, x, 'x');
      close(found.K, K, 'K');
    }
  }
});

for (const region of fixtures.regions) {
  test(`reproduces the region's corners and edge for a=${region.a} b=${region.b}`, () => {
    vertices(region.a, region.b).forEach((corner, k) => corner.forEach((value, m) => close(value, region.vertices[k][m], `corner ${k}`)));
    region.directions.forEach(([dPhi, dPsi], k) => close(edge(region.a, region.b, dPhi, dPsi), region.edges[k], `edge along ray ${k}`));
  });
}

test('largest difference from the reference', (t) => {
  t.diagnostic(`largest relative difference from the reference: ${worst.toExponential(1)}`);
});

test('is fast enough to recompute on every change', (t) => {
  const repeats = 200;
  let began = performance.now();
  for (let r = 0; r < repeats; r++) tabulate(solveOrbit(0.9, 0.09, 0.1 + r * 1e-4, 0.01));
  const orbit = (performance.now() - began) / repeats;
  began = performance.now();
  for (let r = 0; r < 5; r++) boundary(0.6 + r * 1e-3, 0.3);
  const region = (performance.now() - began) / 5;
  t.diagnostic(`an orbit and its table: ${orbit.toFixed(3)} ms; the region's edge on 96 rays: ${region.toFixed(1)} ms`);
  assert.ok(orbit < 5 && region < 100);
});
