import assert from 'node:assert/strict';
import { test } from 'node:test';

import { angleTicks, linearTicks, numberLabel, piLabel, tickSet } from '../src/panel/ticks.js';

const PI = Math.PI;
const values = (ticks) => [...ticks.values.subarray(0, ticks.count)];
const labels = (ticks) => ticks.labels.slice(0, ticks.count);

test('angles are labelled as multiples of π', () => {
  assert.deepEqual([0, PI / 4, PI / 2, (3 * PI) / 4, PI, -2 * PI, 40 * PI, -PI / 2].map(piLabel), ['0', 'π/4', 'π/2', '3π/4', 'π', '−2π', '40π', '−π/2']);
});

test('numbers get the decimals their step needs, and a true minus sign', () => {
  assert.equal(numberLabel(2, 1), '2');
  assert.equal(numberLabel(2.5, 0.5), '2.5');
  assert.equal(numberLabel(-1, 1), '−1');
  assert.equal(numberLabel(-1e-17, 0.2), '0.0');
  assert.equal(numberLabel(150, 50), '150');
});

test('λ\'s ticks are round numbers, including those at the ends', () => {
  // Fig 4.10's window, λ from 1 to 5 to within its fitted start.
  const ticks = linearTicks(1.0000043, 4.9999946, 6, tickSet());
  assert.deepEqual(values(ticks), [1, 2, 3, 4, 5]);
  assert.deepEqual(labels(ticks), ['1', '2', '3', '4', '5']);
  assert.deepEqual(values(linearTicks(63.4, 71.6, 5, tickSet())), [64, 66, 68, 70]);
  assert.deepEqual(values(linearTicks(63.4, 71.6, 4, tickSet())), [65, 70]);
  assert.deepEqual(values(linearTicks(-0.08, 17.48, 5, tickSet())), [0, 5, 10, 15]);
});

test('θ\'s ticks are fractions of π: 0, π/4, π/2, or wider steps on a longer range', () => {
  assert.deepEqual(labels(angleTicks(0, PI / 2, 4, tickSet())), ['0', 'π/4', 'π/2']);
  assert.deepEqual(labels(angleTicks(-PI / 2, PI / 2, 5, tickSet())), ['−π/2', '−π/4', '0', 'π/4', 'π/2']);
  assert.deepEqual(labels(angleTicks(-PI / 2, PI / 2, 4, tickSet())), ['−π/2', '0', 'π/2']);
  assert.deepEqual(labels(angleTicks(0, 30 * PI, 4, tickSet())), ['0', '10π', '20π', '30π']);
});

test('a tick keeps its label while it stays put, so drawing allocates nothing', () => {
  const ticks = tickSet();
  linearTicks(63.4, 71.6, 5, ticks);
  const before = labels(ticks);
  linearTicks(63.5, 71.7, 5, ticks);
  labels(ticks).forEach((label, k) => assert.equal(label, before[k]));
  linearTicks(64.5, 72.7, 5, ticks);
  assert.deepEqual(labels(ticks), ['66', '68', '70', '72']);
});
