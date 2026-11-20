// The demo's entry point. Vite bundles it with the kit, three.js and the CSS into demo.js and demo.css.
//
// For now the momenta are set with plain sliders, and most pairs have no orbit; step 3 of the plan replaces them
// with a pad that shows where orbits exist.

import manifest from 'virtual:manifest';
import { mount, slider } from '@mark-bugden/demo-kit';
import '@mark-bugden/demo-kit/tokens.scoped.css';
import '@mark-bugden/demo-kit/frame.css';
import './demo.css';

import { frequencies, horizon, solveOrbit, tabulate } from './physics/index.js';
import { createView } from './scene/index.js';
import { chiAt, lamAt, refine } from './trail.js';

// Fig 4.10, until step 3 of the plan chooses the opening state from the presets.
const OPENING = { a: 0.9, b: 0.09, Phi: 0.1, Psi: 0.01, trail: 40 };
const RATE = 2; // λ per second: about one oscillation every 1.5 s
const TRAIL = [3, 300]; // the trail's length in λ, from about one oscillation to about 100

const minus = (text) => text.replace('-', '−');
const fixed = (digits) => (value) => minus(value.toFixed(digits));
const signed = (value) => `${value < 0 ? '−' : '+'}${Math.abs(value).toFixed(2)}`;
const degrees = (radians) => ((180 / Math.PI) * radians).toFixed(1);

const controls = {
  a: slider({ label: 'Spin a', min: 0, max: 0.99, step: 0.01, value: OPENING.a, format: fixed(2) }),
  b: slider({ label: 'Spin b', min: 0, max: 0.99, step: 0.01, value: OPENING.b, format: fixed(2) }),
  Phi: slider({ label: 'Momentum Φ', min: -3, max: 3, step: 0.01, value: OPENING.Phi, format: fixed(2) }),
  Psi: slider({ label: 'Momentum Ψ', min: -3, max: 3, step: 0.01, value: OPENING.Psi, format: fixed(2) }),
  trail: slider({
    label: 'Trail length',
    min: Math.log10(TRAIL[0]),
    max: Math.log10(TRAIL[1]),
    step: 0.01,
    value: Math.log10(OPENING.trail),
    format: (value) => `λ ${Math.round(10 ** value)}`,
  }),
};

let view = null;
let current = null; // { key, orbit, table, path }, or { key, orbit } when there is no orbit
let lam = 0; // the photon's λ
let chi = 0; // and its χ, which keeps its place in the oscillation when the orbit changes
let note = ''; // shown in the status line instead of the readout, until the next change
let shown = null;
let painted = null; // the theme the view's background was last painted for

const frame = mount({
  el: document.querySelector(`[data-demo="${manifest.slug}"]`),
  manifest,
  controls,
  stage: 'div',
  animate: true,
  reset() {
    lam = 0;
    chi = 0;
  },
  render: draw,
});

// a + b ≤ 1: beyond that there is no horizon. The slider that would pass it stops at the limit instead.
function holdSpins(moved, other) {
  return () => {
    note = '';
    if (controls.a.value + controls.b.value <= 1 + 1e-9) return;
    controls[moved].value = Math.round(100 * (1 - controls[other].value)) / 100;
    note = 'a + b stops at 1: beyond that, the black hole would have no horizon.';
  };
}
controls.a.onChange(holdSpins('a', 'b'));
controls.b.onChange(holdSpins('b', 'a'));
for (const name of ['Phi', 'Psi', 'trail']) controls[name].onChange(() => (note = ''));

// The orbit for the controls, its table and drawn path, and the semi-axes of its spheroid and the horizon's.
function solve({ a, b, Phi, Psi }, key) {
  const xPlus = horizon(a, b);
  const horizonAxes = [Math.sqrt(xPlus + a * a), Math.sqrt(xPlus + b * b)];
  const orbit = solveOrbit(a, b, Phi, Psi);
  if (!orbit.ok) return { key, orbit, horizonAxes, orbitAxes: null };
  const table = tabulate(orbit);
  return {
    key, orbit, table, path: refine(orbit, table), horizonAxes,
    orbitAxes: [Math.sqrt(orbit.x + a * a), Math.sqrt(orbit.x + b * b)],
  };
}

function readout({ orbit, table }) {
  if (!table) return 'No spherical orbit for these momenta.';
  const [turnsPhi, turnsPsi] = frequencies(orbit, table);
  const low = degrees(Math.asin(Math.sqrt(orbit.y1)));
  const high = degrees(Math.acos(Math.sqrt(orbit.w2)));
  return `r = ${Math.sqrt(orbit.x).toFixed(2)} · θ from ${low}° to ${high}° · per swing: φ ${signed(turnsPhi)} turns, ψ ${signed(turnsPsi)} turns`;
}

function draw({ width, height, pixelRatio, values, delta, theme, token }) {
  if (!view) return;
  view.setSize(width, height, pixelRatio);
  if (theme !== painted) {
    view.setBackground(token('--color-surface'));
    painted = theme;
  }
  const key = [values.a, values.b, values.Phi, values.Psi].join();
  if (key !== current?.key) {
    const first = current === null;
    current = solve(values, key);
    view.setShape(current.orbitAxes, current.horizonAxes);
    if (current.path) lam = lamAt(current.path, chi);
    if (first) view.fit();
  }
  if (current.path) {
    lam += RATE * delta;
    chi = chiAt(current.path, lam);
    view.showTrail(current.path, lam, 10 ** values.trail);
  } else {
    view.hideTrail();
  }
  const status = note || readout(current);
  if (status !== shown) frame.status((shown = status));
  view.render();
}

// static.webp sits beside demo.js in the bundle; `npm run dev` serves it from the page's folder.
function showStill() {
  const base = import.meta.env.DEV ? document.baseURI : import.meta.url;
  const still = document.createElement('img');
  Object.assign(still, {
    src: new URL('static.webp', base).href,
    alt: 'Title card reading \'Plotting the Photons: spherical photon orbits around a five-dimensional rotating black hole\', beside a drawing of a black disc circled by a ring.',
    width: 1200,
    height: 675,
  });
  frame.stage.append(still);
}

try {
  view = createView(frame.stage, {
    onChange: () => frame.redraw(),
    reducedMotion: globalThis.matchMedia?.('(prefers-reduced-motion: reduce)').matches ?? false,
  });
} catch (error) {
  console.error(error);
  frame.pause();
  for (const control of Object.values(controls)) control.disabled = true;
  showStill();
  frame.status('This demo draws with WebGL, which this browser could not start.');
}
