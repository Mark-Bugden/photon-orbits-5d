// The demo's entry point. Vite bundles it with the kit, three.js and the CSS into demo.js and demo.css.
//
// The controls, top to bottom: the presets, the spins, the momentum pad with its two sliders, and the trail's
// length. A preset sets all of them, the photon's place and the camera. The stage holds the 3D view and, beside
// it, the side panel, which plots the same trail.

import manifest from 'virtual:manifest';
import { mount, select, slider } from '@mark-bugden/demo-kit';
import '@mark-bugden/demo-kit/tokens.scoped.css';
import '@mark-bugden/demo-kit/frame.css';
import './demo.css';

import { extremal } from './momenta.js';
import { pad } from './pad.js';
import { createLegend } from './panel/legend.js';
import { createPanel } from './panel/panel.js';
import { frequencies, horizon, solveOrbit, tabulate } from './physics/index.js';
import { OPENING, PRESETS } from './presets.js';
import { createView } from './scene/index.js';
import { spins } from './spins.js';
import { chiAt, lamAt, refine } from './trail.js';

const RATE = 2; // λ per second: about one oscillation every 1.5 s
const TRAIL = [3, 300]; // the trail's length in λ, from about one oscillation to about 100

const EXTREMAL = 'Extremal: a + b = 1. The inner and outer horizons have merged, and the hole spins as fast as a black hole of its mass can.';
const HELD = 'a + b stops at 1, where the hole is extremal: its inner and outer horizons have merged. Beyond that, it would have no horizon.';
const EDGE = (name) => `${name} stops at the edge of the shaded region: past it, no spherical orbit has these momenta.`;

const signed = (value) => `${value < 0 ? '−' : '+'}${Math.abs(value).toFixed(2)}`;
const degrees = (radians) => ((180 / Math.PI) * radians).toFixed(1);

const controls = {
  preset: select({
    label: 'Preset',
    options: [{ label: 'Custom', value: null }, ...PRESETS.map((preset) => ({ label: preset.name, value: preset }))],
    value: PRESETS.find((preset) => preset.name === OPENING.name),
  }),
  spins: spins({ value: [OPENING.a, OPENING.b] }),
  momenta: pad({ spins: [OPENING.a, OPENING.b], value: [OPENING.Phi, OPENING.Psi] }),
  trail: slider({
    label: 'Trail length',
    min: Math.log10(TRAIL[0]),
    max: Math.log10(TRAIL[1]),
    step: 0.01,
    format: (value) => `λ ${Math.round(10 ** value)}`,
  }),
};

let view = null;
let legend = null;
let panel = null;
const viewSize = { width: 0, height: 0 }; // of the 3D view, beside the panel
let current = null; // { key, orbit, table, path, … }, or { key, orbit, … } with no path when there is no orbit
let lam = 0; // the photon's λ
let chi = 0; // and its χ, which keeps its place in the oscillation when the orbit changes
let length = 0; // the trail's, in λ: a preset's exactly, until the slider (which rounds it) moves
let start = OPENING; // the preset last chosen: where Reset puts the photon back
let psi0 = OPENING.psi0; // ψ at χ = 0, which sets the colours and the panel's ψ
let phi0 = 0; // φ at χ = 0, for the panel's φ (the 3D view's is in its camera)
let lam0 = 0; // the demo's λ at the paper's λ = 0, for the panel's λ axis
let aim = OPENING.view; // a viewpoint to turn the camera to once the spheroid is drawn
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
    chi = start.span[1];
    if (current?.path) lam = lamAt(current.path, chi);
  },
  render: draw,
});

// Any change of the user's own makes the state custom, and clears the last note.
function custom(text = '') {
  controls.preset.value = null;
  note = text;
}
controls.spins.onChange(([a, b], { held }) => {
  controls.momenta.setSpins(a, b);
  custom(held ? HELD : extremal(a, b) ? EXTREMAL : '');
});
controls.momenta.onChange((value, { stopped, axis }) => custom(stopped ? EDGE(axis === 'Phi' ? 'Φ' : 'Ψ') : ''));
controls.trail.onChange((value) => {
  length = 10 ** value;
  custom();
});
controls.preset.onChange((preset) => preset && apply(preset));

function apply(preset) {
  const { a, b, Phi, Psi, span } = preset;
  controls.spins.value = [a, b];
  controls.momenta.setSpins(a, b);
  controls.momenta.value = [Phi, Psi];
  current = solve(a, b, Phi, Psi);
  view?.setShape(current.orbitAxes, current.horizonAxes);
  chi = span[1];
  lam = lamAt(current.path, chi);
  length = lam - lamAt(current.path, span[0]);
  controls.trail.value = Math.log10(length);
  start = preset;
  psi0 = preset.psi0;
  phi0 = preset.phi0 ?? 0;
  lam0 = preset.lam0 ?? 0;
  aim = preset.view;
  note = extremal(a, b) ? EXTREMAL : '';
}

// The orbit for these spins and momenta, its table and drawn path, and the semi-axes of its spheroid and the horizon's.
function solve(a, b, Phi, Psi) {
  const key = [a, b, Phi, Psi].join();
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

function draw({ pixelRatio, values, delta, theme, token }) {
  if (!view) return;
  view.setSize(viewSize.width, viewSize.height, pixelRatio);
  if (theme !== painted) {
    view.setBackground(token('--color-surface'));
    painted = theme;
  }
  const [a, b] = values.spins;
  const [Phi, Psi] = values.momenta;
  const key = [a, b, Phi, Psi].join();
  if (key !== current?.key) {
    current = solve(a, b, Phi, Psi);
    view.setShape(current.orbitAxes, current.horizonAxes);
    if (current.path) lam = lamAt(current.path, chi);
  }
  if (aim) {
    view.setShape(current.orbitAxes, current.horizonAxes);
    view.fit(aim);
    aim = null;
  }
  if (current.path) {
    lam += RATE * delta;
    chi = chiAt(current.path, lam);
    legend.set(view.showTrail(current.path, lam, length, psi0));
    panel.show(current.path, lam, length, lam0, phi0, psi0);
  } else {
    view.hideTrail();
    panel.hide();
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

// The stage: the 3D view, with the colour legend over it, and the side panel beside it (below it on phones).
function build() {
  const layout = document.createElement('div');
  layout.className = 'po-layout';
  const box = document.createElement('div');
  box.className = 'po-view';
  layout.append(box);
  frame.stage.append(layout);
  view = createView(box, {
    onChange: () => frame.redraw(),
    reducedMotion: globalThis.matchMedia?.('(prefers-reduced-motion: reduce)').matches ?? false,
  });
  legend = createLegend(box);
  panel = createPanel({ onChange: () => frame.redraw() });
  layout.append(panel.el);
  new ResizeObserver(([entry]) => {
    viewSize.width = entry.contentRect.width;
    viewSize.height = entry.contentRect.height;
    frame.redraw();
  }).observe(box);
}

apply(OPENING);
try {
  build();
} catch (error) {
  console.error(error);
  view = null;
  frame.pause();
  for (const control of Object.values(controls)) control.disabled = true;
  frame.stage.replaceChildren();
  showStill();
  frame.status('This demo draws with WebGL, which this browser could not start.');
}
