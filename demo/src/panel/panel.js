// The side panel, beside the 3D view: two tabs over one canvas. The first plots θ, φ and ψ against λ over the
// trail, as the paper's Figs 4.9 and 4.11 do; the second puts the trail on the (φ, ψ) torus.

import { getTheme, token } from '@mark-bugden/demo-kit';

import { createStrip } from './strip.js';
import { createTorus } from './torus.js';

const TABS = [
  {
    label: 'θ, φ, ψ against λ',
    description: 'θ, φ and ψ against λ over the trail, with the photon at the right-hand end. φ and ψ are unwrapped.',
  },
  {
    label: '(φ, ψ) torus',
    description: 'The trail on the (φ, ψ) torus: φ across and ψ up, each from 0 to 2π, wrapping round at the edges.',
  },
];
const FONT = '500 12px Inter, system-ui, sans-serif';
let panels = 0;

/** { el, show(path, head, length, lam0, phi0, psi0), hide() }. onChange() is called when the panel needs drawing. */
export function createPanel({ onChange }) {
  const el = document.createElement('section');
  el.className = 'po-panel';
  el.setAttribute('aria-label', 'The orbit\'s angles');
  const list = document.createElement('div');
  list.className = 'po-tabs';
  list.setAttribute('role', 'tablist');
  list.setAttribute('aria-label', 'Plot');
  const body = document.createElement('div');
  body.className = 'po-panel-body';
  body.id = `po-panel-${++panels}`;
  body.setAttribute('role', 'tabpanel');
  const canvas = document.createElement('canvas');
  canvas.className = 'po-panel-canvas';
  canvas.setAttribute('role', 'img');
  body.append(canvas);
  const tabs = TABS.map(({ label }, index) => {
    const tab = document.createElement('button');
    tab.type = 'button';
    tab.className = 'po-tab';
    tab.id = `${body.id}-tab-${index}`;
    tab.textContent = label;
    tab.setAttribute('role', 'tab');
    tab.setAttribute('aria-controls', body.id);
    tab.addEventListener('click', () => select(index));
    return tab;
  });
  list.append(...tabs);
  el.append(list, body);

  const views = [createStrip(), createTorus()];
  let active = 0;
  function select(index, focus = false) {
    active = index;
    tabs.forEach((tab, k) => {
      tab.setAttribute('aria-selected', String(k === index));
      tab.tabIndex = k === index ? 0 : -1;
    });
    body.setAttribute('aria-labelledby', tabs[index].id);
    canvas.setAttribute('aria-label', TABS[index].description);
    if (focus) tabs[index].focus();
    onChange();
  }
  // The arrow keys move between the tabs, as in any tab list.
  list.addEventListener('keydown', (event) => {
    const next = { ArrowLeft: active - 1, ArrowRight: active + 1, Home: 0, End: tabs.length - 1 }[event.key];
    if (next === undefined) return;
    event.preventDefault();
    select((next + tabs.length) % tabs.length, true);
  });
  select(0);

  // The canvas matches its box, at up to two device pixels a CSS pixel.
  let width = 0;
  let height = 0;
  let ratio = 1;
  const context = canvas.getContext('2d');
  new ResizeObserver(([entry]) => {
    width = entry.contentRect.width;
    height = entry.contentRect.height;
    ratio = Math.min(globalThis.devicePixelRatio || 1, 2);
    canvas.width = Math.max(1, Math.round(width * ratio));
    canvas.height = Math.max(1, Math.round(height * ratio));
    onChange();
  }).observe(canvas);

  // The theme's colours, read again when the theme changes.
  const colours = { bg: '', text: '', muted: '', border: '', curves: ['', '', ''] };
  let theme = null;
  function prepare() {
    if (getTheme() !== theme) {
      theme = getTheme();
      colours.bg = token(el, '--color-bg');
      colours.text = token(el, '--color-text');
      colours.muted = token(el, '--color-muted');
      colours.border = token(el, '--color-border');
      // The paper's colours for θ, φ and ψ: orange, green and red.
      colours.curves[0] = token(el, '--chart-2');
      colours.curves[1] = token(el, '--chart-3');
      colours.curves[2] = token(el, '--chart-8');
    }
    context.setTransform(ratio, 0, 0, ratio, 0, 0);
    context.clearRect(0, 0, width, height);
  }

  return {
    el,
    /** Draws the trail from λ = head − length to head; the axes read λ − lam0, φ + phi0 and ψ + psi0. */
    show(path, head, length, lam0, phi0, psi0) {
      if (!width || !height) return;
      prepare();
      views[active].draw(context, width, height, ratio, path, head, length, lam0, phi0, psi0, colours);
    },
    /** Shows that there is no orbit to plot. */
    hide() {
      if (!width || !height) return;
      prepare();
      context.font = FONT;
      context.fillStyle = colours.muted;
      context.textAlign = 'center';
      context.textBaseline = 'middle';
      context.fillText('No orbit to plot', 0.5 * width, 0.5 * height);
    },
  };
}
