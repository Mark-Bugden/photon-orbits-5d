// The colour legend: the paper's wheel as a ring in a corner of the 3D view, labelled 0, π/2, π and 3π/2, with a
// needle at the photon's ψ. ψ grows counter-clockwise from 0 on the right, as in the paper's Fig 4.1.
//
// The wheel is Mathematica's Hue at full saturation and value, which is linear in sRGB between red, yellow, green,
// cyan, blue and magenta, so a conic gradient through those six draws it exactly.
//
// The wheel is the b-plane, so the spin b's arrow goes round it: the hole turns towards −ψ (docs/equations.md,
// "Signs"), which is clockwise here. It is drawn as the 3D view's arrow for a is, pale with a dark outline.

import { SURFACE } from '../scene/surfaces.js';

const LABELS = [['0', 0], ['π/2', 90], ['π', 180], ['3π/2', 270]]; // and their angles, counter-clockwise

// The arrow, in the SVG's units, in which the ring's outer radius is 28.
const REACH = 40; // the SVG's half-width
const RADIUS = 34.5; // the arc's, just clear of the needle
const END = (230 * Math.PI) / 180; // where the arc meets the head: lower left
const ARC = (285 * Math.PI) / 180;
const HEAD = [4, 7]; // the head's half-width and length

// The point at angle ψ on the arc, with y down.
const at = (angle, r = RADIUS) => [r * Math.cos(angle), -r * Math.sin(angle)];
const xy = ([x, y]) => `${x.toFixed(2)} ${y.toFixed(2)}`;

function arrow() {
  const [x, y] = at(END);
  const along = [Math.sin(END), Math.cos(END)]; // clockwise, on screen
  const out = [Math.cos(END), -Math.sin(END)];
  const arc = `M${xy(at(END + ARC))} A${RADIUS} ${RADIUS} 0 1 1 ${xy([x, y])}`;
  const head = [
    [x + HEAD[1] * along[0], y + HEAD[1] * along[1]],
    [x + HEAD[0] * out[0], y + HEAD[0] * out[1]],
    [x - HEAD[0] * out[0], y - HEAD[0] * out[1]],
  ].map(xy).join(' ');
  return `<svg class="po-legend-spin" viewBox="${-REACH} ${-REACH} ${2 * REACH} ${2 * REACH}" fill="none" stroke-linecap="round" stroke-linejoin="round">
    <path d="${arc}" stroke="${SURFACE.backdrop}" stroke-width="4.5"/>
    <polygon points="${head}" fill="${SURFACE.ring}" stroke="${SURFACE.backdrop}" stroke-width="2.3" paint-order="stroke"/>
    <path d="${arc}" stroke="${SURFACE.ring}" stroke-width="2.2"/>
  </svg>`;
}

export function createLegend(stage) {
  const el = document.createElement('div');
  el.className = 'po-legend';
  el.setAttribute('aria-hidden', 'true');
  const ring = document.createElement('div');
  ring.className = 'po-legend-ring';
  const needle = document.createElement('div');
  needle.className = 'po-legend-needle';
  const centre = document.createElement('span');
  centre.className = 'po-legend-centre';
  centre.textContent = 'ψ';
  el.append(ring, needle, centre);
  for (const [text, angle] of LABELS) {
    const label = document.createElement('span');
    label.className = 'po-legend-label';
    label.textContent = text;
    label.style.setProperty('--angle', `${-angle}deg`);
    el.append(label);
  }
  el.insertAdjacentHTML('beforeend', arrow());
  const spin = el.lastElementChild;
  stage.append(el);

  let shown = NaN;
  return {
    el,
    /** Points the needle at ψ. */
    set(psi) {
      const degrees = Math.round(((psi * 180) / Math.PI) % 360);
      if (degrees === shown) return;
      shown = degrees;
      needle.style.transform = `rotate(${-degrees}deg) translateX(calc(var(--size) / 2 - var(--thickness) - 0.2rem))`;
    },
    /** Shows the spin b's arrow while the hole spins in the b-plane. */
    setSpin(spinning) {
      spin.style.display = spinning ? '' : 'none';
    },
  };
}
