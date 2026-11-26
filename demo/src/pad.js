// The momentum pad: Φ across and Ψ up, on fixed axes from −3 to 3, with the momenta that have a spherical orbit
// shaded and the current orbit as a dot to drag. The region reshapes as the spins change, and the dot keeps its
// direction from the origin and its fraction of the way to the edge, so the orbit keeps its character. Dragging
// snaps to the axes, where the special cases Φ = 0 and Ψ = 0 are. Sliders for Φ and Ψ below it give fine control,
// and the arrow keys move the dot; both stop at the edge.
//
// A control with the kit's shape, { el, value, disabled, onChange }, where value is [Φ, Ψ]; plus setSpins(a, b).

import { onThemeChange, ramp, slider, token } from '@mark-bugden/demo-kit';

import { clampAlong, locate, place, toRegion } from './momenta.js';
import { boundary, solveOrbit } from './physics/index.js';

const RANGE = 3; // the axes run from −RANGE to RANGE
const SNAP = 6; // CSS pixels from an axis within which a drag lands exactly on it
const STEP = 0.01; // an arrow key's move; with Shift, ten times as far
const MARGIN = { left: 20, right: 6, top: 6, bottom: 20 }; // around the plot, for the axis labels
const FONT = '500 11px Inter, system-ui, sans-serif';
const GRID = 56; // the region's tint is the orbit's radius on a GRID × GRID grid over the pad, smoothed
const RADII = [0.5, 2]; // the radii at the ends of the tint's scale: r runs from about 0.1 to 2 over all spins
const TINT = 0.6; // the tint's opacity over the pad

const minus = (text) => text.replace('-', '−');
const format = (value) => minus(value.toFixed(2));

/** The pad, for the hole with spins a and b, with the dot at value = [Φ, Ψ]. */
export function pad({ spins: [a0, b0], value: [Phi0, Psi0] }) {
  let spins = [a0, b0];
  let value = [Phi0, Psi0];
  let where = locate(...spins, ...value);
  let region = boundary(...spins);
  let radii = radiiFor(...spins);
  let tint = null; // { id, canvas, strip }: the radii painted with the theme's spectrum, and the key's gradient
  const listeners = new Set();

  const canvas = document.createElement('canvas');
  canvas.className = 'po-pad-canvas';
  canvas.tabIndex = 0;
  canvas.setAttribute('role', 'application');
  canvas.setAttribute('aria-label', 'Momenta: Φ across, Ψ up. The shaded region has orbits. Arrow keys move the dot.');
  const sliders = {
    Phi: slider({ label: 'Φ', min: -RANGE, max: RANGE, step: STEP, value: Phi0, format }),
    Psi: slider({ label: 'Ψ', min: -RANGE, max: RANGE, step: STEP, value: Psi0, format }),
  };
  sliders.Phi.el.querySelector('input').setAttribute('aria-label', 'Momentum Φ');
  sliders.Psi.el.querySelector('input').setAttribute('aria-label', 'Momentum Ψ');

  const heading = document.createElement('span');
  heading.className = 'dk-control-label';
  heading.textContent = 'Momenta';
  const pair = document.createElement('div');
  pair.className = 'po-pair';
  pair.append(sliders.Phi.el, sliders.Psi.el);
  const el = document.createElement('div');
  el.className = 'dk-control po-group po-pad';
  el.setAttribute('role', 'group');
  el.setAttribute('aria-label', 'Momenta');
  el.append(heading, canvas, pair);

  // Sets the dot without telling listeners, and shows it.
  function show(next, remember = true) {
    value = next;
    if (remember) where = locate(...spins, ...value);
    sliders.Phi.value = value[0];
    sliders.Psi.value = value[1];
    draw();
  }
  function change(next) {
    show(next);
    for (const listener of listeners) listener(value, { stopped: false });
  }

  // The sliders stop at the edge, along their own axis.
  for (const [name, axis] of [['Phi', 0], ['Psi', 1]]) {
    sliders[name].onChange((input) => {
      const target = axis === 0 ? [input, value[1]] : [value[0], input];
      const next = clampAlong(...spins, value, target);
      show(next);
      const stopped = next[axis] !== input;
      for (const listener of listeners) listener(value, { stopped, axis: name });
    });
  }

  // Dragging. Points outside the region go to its edge along their ray from the origin.
  let size = 0;
  const plot = () => size - MARGIN.left - MARGIN.right;
  const toMomenta = (x, y) => [
    ((x - MARGIN.left) / plot()) * 2 * RANGE - RANGE,
    RANGE - ((y - MARGIN.top) / plot()) * 2 * RANGE,
  ];
  const perPixel = () => (2 * RANGE) / plot();
  function drag(event) {
    const box = canvas.getBoundingClientRect();
    let [Phi, Psi] = toMomenta(event.clientX - box.left, event.clientY - box.top);
    if (Math.abs(Phi) < SNAP * perPixel()) Phi = 0;
    if (Math.abs(Psi) < SNAP * perPixel()) Psi = 0;
    if (Phi === 0 && Psi === 0) change([0, 0]);
    else change(toRegion(...spins, Phi, Psi));
  }
  // The browser focuses the canvas on a press, as it does any element: a script's focus() would count as
  // :focus-visible in Chromium, and show the focus ring after a click.
  canvas.addEventListener('pointerdown', (event) => {
    if (el.hasAttribute('data-disabled') || event.button !== 0) return;
    canvas.setPointerCapture(event.pointerId);
    drag(event);
  });
  canvas.addEventListener('pointermove', (event) => {
    if (canvas.hasPointerCapture(event.pointerId)) drag(event);
  });

  const KEYS = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowDown: [0, -1], ArrowUp: [0, 1] };
  canvas.addEventListener('keydown', (event) => {
    const move = KEYS[event.key];
    if (!move || event.altKey || event.ctrlKey || event.metaKey || el.hasAttribute('data-disabled')) return;
    event.preventDefault();
    const step = event.shiftKey ? 10 * STEP : STEP;
    // Round to the step, so the keys land on the sliders' values, and on 0.
    const target = value.map((v, i) => Math.round((v + move[i] * step) / STEP) * STEP);
    change(clampAlong(...spins, value, [move[0] ? target[0] : value[0], move[1] ? target[1] : value[1]]));
  });
  // The focus ring shows while the canvas matches :focus-visible, as the kit's outlines do: after Tab, not after a
  // press. A key can turn :focus-visible on, so each key redraws if the ring should change.
  let ring = false;
  canvas.addEventListener('keydown', () => {
    if (canvas.matches(':focus-visible') !== ring) draw();
  });
  canvas.addEventListener('focus', () => draw());
  canvas.addEventListener('blur', () => draw());

  function draw() {
    const box = canvas.getBoundingClientRect();
    if (box.width === 0) return;
    size = box.width;
    const ratio = Math.min(globalThis.devicePixelRatio || 1, 2);
    const pixels = Math.round(size * ratio);
    if (canvas.width !== pixels || canvas.height !== pixels) {
      canvas.width = pixels;
      canvas.height = pixels;
    }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.clearRect(0, 0, size, size);
    const colour = (name) => token(el, name);
    const side = plot();
    const X = (Phi) => MARGIN.left + ((Phi + RANGE) / (2 * RANGE)) * side;
    const Y = (Psi) => MARGIN.top + ((RANGE - Psi) / (2 * RANGE)) * side;

    ctx.fillStyle = colour('--color-bg');
    ctx.fillRect(MARGIN.left, MARGIN.top, side, side);
    ctx.lineWidth = 1;
    ctx.strokeStyle = colour('--color-border');
    ctx.beginPath();
    for (let k = -RANGE + 1; k < RANGE; k++) {
      if (k === 0) continue;
      ctx.moveTo(Math.round(X(k)) + 0.5, MARGIN.top);
      ctx.lineTo(Math.round(X(k)) + 0.5, MARGIN.top + side);
      ctx.moveTo(MARGIN.left, Math.round(Y(k)) + 0.5);
      ctx.lineTo(MARGIN.left + side, Math.round(Y(k)) + 0.5);
    }
    ctx.stroke();
    ctx.strokeRect(MARGIN.left + 0.5, MARGIN.top + 0.5, side - 1, side - 1);

    // The region, tinted by the orbit's radius, lighter for larger.
    const spectrum = [1, 2, 3, 4].map((k) => colour(`--chart-spectrum-${k}`));
    if (tint?.id !== spectrum.join() + radii.id) tint = paint(radii, spectrum);
    ctx.save();
    ctx.beginPath();
    region.forEach(([Phi, Psi], k) => (k ? ctx.lineTo(X(Phi), Y(Psi)) : ctx.moveTo(X(Phi), Y(Psi))));
    ctx.closePath();
    ctx.clip();
    ctx.globalAlpha = TINT;
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(tint.canvas, MARGIN.left, MARGIN.top, side, side);
    ctx.restore();
    ctx.globalAlpha = 0.7;
    ctx.lineWidth = 1.25;
    ctx.strokeStyle = colour('--color-text');
    ctx.stroke();
    ctx.globalAlpha = 1;

    // The tint's key, in the top-left corner, which the region never reaches.
    const key = { x: MARGIN.left + 8, y: MARGIN.top + 8, width: Math.min(56, side / 4), height: 5 };
    ctx.globalAlpha = TINT;
    ctx.drawImage(tint.strip, key.x, key.y + 14, key.width, key.height);
    ctx.globalAlpha = 1;

    // The axes, brighter when the dot is on one: the special cases.
    for (const [onAxis, x0, y0, x1, y1] of [
      [value[1] === 0, X(-RANGE), Y(0), X(RANGE), Y(0)],
      [value[0] === 0, X(0), Y(-RANGE), X(0), Y(RANGE)],
    ]) {
      ctx.lineWidth = onAxis ? 1.5 : 1;
      ctx.strokeStyle = colour(onAxis ? '--color-text' : '--color-curve');
      ctx.beginPath();
      ctx.moveTo(x0, y0);
      ctx.lineTo(x1, y1);
      ctx.stroke();
    }

    // The labels: Φ and Ψ, and the halves with and against each spin. Co-rotating momenta are negative.
    ctx.font = FONT;
    ctx.fillStyle = colour('--color-muted');
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'left';
    ctx.fillText(`r ${RADII[0]}`, key.x, key.y + 5);
    ctx.textAlign = 'right';
    ctx.fillText(`${RADII[1]}`, key.x + key.width, key.y + 5);
    const below = MARGIN.top + side + MARGIN.bottom / 2;
    ctx.textAlign = 'center';
    ctx.fillText('Φ', X(0), below);
    if (spins[0] > 0) {
      ctx.textAlign = 'left';
      ctx.fillText('with a', MARGIN.left, below);
      ctx.textAlign = 'right';
      ctx.fillText('against a', MARGIN.left + side, below);
    }
    ctx.textAlign = 'center';
    ctx.fillText('Ψ', MARGIN.left / 2, Y(0));
    ctx.save();
    ctx.translate(MARGIN.left / 2, 0);
    ctx.rotate(-Math.PI / 2);
    if (spins[1] > 0) {
      ctx.textAlign = 'left';
      ctx.fillText('with b', -(MARGIN.top + side), 0);
      ctx.textAlign = 'right';
      ctx.fillText('against b', -MARGIN.top, 0);
    }
    ctx.restore();

    // The dot.
    const [x, y] = [X(value[0]), Y(value[1])];
    ring = canvas.matches(':focus-visible');
    if (ring) {
      ctx.lineWidth = 2;
      ctx.strokeStyle = colour('--color-focus');
      ctx.beginPath();
      ctx.arc(x, y, 11, 0, 2 * Math.PI);
      ctx.stroke();
    }
    ctx.fillStyle = colour('--color-bg');
    ctx.beginPath();
    ctx.arc(x, y, 7.5, 0, 2 * Math.PI);
    ctx.fill();
    ctx.fillStyle = colour('--color-text');
    ctx.beginPath();
    ctx.arc(x, y, 5.5, 0, 2 * Math.PI);
    ctx.fill();
  }

  new ResizeObserver(() => draw()).observe(canvas);
  onThemeChange(() => draw());

  return {
    el,
    get value() {
      return value;
    },
    /** Moves the dot without telling listeners, e.g. for a preset. */
    set value(next) {
      show(next);
    },
    get disabled() {
      return el.hasAttribute('data-disabled');
    },
    set disabled(state) {
      el.toggleAttribute('data-disabled', Boolean(state));
      sliders.Phi.disabled = state;
      sliders.Psi.disabled = state;
      canvas.tabIndex = state ? -1 : 0;
    },
    onChange(listener) {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
    /** Reshapes the region for new spins, and moves the dot to keep its direction and fraction of the way. */
    setSpins(a, b) {
      if (a === spins[0] && b === spins[1]) return;
      spins = [a, b];
      region = boundary(a, b);
      radii = radiiFor(a, b);
      show(place(a, b, where), false);
    },
  };
}

// The orbit's radius at the centre of each cell of a GRID × GRID grid over the pad, row by row from the top, or NaN
// where there is none. Just outside the region the radius still exists (the orbit lacks only a θ-band), so the
// smoothed tint runs on smoothly to the region's edge.
function radiiFor(a, b) {
  const values = new Float32Array(GRID * GRID);
  for (let row = 0; row < GRID; row++) {
    const Psi = RANGE - ((row + 0.5) / GRID) * 2 * RANGE;
    for (let column = 0; column < GRID; column++) {
      const Phi = ((column + 0.5) / GRID) * 2 * RANGE - RANGE;
      const { x } = solveOrbit(a, b, Phi, Psi);
      values[row * GRID + column] = x === null ? NaN : Math.sqrt(x);
    }
  }
  return Object.assign(values, { id: `${a},${b}` });
}

// The radii as an image, coloured through the spectrum's four stops, and the key's gradient as a strip.
function paint(radii, spectrum) {
  const colour = ramp(spectrum);
  const rgb = (r) => {
    const hex = colour(Math.min(1, Math.max(0, (r - RADII[0]) / (RADII[1] - RADII[0]))));
    return [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16));
  };
  const canvas = document.createElement('canvas');
  canvas.width = GRID;
  canvas.height = GRID;
  const ctx = canvas.getContext('2d');
  const image = ctx.createImageData(GRID, GRID);
  radii.forEach((r, k) => {
    if (Number.isNaN(r)) return;
    image.data.set([...rgb(r), 255], 4 * k);
  });
  ctx.putImageData(image, 0, 0);
  const strip = document.createElement('canvas');
  strip.width = 64;
  strip.height = 1;
  const row = strip.getContext('2d').createImageData(64, 1);
  for (let k = 0; k < 64; k++) row.data.set([...rgb(RADII[0] + ((k + 0.5) / 64) * (RADII[1] - RADII[0])), 255], 4 * k);
  strip.getContext('2d').putImageData(row, 0, 0);
  return { id: spectrum.join() + radii.id, canvas, strip };
}
