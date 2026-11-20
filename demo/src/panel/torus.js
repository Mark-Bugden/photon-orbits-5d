// The (φ, ψ) torus: φ across and ψ up, each from 0 to 2π, with the trail wrapping round at the edges. It shows ψ
// by position rather than by colour, and how the path winds: a closed orbit is a closed curve, and any other
// slowly covers the square. The trail fades towards its tail, as in the 3D view, and the wheel along the ψ axis
// gives each height its colour there.
//
// A point closer than about half a pixel to the last one drawn is skipped, so a long trail stays quick to draw.
// Nothing is allocated per frame.

import { hue } from '../colour.js';
import { trace } from '../trail.js';

const PI = Math.PI;
const TWO_PI = 2 * PI;
const LEFT = 50; // CSS pixels: ψ's name, ticks and the wheel
const RIGHT = 12;
const TOP = 8;
const BOTTOM = 34; // φ's ticks and name
const BAR = 6; // the wheel's width
const FADED = 0.15; // the opacity left at the tail, as in the 3D view
const FULL = 0.5; // the newest half of the trail is at full opacity
const STEPS = 8; // the fade, in steps of opacity over the older half
const NEAR = 0.6; // CSS pixels
const MOST_POINTS = 20000; // of the trail: Fig 4.10's longest has 18,000; more, near the extremal edge, are thinned
const FONT = '500 11px Inter, system-ui, sans-serif';
const NAME_FONT = '500 14px Inter, system-ui, sans-serif';
const TICKS = ['0', 'π', '2π'];
const WHEEL_STOPS = [0, 1, 2, 3, 4, 5, 6].map((k) => `rgb(${hue(k / 6).map((c) => Math.round(255 * c)).join()})`);

const opacity = (step) => (step >= STEPS ? 1 : FADED + ((1 - FADED) * (step + 0.5)) / STEPS);

export function createTorus() {
  let wheel = null; // the ψ axis's gradient, for the square at wheelTop and wheelSide
  let wheelTop = NaN;
  let wheelSide = NaN;

  // What visit() needs, set before each trace.
  let ctx = null;
  let x0 = 0; // the square's left and top, and its side, in CSS pixels
  let y0 = 0;
  let side = 1;
  let tail = 0;
  let length = 1;
  let phi0 = 0;
  let psi0 = 0;
  let started = false;
  let step = 0; // of the fade
  let cellU = 0; // which copy of the square the trail is in: ⌊(φ + φ₀)/2π⌋, ⌊(ψ + ψ₀)/2π⌋
  let cellV = 0;
  let lastU = 0; // the last point visited, unwrapped
  let lastV = 0;
  let drawnX = 0; // and the last drawn, in CSS pixels
  let drawnY = 0;

  const X = (u) => x0 + ((u - TWO_PI * cellU) / TWO_PI) * side;
  const Y = (v) => y0 + side - ((v - TWO_PI * cellV) / TWO_PI) * side;
  function moveTo(u, v) {
    drawnX = X(u);
    drawnY = Y(v);
    ctx.moveTo(drawnX, drawnY);
  }
  function lineTo(u, v) {
    drawnX = X(u);
    drawnY = Y(v);
    ctx.lineTo(drawnX, drawnY);
  }

  function visit(lam, theta, phi, psi) {
    const u = phi + phi0;
    const v = psi + psi0;
    const s = Math.max(0, Math.min(STEPS, Math.floor(((lam - tail) / length) * (STEPS / FULL))));
    if (!started) {
      started = true;
      step = s;
      cellU = Math.floor(u / TWO_PI);
      cellV = Math.floor(v / TWO_PI);
      ctx.globalAlpha = opacity(s);
      ctx.beginPath();
      moveTo(u, v);
      lastU = u;
      lastV = v;
      return;
    }
    // Where the step leaves the square, finish it at that edge and carry on from the opposite one.
    for (let wraps = 0; wraps < 2; wraps++) {
      const nextU = Math.floor(u / TWO_PI);
      const nextV = Math.floor(v / TWO_PI);
      if (nextU === cellU && nextV === cellV) break;
      const tu = nextU === cellU ? Infinity : (TWO_PI * (nextU > cellU ? cellU + 1 : cellU) - lastU) / (u - lastU);
      const tv = nextV === cellV ? Infinity : (TWO_PI * (nextV > cellV ? cellV + 1 : cellV) - lastV) / (v - lastV);
      const t = Math.min(tu, tv);
      const edgeU = lastU + t * (u - lastU);
      const edgeV = lastV + t * (v - lastV);
      lineTo(edgeU, edgeV);
      if (tu <= tv) cellU += nextU > cellU ? 1 : -1;
      if (tv <= tu) cellV += nextV > cellV ? 1 : -1;
      moveTo(edgeU, edgeV);
      lastU = edgeU;
      lastV = edgeV;
    }
    if (s !== step) {
      lineTo(u, v);
      ctx.stroke();
      step = s;
      ctx.globalAlpha = opacity(s);
      ctx.beginPath();
      moveTo(u, v);
    } else if (Math.abs(X(u) - drawnX) + Math.abs(Y(v) - drawnY) >= NEAR) {
      lineTo(u, v);
    }
    lastU = u;
    lastV = v;
  }

  /** Draws the trail from λ = head − length to head; the same arguments as the strip's draw. */
  function draw(context, width, height, ratio, path, head, lengthNow, lam0, phi0Now, psi0Now, colours) {
    side = Math.max(1, Math.floor(Math.min(width - LEFT - RIGHT, height - TOP - BOTTOM)));
    x0 = Math.round(LEFT + 0.5 * (width - LEFT - RIGHT - side));
    y0 = Math.round(TOP + 0.5 * (height - TOP - BOTTOM - side));
    const right = x0 + side;
    const bottom = y0 + side;
    const middleX = x0 + 0.5 * side;
    const middleY = y0 + 0.5 * side;

    context.fillStyle = colours.bg;
    context.fillRect(x0, y0, side, side);
    context.lineWidth = 1;
    context.strokeStyle = colours.border;
    context.beginPath();
    context.moveTo(Math.round(middleX) + 0.5, y0);
    context.lineTo(Math.round(middleX) + 0.5, bottom);
    context.moveTo(x0, Math.round(middleY) + 0.5);
    context.lineTo(right, Math.round(middleY) + 0.5);
    context.stroke();
    context.strokeRect(x0 + 0.5, y0 + 0.5, side - 1, side - 1);

    // The wheel along the ψ axis, so each height has its colour.
    if (wheelTop !== y0 || wheelSide !== side) {
      wheel = context.createLinearGradient(0, bottom, 0, y0);
      WHEEL_STOPS.forEach((colour, k) => wheel.addColorStop(k / 6, colour));
      wheelTop = y0;
      wheelSide = side;
    }
    context.fillStyle = wheel;
    context.fillRect(x0 - BAR - 4, y0, BAR, side);

    // Ticks and names.
    context.font = FONT;
    context.fillStyle = colours.muted;
    context.textAlign = 'center';
    context.textBaseline = 'top';
    for (let k = 0; k < 3; k++) context.fillText(TICKS[k], x0 + 0.5 * k * side, bottom + 4);
    context.textAlign = 'right';
    context.textBaseline = 'middle';
    for (let k = 0; k < 3; k++) context.fillText(TICKS[k], x0 - BAR - 8, bottom - 0.5 * k * side);
    context.font = NAME_FONT;
    context.fillStyle = colours.curves[1];
    context.textAlign = 'center';
    context.textBaseline = 'alphabetic';
    context.fillText('φ', middleX, bottom + BOTTOM - 3);
    context.fillStyle = colours.curves[2];
    context.textAlign = 'left';
    context.textBaseline = 'middle';
    context.fillText('ψ', 4, middleY);

    // The trail, wrapping at the edges, and the photon. The clip leaves room for a trail along an edge, at φ = 0
    // or ψ = 0, to show its whole width.
    ctx = context;
    tail = head - lengthNow;
    length = lengthNow;
    phi0 = phi0Now;
    psi0 = psi0Now;
    started = false;
    context.save();
    context.beginPath();
    context.rect(x0 - 1, y0 - 1, side + 2, side + 2);
    context.clip();
    context.lineWidth = 1.25;
    context.lineJoin = 'round';
    context.strokeStyle = colours.text;
    trace(path, head, lengthNow, visit, MOST_POINTS);
    lineTo(lastU, lastV);
    context.stroke();
    context.globalAlpha = 1;
    context.restore();
    context.fillStyle = colours.text;
    context.strokeStyle = colours.bg;
    context.lineWidth = 2;
    context.beginPath();
    context.arc(X(lastU), Y(lastV), 4, 0, TWO_PI);
    context.stroke();
    context.fill();
    ctx = null;
  }

  return { draw };
}
