// θ, φ and ψ against λ over the trail, as in the paper's Figs 4.9 and 4.11: three rows sharing the λ axis, from
// the trail's tail on the left to the photon on the right. θ keeps a fixed scale, with dashed lines at the band's
// edges; φ and ψ are unwrapped, so each of their rows scales to its stretch. All three are in radians. A band
// under ψ shows the trail's colour at each λ.
//
// Each curve is drawn through at most four values a pixel column (the first, lowest, highest and last there), so
// a trail of a hundred oscillations costs about the same to draw as one of a few. Nothing is allocated per frame.

import { hue } from '../colour.js';
import { trace } from '../trail.js';
import { angleTicks, linearTicks, tickSet } from './ticks.js';

const PI = Math.PI;
const MAX_COLUMNS = 4096; // device pixels across the plot
const MOST_POINTS = 20000; // of the trail: Fig 4.10's longest has 18,000; more, near the extremal edge, are thinned
const LEFT = 46; // CSS pixels: each row's name and its ticks
const RIGHT = 30; // the names of the band's edges
const TOP = 6;
const GAP = 12; // between the rows
const BAND = 6; // the colour band under ψ
const BAND_GAP = 3;
const BOTTOM = 18; // λ's ticks
const FONT = '500 11px Inter, system-ui, sans-serif';
const NAME_FONT = '500 14px Inter, system-ui, sans-serif';
const NAMES = ['θ', 'φ', 'ψ'];

// The paper's wheel in sRGB bytes, for the colour band.
const WHEEL_BYTES = new Uint8ClampedArray(256 * 3);
for (let k = 0; k < 256; k++) WHEEL_BYTES.set(hue(k / 256).map((c) => Math.round(255 * c)), 3 * k);

export function createStrip() {
  // Per curve, four values a column: first, lowest, highest, last; NaN where no point falls.
  const curves = [0, 1, 2].map(() => new Float64Array(4 * MAX_COLUMNS));
  const hues = new Float64Array(MAX_COLUMNS);
  const range = new Float64Array(6); // the lowest and highest θ, φ and ψ
  const lamTicks = tickSet();
  const rowTicks = [tickSet(), tickSet(), tickSet()];
  const edgeValues = new Float64Array(2);
  const edgeNames = ['', ''];
  const band = new ImageData(MAX_COLUMNS, 1);
  const bandCanvas = document.createElement('canvas');
  bandCanvas.width = MAX_COLUMNS;
  bandCanvas.height = 1;
  const bandContext = bandCanvas.getContext('2d');

  // What collect() needs, set before each trace.
  let tail = 0;
  let scale = 1; // columns per unit of λ
  let columns = 1;
  let phi0 = 0;
  let psi0 = 0;
  let pole = false;
  let equator = false;

  function put(row, c, value) {
    const curve = curves[row];
    const i = 4 * c;
    if (Number.isNaN(curve[i])) {
      curve[i] = value;
      curve[i + 1] = value;
      curve[i + 2] = value;
    } else if (value < curve[i + 1]) {
      curve[i + 1] = value;
    } else if (value > curve[i + 2]) {
      curve[i + 2] = value;
    }
    curve[i + 3] = value;
    if (value < range[2 * row]) range[2 * row] = value;
    if (value > range[2 * row + 1]) range[2 * row + 1] = value;
  }

  function collect(lam, theta, phi, psi) {
    let c = Math.floor((lam - tail) * scale);
    if (c < 0) c = 0;
    else if (c >= columns) c = columns - 1;
    put(0, c, theta);
    put(1, c, phi + phi0);
    put(2, c, psi + psi0);
    hues[c] = psi + psi0;
  }

  // The range a row shows, into range[2·row] and range[2·row + 1]. θ's is fixed by how the orbit meets the axes:
  // 0 to π/2, or on through the pole to −π/2, or through the equator to π. Through both, θ = χ/2 grows, and its
  // row fits its curve, as φ's and ψ's do.
  function limits(row) {
    if (row === 0 && !(pole && equator)) {
      range[0] = pole ? -0.5 * PI : 0;
      range[1] = equator ? PI : 0.5 * PI;
      return;
    }
    const span = Math.max(range[2 * row + 1] - range[2 * row], 0.5 * PI);
    const middle = 0.5 * (range[2 * row] + range[2 * row + 1]);
    range[2 * row] = middle - 0.54 * span;
    range[2 * row + 1] = middle + 0.54 * span;
  }

  // The band's edges, which the θ row marks with dashed lines: θ_min and θ_max, or, where θ runs on through an
  // axis, the extremes it swings between. Returns how many.
  function edges(path) {
    const low = path.band[0];
    const high = path.band[1];
    if (pole && equator) return 0;
    edgeValues[0] = low;
    edgeNames[0] = 'θ₋';
    edgeValues[1] = high;
    edgeNames[1] = 'θ₊';
    if (pole) {
      edgeValues[0] = -high;
      edgeNames[0] = '−θ₊';
    } else if (equator) {
      edgeValues[1] = PI - low;
      edgeNames[1] = 'π−θ₋';
    }
    return 2;
  }

  // Strokes a row's curve through its columns: value lo at `bottom`, rising k CSS pixels per unit.
  function stroke(ctx, curve, ratio, bottom, lo, k) {
    ctx.beginPath();
    let started = false;
    for (let c = 0; c < columns; c++) {
      const i = 4 * c;
      if (Number.isNaN(curve[i])) continue;
      const x = LEFT + (c + 0.5) / ratio;
      if (started) ctx.lineTo(x, bottom - k * (curve[i] - lo));
      else ctx.moveTo(x, bottom - k * (curve[i] - lo));
      started = true;
      if (curve[i + 2] > curve[i + 1]) {
        ctx.lineTo(x, bottom - k * (curve[i + 1] - lo));
        ctx.lineTo(x, bottom - k * (curve[i + 2] - lo));
      }
      ctx.lineTo(x, bottom - k * (curve[i + 3] - lo));
    }
    ctx.stroke();
  }

  /**
   * Draws the trail from λ = head − length to head on ctx, which is width × height CSS pixels at `ratio` device
   * pixels each. The axes read λ − lam0, φ + phi0 and ψ + psi0: the presets set these to match the paper's plots.
   */
  function draw(ctx, width, height, ratio, path, head, length, lam0, phi0Now, psi0Now, colours) {
    const plot = width - LEFT - RIGHT;
    const right = LEFT + plot;
    columns = Math.max(1, Math.min(MAX_COLUMNS, Math.round(plot * ratio)));
    tail = head - length;
    scale = columns / length;
    phi0 = phi0Now;
    psi0 = psi0Now;
    pole = path.pole;
    equator = path.equator;
    for (let row = 0; row < 3; row++) {
      curves[row].fill(NaN, 0, 4 * columns);
      range[2 * row] = Infinity;
      range[2 * row + 1] = -Infinity;
    }
    hues.fill(NaN, 0, columns);
    trace(path, head, length, collect, MOST_POINTS);

    const rowHeight = (height - TOP - BOTTOM - BAND - BAND_GAP - 2 * GAP) / 3;
    const perLam = plot / length;
    const lamLeft = tail - lam0; // the λ the axis reads at its left end
    linearTicks(lamLeft, head - lam0, Math.max(2, Math.floor(plot / 50)), lamTicks);
    ctx.font = FONT;
    ctx.lineWidth = 1;

    for (let row = 0; row < 3; row++) {
      const top = Math.round(TOP + row * (rowHeight + GAP));
      const bottom = Math.round(top + rowHeight);
      limits(row);
      const lo = range[2 * row];
      const k = (bottom - top) / (range[2 * row + 1] - lo);
      // θ's ticks are fractions of π: the pole at 0 and the equator at π/2. φ's and ψ's, which can run to
      // hundreds, are plain numbers, as on the paper's axes.
      const ticks = rowTicks[row];
      const most = Math.max(3, Math.min(5, Math.floor(rowHeight / 28)));
      if (row === 0) angleTicks(lo, range[1], most, ticks);
      else linearTicks(lo, range[2 * row + 1], most - 1, ticks);

      // The row, its grid, and its ticks.
      ctx.fillStyle = colours.bg;
      ctx.fillRect(LEFT, top, plot, bottom - top);
      ctx.strokeStyle = colours.border;
      ctx.beginPath();
      for (let t = 0; t < ticks.count; t++) {
        const at = Math.round(bottom - k * (ticks.values[t] - lo)) + 0.5;
        ctx.moveTo(LEFT, at);
        ctx.lineTo(right, at);
      }
      for (let t = 0; t < lamTicks.count; t++) {
        const at = Math.round(LEFT + (lamTicks.values[t] - lamLeft) * perLam) + 0.5;
        ctx.moveTo(at, top);
        ctx.lineTo(at, bottom);
      }
      ctx.stroke();
      ctx.strokeRect(LEFT + 0.5, top + 0.5, plot - 1, bottom - top - 1);
      ctx.fillStyle = colours.muted;
      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      for (let t = 0; t < ticks.count; t++) ctx.fillText(ticks.labels[t], LEFT - 5, bottom - k * (ticks.values[t] - lo));
      ctx.font = NAME_FONT;
      ctx.fillStyle = colours.curves[row];
      ctx.textAlign = 'left';
      ctx.fillText(NAMES[row], 4, 0.5 * (top + bottom));
      ctx.font = FONT;

      // θ's band, dashed and named in the right margin.
      const count = row === 0 ? edges(path) : 0;
      if (count) {
        ctx.setLineDash([4, 3]);
        ctx.strokeStyle = colours.muted;
        ctx.beginPath();
        for (let e = 0; e < count; e++) {
          const at = Math.round(bottom - k * (edgeValues[e] - lo)) + 0.5;
          ctx.moveTo(LEFT, at);
          ctx.lineTo(right, at);
        }
        ctx.stroke();
        ctx.setLineDash([]);
        // A circular orbit's edges coincide: name the line once.
        ctx.fillStyle = colours.muted;
        ctx.textAlign = 'left';
        let named = Infinity;
        for (let e = 0; e < count; e++) {
          const at = bottom - k * (edgeValues[e] - lo);
          if (Math.abs(at - named) < 9) continue;
          ctx.fillText(edgeNames[e], right + 6, at);
          named = at;
        }
      }

      // The curve, and the photon at its right-hand end.
      ctx.save();
      ctx.beginPath();
      ctx.rect(LEFT, top - 5, plot + 5, bottom - top + 10);
      ctx.clip();
      ctx.lineWidth = 1.75;
      ctx.lineJoin = 'round';
      ctx.strokeStyle = colours.curves[row];
      stroke(ctx, curves[row], ratio, bottom, lo, k);
      ctx.fillStyle = colours.curves[row];
      ctx.strokeStyle = colours.bg;
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(LEFT + (columns - 0.5) / ratio, bottom - k * (curves[row][4 * columns - 1] - lo), 3.5, 0, 2 * PI);
      ctx.stroke();
      ctx.fill();
      ctx.restore();
      ctx.lineWidth = 1;
    }

    // The colour band: the trail's colour at each λ.
    let last = NaN;
    for (let c = 0; c < columns; c++) {
      if (!Number.isNaN(hues[c])) last = hues[c];
      const turns = last / (2 * PI);
      const k = 3 * (Math.round((turns - Math.floor(turns)) * 256) % 256);
      band.data[4 * c] = WHEEL_BYTES[k];
      band.data[4 * c + 1] = WHEEL_BYTES[k + 1];
      band.data[4 * c + 2] = WHEEL_BYTES[k + 2];
      band.data[4 * c + 3] = 255;
    }
    bandContext.putImageData(band, 0, 0, 0, 0, columns, 1);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(bandCanvas, 0, 0, columns, 1, LEFT, Math.round(TOP + 3 * rowHeight + 2 * GAP + BAND_GAP), plot, BAND);

    // λ's ticks, under the band.
    ctx.fillStyle = colours.muted;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'alphabetic';
    for (let t = 0; t < lamTicks.count; t++) {
      const x = LEFT + (lamTicks.values[t] - lamLeft) * perLam;
      if (x > LEFT - 0.5 && x < right + 0.5) ctx.fillText(lamTicks.labels[t], x, height - 4);
    }
    ctx.textAlign = 'right';
    ctx.fillText('λ', width - 2, height - 4);
  }

  return { draw };
}
