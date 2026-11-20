// Ticks for the panel's axes: λ in round numbers, and the angles in multiples of π. A tick set is reused frame
// after frame, and a tick's label is only written again when the tick moves, so drawing allocates nothing per
// frame while the axes stay put.

const PI = Math.PI;
const NICE = [1, 2, 5, 10];
const QUARTERS = [1, 2, 4, 8, 20, 40, 80, 200, 400, 800, 2000, 4000, 8000, 20000]; // steps, in multiples of π/4
const MOST = 12;

/** An empty tick set: { count, step, values, labels }. */
export function tickSet() {
  return { count: 0, step: 0, angles: false, values: new Float64Array(MOST), labels: new Array(MOST).fill('') };
}

/** Round-number ticks from lo to hi, about `count` of them, e.g. 0, 2, 4, 6. */
export function linearTicks(lo, hi, count, ticks) {
  const rough = (hi - lo) / Math.max(1, count);
  const power = 10 ** Math.floor(Math.log10(rough));
  let k = 0;
  while (k < NICE.length - 1 && NICE[k] * power < rough) k++;
  return place(lo, hi, NICE[k] * power, false, ticks);
}

/** Ticks at multiples of π/4, π/2, π, 2π, 5π, …, from lo to hi: as close together as fit `count` of them. */
export function angleTicks(lo, hi, count, ticks) {
  let k = 0;
  while (k < QUARTERS.length - 1 && Math.floor((4 * hi) / (QUARTERS[k] * PI)) - Math.ceil((4 * lo) / (QUARTERS[k] * PI)) + 1 > count) k++;
  return place(lo, hi, (QUARTERS[k] * PI) / 4, true, ticks);
}

function place(lo, hi, step, angles, ticks) {
  const relabel = step !== ticks.step || angles !== ticks.angles;
  ticks.step = step;
  ticks.angles = angles;
  // A tick within rounding of an end counts as in range, e.g. λ = 5 at the end of the trail Fig 4.11 shows.
  let n = 0;
  for (let k = Math.ceil(lo / step - 1e-4); k * step <= hi + 1e-4 * step && n < MOST; k++) {
    const value = k * step + 0; // + 0 turns −0 into 0
    if (relabel || n >= ticks.count || ticks.values[n] !== value) {
      ticks.values[n] = value;
      ticks.labels[n] = angles ? piLabel(value) : numberLabel(value, step);
    }
    n++;
  }
  ticks.count = n;
  return ticks;
}

/** A multiple of π/4 as text: 0, π/4, π/2, 3π/4, π, −2π, 40π. */
export function piLabel(value) {
  const quarters = Math.round((4 * value) / PI);
  if (quarters === 0) return '0';
  const sign = quarters < 0 ? '−' : '';
  const q = Math.abs(quarters);
  const [top, bottom] = q % 4 === 0 ? [q / 4, 1] : q % 2 === 0 ? [q / 2, 2] : [q, 4];
  return `${sign}${top === 1 ? '' : top}π${bottom === 1 ? '' : `/${bottom}`}`;
}

/** A number as text, with as many decimals as its step needs: 0, 2.5, −1. */
export function numberLabel(value, step) {
  const decimals = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return (Math.abs(value) < step / 1e6 ? 0 : value).toFixed(decimals).replace('-', '−');
}
