// The paper's colour wheel for ψ (Fig 4.1): hue ψ/2π at full saturation and value, as Mathematica's Hue. Red at 0,
// yellow at π/3, green at 2π/3, cyan at π, blue at 4π/3 and magenta at 5π/3.

const TWO_PI = 2 * Math.PI;
export const WHEEL_SIZE = 256;

/** [r, g, b] in sRGB, from 0 to 1, for the hue h (in turns) at full saturation and value. */
export function hue(h) {
  const channel = (n) => {
    const k = (((n + 6 * h) % 6) + 6) % 6;
    return 1 - Math.max(0, Math.min(k, 4 - k, 1));
  };
  return [channel(5), channel(3), channel(1)];
}

/** An sRGB channel, from 0 to 1, in linear light, as three.js wants vertex colours. */
export const toLinear = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);

/** The wheel as WHEEL_SIZE linear-RGB triples, entry k at ψ = 2πk/WHEEL_SIZE. */
export const WHEEL = new Float32Array(WHEEL_SIZE * 3);
for (let k = 0; k < WHEEL_SIZE; k++) WHEEL.set(hue(k / WHEEL_SIZE).map(toLinear), 3 * k);

/** The entry of WHEEL for the angle ψ, for any real ψ: the nearest entry, so 0 and 2π give the same one. */
export function wheelIndex(psi) {
  const turns = psi / TWO_PI;
  return Math.round((turns - Math.floor(turns)) * WHEEL_SIZE) % WHEEL_SIZE;
}
