// The colour legend: the paper's wheel as a ring in a corner of the 3D view, labelled 0, π/2, π and 3π/2, with a
// needle at the photon's ψ. ψ grows counter-clockwise from 0 on the right, as in the paper's Fig 4.1.
//
// The wheel is Mathematica's Hue at full saturation and value, which is linear in sRGB between red, yellow, green,
// cyan, blue and magenta, so a conic gradient through those six draws it exactly.

const LABELS = [['0', 0], ['π/2', 90], ['π', 180], ['3π/2', 270]]; // and their angles, counter-clockwise

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
  };
}
