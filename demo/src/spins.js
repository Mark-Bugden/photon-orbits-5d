// The black hole's spins, a and b, side by side, held to a + b ≤ 1: beyond that there is no horizon. At a + b = 1
// the hole is extremal, and a label says so beside the sliders.
//
// A control with the kit's shape, { el, value, disabled, onChange }, where value is [a, b]. Listeners get
// (value, { held }), where held says the slider just moved was stopped at a + b = 1.

import { slider } from '@mark-bugden/demo-kit';

import { extremal } from './momenta.js';

const format = (value) => value.toFixed(2);

export function spins({ value: [a0, b0] }) {
  const sliders = {
    a: slider({ label: 'a', min: 0, max: 0.99, step: 0.01, value: a0, format }),
    b: slider({ label: 'b', min: 0, max: 0.99, step: 0.01, value: b0, format }),
  };
  sliders.a.el.querySelector('input').setAttribute('aria-label', 'Spin a');
  sliders.b.el.querySelector('input').setAttribute('aria-label', 'Spin b');
  const listeners = new Set();

  const label = document.createElement('span');
  label.className = 'dk-control-label';
  label.textContent = 'Spins';
  const badge = document.createElement('span');
  badge.className = 'po-badge';
  badge.textContent = 'Extremal';
  const head = document.createElement('div');
  head.className = 'dk-control-head';
  head.append(label, badge);
  const pair = document.createElement('div');
  pair.className = 'po-pair';
  pair.append(sliders.a.el, sliders.b.el);
  const el = document.createElement('div');
  el.className = 'dk-control po-group';
  el.setAttribute('role', 'group');
  el.setAttribute('aria-label', 'Spins');
  el.append(head, pair);

  const read = () => [sliders.a.value, sliders.b.value];
  const showBadge = () => (badge.hidden = !extremal(...read()));
  showBadge();

  for (const [moved, other] of [['a', 'b'], ['b', 'a']]) {
    sliders[moved].onChange(() => {
      const held = sliders.a.value + sliders.b.value > 1 + 1e-9;
      if (held) sliders[moved].value = Math.round(100 * (1 - sliders[other].value)) / 100;
      showBadge();
      for (const listener of listeners) listener(read(), { held });
    });
  }

  return {
    el,
    get value() {
      return read();
    },
    set value([a, b]) {
      sliders.a.value = a;
      sliders.b.value = b;
      showBadge();
    },
    get disabled() {
      return sliders.a.disabled;
    },
    set disabled(state) {
      sliders.a.disabled = state;
      sliders.b.disabled = state;
    },
    onChange(listener) {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
  };
}
