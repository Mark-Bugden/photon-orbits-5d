// The presets: the paper's figures, each matched to the figure itself, and two special holes.
//
// Each sets the spins and momenta, and:
// - span: the figure's stretch of path, as [from, to] in the angle variable χ, which is 0 at θ_min and grows by 2π
//   each oscillation. The photon starts at `to`. The trail's length is the visitor's, so `from` only records the
//   figure's, except for the opening state, whose trail reaches back to it.
// - psi0: ψ at χ = 0, which sets the colours and the side panel's ψ.
// - view: the direction the camera looks from, as [azimuth, elevation] in degrees about the hole, with z up. The
//   zoom is the visitor's.
// - lam0, phi0 (optional, 0 if left out): λ at the paper's λ = 0, and φ at χ = 0, so that the side panel's axes
//   read as the paper's plots against λ. Only Figs 4.8 and 4.10 have them, from Figs 4.9 and 4.11. The 3D view
//   needs no φ₀: it is in the camera's azimuth.
//
// The figures' spans, colours and viewpoints were fitted to the paper's images: the outline and equator ring give
// the elevation, and the path's pixels and their colours the rest (docs/research/scripts/presets/). The fitted
// azimuth includes the start's φ and the constant twist arctan(a/r) of the paper's projection, which the demo
// leaves out: turning the path about z is the same as turning the camera. For Fig 4.10, the start and colours come
// from Fig 4.11 instead, which plots the same stretch against λ: λ from 1 to 5.

const PI = Math.PI;
const DEFAULT_VIEW = [-61.56, 36.2]; // Mathematica's default viewpoint, (1.3, −2.4, 2)

export const PRESETS = [
  { name: 'Fig 4.2: zero momenta', a: 0.5, b: 0, Phi: 0, Psi: 0, span: [PI, 5 * PI], psi0: 0, view: [287.6, 18.1] },
  { name: 'Fig 4.3: Φ = 2', a: 0.5, b: 0, Phi: 2, Psi: 0, span: [4.765, 318.893], psi0: 0, view: [231, 12.3] },
  { name: 'Fig 4.4: Φ = 0, Ψ = 0.5', a: 0.8, b: 0, Phi: 0, Psi: 0.5, span: [1.656, 14.67], psi0: 0.389, view: [19.2, 33.6] },
  { name: 'Fig 4.5: the same, longer', a: 0.8, b: 0, Phi: 0, Psi: 0.5, span: [3.09, 170.017], psi0: 1.688, view: [254.9, 33.5] },
  { name: 'Fig 4.6: the “football”', a: 0, b: 0.99, Phi: 0, Psi: 0, span: [0, 4 * PI], psi0: 2.725, view: [21.2, 0] },
  { name: 'Fig 4.7: two spins, a band', a: 0.6, b: 0.3, Phi: 1, Psi: -1, span: [3.523, 136.789], psi0: 3.434, view: [303.8, 22.8] },
  {
    name: 'Fig 4.8: cusps, extremal', a: 0.9, b: 0.1, Phi: 0.8, Psi: 0.8, span: [0, 33.696], psi0: 3.683, view: [104.1, 29.1],
    lam0: 0.076, phi0: 1.6122,
  },
  {
    name: 'Fig 4.10: φ and ψ turn back', a: 0.9, b: 0.09, Phi: 0.1, Psi: 0.01, span: [4.2418, 13.262], psi0: -0.0249, view: [117.1, 28.5],
    lam0: 0.7509, phi0: -1.5887,
  },
  // Fig 4.12 is on the edge of the region, a circular orbit at θ = 45°: one turn of φ.
  { name: 'Fig 4.12: constant θ', a: 0.5, b: 0.25, Phi: -0.929241, Psi: 0.903667, span: [0, 9.8301], psi0: 1.099, view: [232, 32] },
  { name: 'Equal spins (a = b)', a: 0.4, b: 0.4, Phi: 1.2, Psi: -0.5, span: [0, 12 * PI], psi0: 0, view: DEFAULT_VIEW },
  { name: 'No spin', a: 0, b: 0, Phi: 0.8, Psi: 0.6, span: [0, 4 * PI], psi0: 0, view: DEFAULT_VIEW },
];

// The opening state (Mark, 2026-11-20): Fig 4.5's orbit, moment, colours and viewpoint, with only the last three
// oscillations of trail rather than the figure's 27.
const PETALS = PRESETS.find((preset) => preset.name.startsWith('Fig 4.5'));
export const OPENING = { ...PETALS, span: [PETALS.span[1] - 6 * PI, PETALS.span[1]] };
