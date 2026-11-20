# Plan: Plotting the Photons

This is the plan for building this project's interactive demo for mark-bugden.github.io. Claude wrote it on 2026-11-20 from two research reports, and Mark settled its decisions the same day (see "Decisions" at the end). The site's requirements are in `site-contract-photon-orbits-5d.md`, which lives outside this repo, and in the [demo-kit README](https://github.com/Mark-Bugden/demo-kit#readme). This document covers what we build on top of them.

The research behind it:

- [`research/paper-review.md`](research/paper-review.md): a line-by-line check of arXiv:1801.03389. It has the errata, the formulas in the form the demo needs, verified presets for every figure, and unit-test values.
- [`research/prior-art.md`](research/prior-art.md): what already exists (Leo Stein's Kerr tool and others) and how it is built. It also has measured bundle sizes, ways to draw the path, and a comparison of colour maps.
- [`research/scripts/`](research/scripts/): the Python checks behind the paper review, including a reference implementation (`mp5d.py`, `reference_impl.py`).

Where the reports and this plan differ, this plan wins. For example, the paper review recommends coordinate time as the clock, but the demo uses λ.

## The goal

The demo shows spherical photon orbits around a five-dimensional Myers–Perry black hole, after Mark's paper, in the spirit of Leo Stein's [Kerr spherical photon orbits](https://duetosymmetry.com/tool/kerr-circular-photon-orbits/) tool. A visitor:

1. sets the black hole's two spins, a and b, and watches the horizon and the orbit's spheroid change shape: oblate, a sphere, or prolate
2. picks the photon's two angular momenta, Φ and Ψ, on a pad that shows which pairs have a spherical orbit, where each point is exactly one orbit
3. watches the photon trace its path over the spheroid, with the fourth spatial direction shown as colour, and chooses how much of its past path stays drawn
4. jumps to presets that reproduce the paper's Figures 4.2 to 4.12.

Nothing like it exists for five dimensions. Of the 9 works that cite the paper, none redraws 5D orbits, and no Myers–Perry orbit tool is interactive ([prior art §3](research/prior-art.md#3-myersperry-and-higher-dimensional-black-holes)).

## What changes from the brief

The portfolio brief already sketches this demo. These things change:

| Brief | Plan | Why |
| --- | --- | --- |
| Working title "Circling the Sphere" | **Plotting the Photons**. The slug stays `photon-orbits-5d`. | Mark's choice (2026-11-20). |
| Controls Φ and Ψ, plus a "find spherical orbit" solver for the radius | Φ and Ψ on a 2D pad that draws the region where orbits exist, with two sliders for fine control. The solver runs on every change. | Outside a convex "curved diamond" in the (Φ, Ψ) plane, no spherical orbit exists. Inside it, there is exactly one, so with two plain sliders most combinations would just fail. The diamond's shape changes with a and b. |
| Adaptive RK45 in a Web Worker | No ODE solver and no worker. One latitude oscillation is computed by quadrature, then repeated. | At fixed r the motion reduces to integrals over one oscillation. Measured: accurate to 1e-7 even next to the axes, in microseconds. Stein-style fixed-step RK4 is off by 1.8 rad in ψ on Fig 4.10. |
| "The pseudo-Cartesian projection from the Kerr–Schild form" | The paper's x₁, y₁, z₁ without the constant arctan(a/r) twist | It is the same picture. It is the exact flat-space embedding of ℝ⁴, and at fixed r the offset is a constant rotation. |

Unchanged from the brief:
- the paper's colour wheel;
- θ, φ, ψ plotted against λ;
- presets for Figs 4.2–4.12;
- three.js;
- Mark reviewing the equations before any drawing.

## The physics

The full derivations, checks and special cases are in [paper review, Part B](research/paper-review.md#part-b--what-the-implementer-needs). In short:

- **Units and normalisation.** μ = 1 and E = 1, so Φ, Ψ mean Φ/E, Ψ/E, and K means K/E².
- **Coordinates.** x = r², y = sin²θ, θ ∈ [0, π/2].
- **The metric** is eq. (2.1) with its typo corrected (see "Errata in the paper" below).
- **Signs.** For a > 0 the hole turns towards −φ, so co-rotating photons have Φ < 0. The same holds for b and Ψ.

```
Horizon      x₊ = ½(1 − a² − b² + √((1 − a² − b²)² − 4a²b²));  a horizon exists iff a + b ≤ 1  (a, b ≥ 0)
Radial       𝒳(x) = P(x) − K·Δ(x),   Δ(x) = x² + (a² + b² − 1)x + a²b²
             P(x) = x³ + (a² + b²)x² + [(a² − b²)(Φ² − Ψ²) + a²b² + a² + b² + 2aΦ + 2bΨ]x
                    + (ab + bΦ + aΨ)² + (a² − b²)(b²Φ² − a²Ψ²)
Orbit        x₀ = the root of F = P′Δ − PΔ′ above x₊ (bisection);  K = P(x₀)/Δ(x₀)
Latitude     (dy/dσ)² = 4Q(y),  Q(y) = (a² − b²)y³ + (b² − 2a² − K)y² + (Φ² − Ψ² + a² + K)y − Φ²
             an orbit exists iff Q ≥ 0 somewhere in (0, 1); it oscillates in the band [y₁, y₂]
Angle        y = ȳ − Δy·cos χ  (ȳ, Δy: the band's middle and half-width);  no turning points in χ
             dσ/dχ = 1/(2√G(y)),   G(y) = (b² − a²)(y + y₁ + y₂) − (b² − 2a² − K) > 0
Motion       d(σ, λ, φ, ψ)/dχ = (1, ρ², Φ/y − A_φ, Ψ/(1 − y) − A_ψ) / (2√G)
             ρ² = x + a² + (b² − a²)y;  A_φ, A_ψ are constants of the orbit (paper review §B3)
Drawing      X = √(x + a²) sinθ cos φ,   Y = √(x + a²) sinθ sin φ,   Z = √(x + b²) cosθ,   colour = ψ
```

- **One orbit per (Φ, Ψ).** In 200,000 random samples and 16 dense grids, F had exactly one root above the horizon. So the momenta fix the radius and K whenever an orbit exists. This is numerical evidence, not a proof. The solver still checks for a second root and logs a warning if it ever finds one.
- **The allowed region 𝒜(a, b)** is a convex "curved diamond" around the origin.
  - Its corners on the axes are the paper's equatorial and polar circular orbits, eqs (3.9)–(3.15).
  - Its edges are circular orbits of constant θ; Fig 4.12 is one of them.
  - The boundary is found by bisection along 96 rays from the origin. That works because the region is star-shaped.
- **Every orbit is radially unstable.** Nudged, r grows by about e⁴ per oscillation. A full integration that includes r leaves the sphere within about 6 oscillations, at any tolerance. So the demo never integrates r: it holds r fixed, as the paper does, and "How this works" says so.
- **One oscillation (χ from 0 to 2π)** is tabulated with Simpson's rule, at N = 256 points, or 512 when the band comes within 10⁻⁴ of an axis. The 1/y and 1/(1 − y) spikes near the axes are subtracted and integrated in closed form (paper review §B3, step 4).
- **Longer paths repeat the table,** shifting φ, ψ and λ by their advance per oscillation. A long trail costs nothing extra and never drifts.
- **The clock is λ, the affine parameter** (Mark, 2026-11-20), as in the paper. One parameter is used throughout: the animation, the trail length and the side panel. Coordinate time t isn't shown. Its formula is in the paper review, if it's ever needed.
- **Special cases are drawn the paper's way** (Mark, 2026-11-20). Both use signed coordinates (paper review §B3, "Special handling"), as the paper's figures do.
  - **Φ = 0:** the orbit passes straight through the pole. A nonzero Φ keeps it away from the pole.
  - **Ψ = 0:** the orbit crosses the equator onto the lower half of the spheroid and keeps its colour. A nonzero Ψ keeps it away from the equator: it turns back just short of the equator, with a quick change of colour.

  So the picture changes at exactly Ψ = 0, by design, because Ψ = 0 is a special case. The pad snaps to both axes, so the special cases are easy to land on.
- **Exact extremality, a + b = 1, is allowed** (Mark, 2026-11-20); Fig 4.8 uses a = 0.9, b = 0.1. At exact extremality the co-rotating corner of the region touches the horizon, so near that corner the solver can meet 0/0. Two guards make it safe:
  - the boundary search treats a 0/0 at the horizon as "no orbit here";
  - the solver's bracket starts just above x₊, and square roots are clamped at 0.

  With them, a stress test passed every time. It used 96 directions × 8 fractions of the way to the edge, on a = 0.9, b = 0.1; a = 0.99, b = 0.01; a = 0.7, b = 0.3; a = b = 0.5; and three cases with a + b = 0.999. Even at 0.999 of the way to the edge, the orbits stay at least 0.002 in x outside the horizon. a = 1, b = 0 (x₊ = 0) stays excluded by the sliders' 0.99 maximum.

  **Update from milestone 1** ([equations §8](equations.md#extremal-holes)):
  - It is not just the corner. At exact extremality the whole co-rotating edge is the straight segment Φ + Ψ = −1, where 𝓔(x₊) = 0, and the orbits close on the horizon along it.
  - The research solver stopped about 3·10⁻⁴ short of that edge. The solver now works in u = x − x₊, which reaches it exactly.
  - With the true edge, the stress test's orbits stay at least 0.0019 in x outside the horizon.
  - Near that edge, φ winds tens of turns per oscillation.
- **Near-zero momenta are snapped to 0.** A |Φ| or |Ψ| below 10⁻¹² counts as exactly 0. Otherwise a value like Φ = 10⁻¹⁶ (for example, from the ray straight up the Ψ axis) puts the band edge at y₁ ≈ 0 while the formulas still divide by y. The reference script `reference_impl.py` divides by zero there, at any spins.

### Errata in the paper

The demo uses the corrected forms. Details are in [paper review §A4](research/paper-review.md#a4-list-of-mistakes-and-issues).

1. **Eq. (2.1):** `(x+b)² cos²θ dψ²` should be `(x + b²) cos²θ dψ²`.
2. **The mass formula (p. 2):** `M = 3μ/(8√π G)` should be `M = 3πμ/(8G)`.
3. **§3.1.1, the a = b case:** the maximum of Θ is a² + K − (|Φ| + |Ψ|)². The paper's (Φ + Ψ)² is wrong when Φ and Ψ have opposite signs.
4. **§3.4:** the constant term of the equatorial cubic repeats the x¹ coefficient. It should be b²[(Φ + a)² + b²(a² − Φ²)]. The solutions (3.9)–(3.15) are still right.
5. **Figs 4.9 and 4.11:** the curve labelled "x" is r = √x.

Everything else checks out exactly, and the figure parameters reproduce every figure.

## The experience

| | |
| --- | --- |
| Spins | Sliders for a and b, each from 0 to 0.99, held to a + b ≤ 1: beyond that there is no horizon. When a slider stops, the status line says why. Negative spins only mirror the picture, so they're left out. |
| Extremal black holes | When a + b = 1, an "Extremal" label appears beside the spin sliders, and the status line says why it is special: the inner and outer horizons have merged, and the hole spins as fast as a black hole of its mass can. The Fig 4.8 preset is extremal. |
| Momenta | A square pad in the controls rail: Φ across and Ψ up, on fixed axes from −3 to 3, so the region visibly reshapes as the spins change. The allowed region is shaded, optionally tinted by the orbit's radius, and the current orbit is a dot to drag. With spin, the halves are labelled "with the spin" and "against it". Below the pad, sliders for Φ and Ψ give fine and keyboard control; each stops at the region's edge. When a or b changes, the dot keeps its direction from the origin and its fraction of the way to the edge, so the orbit keeps its character as the region changes shape. |
| Trail | A "Trail length" slider in λ, on a log scale from about one oscillation (λ ≈ 3) to about 100 (λ ≈ 300). The trail fades towards its tail. |
| Animation | Play/Pause and Reset come from the kit. The photon moves at a steady rate in λ, about one oscillation every 1.5 s. A change to any control keeps the photon's place in its oscillation and redraws its past path for the new orbit at once, so the trail is never empty. The animation starts paused under reduced motion. |
| Presets | The paper's figures: 4.2, 4.3, 4.4/4.5, 4.6, 4.7, 4.8/4.9, 4.10/4.11 and 4.12. Also "Equal spins (a = b)", where every orbit's radius depends only on Φ + Ψ, and "No spin", the Tangherlini case, where every orbit is a great circle with r = √2. Each preset sets a, b, Φ, Ψ, the trail length, the starting point and the camera, from [paper review §B5](research/paper-review.md#b5-presets-reproducing-figs-42412). |
| Opening state | One of the presets, chosen in step 3 for the most striking first view. Candidates: Fig 4.10's flower, or Fig 4.4/4.5's petals. |
| 3D view | Drag to rotate. Zoom with buttons and pinch, while the mouse wheel keeps scrolling the page. A "Reset view" button. |
| Side panel | θ, φ and ψ against λ over the visible trail, as in Figs 4.9 and 4.11. A cursor marks the photon, and dashed lines mark the band's edges, θ₋ and θ₊. φ and ψ are unwrapped, as in the paper. A second tab shows the (φ, ψ) torus: a square in which the path wraps around. It shows ψ without colour and makes quasiperiodicity visible. |
| Readouts | The status line gives r, the θ band and the advance per oscillation in turns, for example "r = 1.32 · θ from 37.6° to 51.5° · per swing: φ +0.34 turns, ψ −0.61 turns". If the dot ever lands outside the region, it says "No spherical orbit for these momenta". |
| Copy | Plain. "How this works" is drafted for Mark to rewrite: on the site, Mark writes the English copy himself. |

## Look

- **The paper's look.**
  - A dark orbit spheroid with a pale equator ring, and the black horizon spheroid inside it.
  - The coloured path lies on the spheroid's surface.
  - The spheroid is dark in both themes, so the path's colours read the same in both. In the light theme it sits on the white page like the paper's figures.
- **Depth.** The spheroid is slightly translucent, so the far side of the path shows through, dimmed. Faint fog completes the depth cue.
- **The path.**
  - A thick, anti-aliased line, unlit, so that its lightness carries only ψ.
  - A thin dark casing around it, so it stays legible where it crosses the spheroid's silhouette.
  - A small bright dot marks the photon.
- **The colour wheel is the paper's** (Mark, 2026-11-20): hue = ψ/2π at full saturation and value, as Mathematica's `Hue` in Fig 4.1. That is red at 0, yellow at π/3, green at 2π/3, cyan at π, blue at 4π/3 and magenta at 5π/3. It is stored as a 256-entry table.
  - **Known limits** (measured for the comparison page, <https://claude.ai/artifact/H3VwXgAEs8gpsK1A31fqcv>):
    - about 12% of the wheel, the deep blue around ψ = 4π/3, is dim against the black sphere. Check in step 2 whether the path needs a faint light edge there.

      **Measured in step 2:** against the drawn spheroid, pure blue is 1.6–1.7:1, below WCAG's 3:1 for lines; red is 3.4:1, magenta 4.3:1, and green, cyan and yellow 10:1 or more. The blue strands still stand out by hue, so step 2 adds no light edge. Mark to decide.
    - Opposite angles are less distinct for red-green colour-blind viewers. The torus tab shows ψ by position, so colour is never the only way to read it.
- **Colour legend.** A ring in a corner of the 3D view, labelled 0, π/2, π and 3π/2, with a needle at the photon's current ψ.
- **The page's colours** come from the kit v0.5.0 tokens.

## Tech

- **demo-kit v0.5.0,** the latest, which brings dark-by-default and the site's current palette (Mark, 2026-11-20). The change (demo version 0.1.2):
  - the tag in `demo/package.json`;
  - `kitVersion` in the manifest;
  - `<meta name="color-scheme" content="dark">` in `index.html`, as in the kit's template.
- **three.js r187, with named imports only:**
  - `WebGLRenderer`;
  - `OrbitControls`;
  - `Line2`, `LineMaterial` and `LineGeometry` for the thick path;
  - plain meshes for the spheroids.

  This is about 580 KB minified, of the 2 MB budget.
- **No Observable Plot.** The side panel, pad and legend are drawn on 2D canvases, which saves about 275 KB. There are no other runtime dependencies.
- **Manifest.** Set `needsWebGL: true`. If WebGL fails to start, the still image stays and the status line says why.
- **Performance.**
  - **Main thread only.** Recomputing an orbit costs well under a millisecond, and redrawing the pad's region a few milliseconds.
  - **Changes.** At most one recompute per frame, and nothing allocated per frame.
  - **The path's buffer** is allocated once, for 65,536 segments. Each frame rewrites only the visible window. A trail of λ = 300 on Fig 4.10 is about 18,000 segments, filled in about 2 ms. If a trail ever needs more, near the extremal edge, its oldest part is left out.
  - **Detail near the axes.** One oscillation is drawn in 128 steps, and any step in which φ or ψ turns by more than 10° is halved, as often as needed, using the state between table nodes ([equations §6](equations.md#6-one-oscillation)).
  - **Pixel ratio** is capped at 2.
- **Reuse.** The licence of Stein's page code is unclear, so we reimplement rather than copy. Ricky Reusser's "Geodesics in Kerr Spacetime" notebook (MIT) is the bar for visual polish.
- **Two implementations.**
  - A Python reference in `src/photon_orbits/`, built from the research scripts.
  - The JS port in `demo/src/physics/`.
  - JS tests replay fixtures that the Python writes, as in zombie-outbreak.

## Repo layout (target)

```
demo/                  the site bundle, built with demo-kit (see the site contract)
  src/main.js          frame, controls and wiring
  src/physics/         pure functions: horizon, solver, band, region, one-oscillation table, frequencies
  src/scene/           three.js: spheroids, path, photon, camera
  src/pad.js           the momentum pad, a control with the kit's shape ({ el, value, onChange })
  src/panel/           the side panel and the colour legend (2D canvas)
  src/colour.js        the paper's colour wheel
  src/trail.js         the path ready to draw: one oscillation refined, and the trail's window of it
  test/                node:test, including the fixture replay
src/photon_orbits/     Python reference model
tests/                 pytest: the review's checks, and the fixture writer
scripts/stills.py      renders static, thumb and og images from the running demo
docs/                  this plan, equations.md, research/
pyproject.toml         a uv project; commit uv.lock
```

## Milestones

### 0. Research and plan (PR)
- [x] Paper review, prior-art survey, and this plan.
- [x] Mark settled the decisions (2026-11-20).
- [x] Kit v0.5.0 and the new title. The stub's title card, `static.webp`, `thumb.webp` and `og.jpg`, now reads "Plotting the Photons", redrawn in the original's font. Demo version 0.1.2.
- [x] Mark merges the PR, which releases v0.1.2, still a stub.
- [ ] Mark updates the title in the site's write-up, `content/projects/photon-orbits-5d.md` in Mark-Bugden.github.io, which still says "Circling the Sphere".

### 1. Physics core (PR). Gate: Mark reviews the equations before any drawing
- [x] Set up the uv project, and turn `research/scripts/mp5d.py` and `reference_impl.py` into the Python reference module, with the two guards from "The physics". pytest covers:
  - the unit-test table, region vertices, out-of-region cases, invariances and the closed orbit ([paper review §B7](research/paper-review.md#b7-unit-test-values));
  - extremal black holes, and momenta of 10⁻¹⁶;
  - the Hamiltonian cross-check against the full metric.
- [x] `docs/equations.md`: the corrected equations exactly as implemented, for Mark's review.
- [x] The JS port in `demo/src/physics/`, with node:test replaying the Python fixtures to 1e-9.
- [x] A README for the repo.
- [x] Bump the version, since `demo/` changes (0.1.3). The site's stub doesn't change.
- [x] Mark reviews `docs/equations.md` (confirmed 2026-11-20). Milestone 2 adds one function, `between`, the state between two table nodes, for the drawing (equations §6, "Between the nodes").

### 2. The 3D view (PR, then release `v0.2.0` with status `in-progress`)
- [x] The three.js scene:
  - the spheroids and equator ring;
  - the path coloured by ψ, the photon, the trail, and the animation clock.

  The spheroid is drawn as two shells either side of the path: an opaque backdrop just outside, for the far half, and a translucent shaded shell just inside, over the near half. So nothing behind the spheroid shows through it, and it looks the same in both themes. The casing is drawn only where the coloured line isn't, using the stencil buffer. A depth offset left the line looking dashed, because each step's round end caps are flat in depth.
- [x] Temporary sliders for a, b, Φ and Ψ, and the "Trail length" slider. The opening state is Fig 4.10 until step 3 chooses one.
- [x] Camera controls, both themes, fitting the window, `needsWebGL`. Drag, or the arrow keys, to turn; the buttons, pinching, or + and − to zoom; "Reset view" goes back to the paper's viewpoint, fitted to the current spheroid. The zoom stays put when the orbit changes, so the spheroid visibly grows and shrinks.
- [x] Check in Edge (headless, through Playwright), with no console errors:
  - standalone, at phone width, and embedded as the site does it, in both themes;
  - Figs 4.2, 4.3, 4.4, 4.6, 4.7, 4.8 and 4.10, Tangherlini, near the extremal edge, and momenta with no orbit. Every readout matches the paper review's values;
  - `npm run build` and the validator: valid;
  - the bundle's size: `demo.js` 599 KB, the whole bundle 677 KB of the 2 MB;
  - the longest trail plays at 60 frames a second.
- [x] **Phones, early:** with emulated touch, a swipe on the 3D view turns it and doesn't scroll the page; a swipe anywhere else scrolls. On phones the 3D view now comes first, above the controls, and is square, at most 70% of the window's height, so there is always page to scroll by. Still to try on a real phone.

### 3. Momentum pad and presets (PR)
- [ ] The pad:
  - the region by ray bisection;
  - direction-and-fraction memory when the spins change;
  - snapping to the Φ = 0 and Ψ = 0 axes;
  - linked Φ and Ψ sliders, and arrow keys to move the dot.
- [ ] Presets for every figure. Check each side by side with the paper and with `research/scripts/figures_repro.png`.
- [ ] Choose the opening state.
- [ ] The status-line readouts, the out-of-region message, and the extremal label.
- [ ] The colour legend ring.

### 4. Side panel (PR)
- [ ] The θ, φ, ψ strip against λ. The (φ, ψ) torus tab.
- [ ] On wide frames, the panel sits beside the 3D view; on phones, below it or folded away.

### 5. Polish and release (PR, then release `v1.0.0`)
- [ ] **Accessibility:**
  - every control works from the keyboard;
  - a live text description of the orbit;
  - reduced motion;
  - colour-blind checks with DevTools emulation;
  - contrast.
- [ ] **Performance on a phone,** with no console errors.
- [ ] **Stills** from the demo, with `scripts/stills.py` and Playwright: `static`, `thumb` and `og`, plus the optional light versions.
- [ ] **Manifest:** summary, tags, `about` in Mark's wording, and `built`. Mark decides when the status becomes `live`.
- [ ] **Site:** Mark writes the write-up.

### Later (stretch, new work)
- **Resonance hunt.**
  - The advances per oscillation, ν_φ and ν_ψ, come free from the table.
  - An orbit closes only when both are rational, and damped Newton steps find such orbits ([paper review §B6](research/paper-review.md#b6-frequencies-and-resonances-stretch-goal)).
  - On the torus tab, a closed orbit is a closed curve.
- **"Swap planes."** Draw the b-plane as the horizontal and colour by φ. This shows the a ↔ b symmetry: Figs 4.2 and 4.6 are the same kind of orbit.
- **A stereographic view** of the 3-sphere, as an advanced view.
- **"Nudge the photon."** Integrate the full equations from a slightly perturbed start, to show the instability.

## Decisions

Settled with Mark on 2026-11-20. The comparison page, <https://claude.ai/artifact/H3VwXgAEs8gpsK1A31fqcv>, which is private to Mark, showed him the colour wheels and the Ψ = 0 pictures.

| Question | Decision |
| --- | --- |
| Title | "Plotting the Photons". The slug stays `photon-orbits-5d`. |
| Choosing the momenta | The pad, with two fine sliders under it. |
| Colour wheel | The paper's wheel, as in Fig 4.1. |
| Special cases | The paper's way: Φ = 0 orbits pass through the pole, and Ψ = 0 orbits cross the equator. |
| Extremal spins | Allowed (a + b = 1), with an "Extremal" label. |
| Clock | λ only: the animation, the trail length and the side panel all use it. |
| demo-kit | v0.5.0. |
| Errata | Whether to post a corrected arXiv version is decided later. |

## Risks

- **Phones.** Dragging to rotate competes with scrolling the page. Tested in step 2 with emulated touch: the page scrolls everywhere but on the 3D view, which is at most 70% of the window's height. Still to try on a real phone. The fallback is two-finger rotation on touch screens.
- **The pad on narrow screens.** On phones, the controls rail stacks below the stage (since step 2), and the pad must stay usable at about 280 px.
- **Orbits near an axis.** φ or ψ can swing by about π within one step, which would smear the colour. Subdivision handles it, and Fig 4.10 is the test case.
- **Orbits near the extremal co-rotating edge** wind φ and ψ by tens of turns per oscillation: about 21 at a = b = 0.5, 1 + Φ + Ψ = 10⁻⁴. The path must subdivide by φ as well as ψ. The pad might stop the dot just short of that edge.
- **The one-orbit result** is numerical. The solver is written defensively, and logs if it ever finds a second root.
- **Size.** About 0.6 MB of the 2 MB budget, so there is room.
