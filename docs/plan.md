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
| Momenta | A square pad in the controls rail: Φ across and Ψ up, on fixed axes from −3 to 3, so the region visibly reshapes as the spins change. The allowed region is tinted by the orbit's radius, darker for orbits closer in, with a key in the pad's corner (Mark, 2026-11-20), and the current orbit is a dot to drag. With spin, the halves are labelled "with a" and "against a", and "with b" and "against b". Below the pad, sliders for Φ and Ψ give fine and keyboard control; each stops at the region's edge. When a or b changes, the dot keeps its direction from the origin and its fraction of the way to the edge, so the orbit keeps its character as the region changes shape. |
| Trail | A "Trail length" slider in λ, on a log scale from about one oscillation (λ ≈ 3) to about 100 (λ ≈ 300). The trail fades towards its tail. Its length is the visitor's: presets leave it alone (Mark, 2026-11-20). |
| Animation | Play/Pause and Reset come from the kit. The photon moves at a steady rate in λ, about one oscillation every 1.5 s. A change to any control keeps the photon's place in its oscillation and redraws its past path for the new orbit at once, so the trail is never empty. The animation starts paused under reduced motion. |
| Presets | A drop-down: the paper's Figs 4.2, 4.3, 4.4 and 4.5 (one preset: the same orbit), 4.6, 4.7, 4.8 and 4.10 (4.9 and 4.11 plot the same orbits against λ) and 4.12. Also "Equal spins (a = b)", where every orbit's radius depends only on Φ + Ψ, and "No spin", the Tangherlini case, where every orbit is a great circle with r = √2. Each preset sets a, b, Φ, Ψ, the starting point, the colours and the direction the camera looks from, fitted to the figure's image ([research/scripts/presets/](research/scripts/presets/)). The trail's length and the zoom stay as the visitor left them. Any change of the visitor's own shows "Custom", and measures ψ, φ and λ from 0 at θ_min, so an orbit with Ψ = 0 around a hole with b = 0 is red, as in pictures of Kerr orbits (Mark, 2026-11-20). Reset puts the photon back where the last preset started it. |
| Opening state | Fig 4.5's orbit, moment, colours and viewpoint, with three oscillations of trail rather than the figure's 27 (Mark, 2026-11-20). The drop-down shows Figs 4.4 and 4.5's preset. |
| 3D view | Drag to rotate. Zoom with buttons and pinch, while the mouse wheel keeps scrolling the page. The zoom is the visitor's, and stays put through presets and every control, so an orbit's size on screen is its size (Mark, 2026-11-20). The default zoom fits every orbit's spheroid, and a "Reset view" button goes back to it and to the last preset's viewpoint. An arrow below the spheroid shows which way the hole turns in the a-plane: towards −φ, clockwise seen from above. It is below because most paths are on the upper half (Mark, 2026-11-20), and is hidden when a = 0. The b-plane's turning shows only in ψ, so its arrow goes round the colour wheel, towards −ψ, and is hidden when b = 0. |
| Side panel | θ, φ and ψ against λ over the visible trail, as in Figs 4.9 and 4.11: three rows sharing the λ axis, with the photon marked at the right-hand end. Dashed lines mark the band's edges, θ₋ and θ₊. φ and ψ run from 0 to 2π, as on the torus, and wrap at the edges, so no row's axis shifts (Mark, 2026-11-20). A band under ψ shows the trail's colour at each λ. Through an axis, θ, φ and ψ run on smoothly in the 3D view's signed coordinates (step 4). A second tab shows the (φ, ψ) torus: a square in which the path wraps around. It shows ψ without colour and makes quasiperiodicity visible. |
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

      **Measured in step 2:** against the drawn spheroid, pure blue is 1.6–1.7:1, below WCAG's 3:1 for lines; red is 3.4:1, magenta 4.3:1, and green, cyan and yellow 10:1 or more. The blue strands still stand out by hue. Mark looked at a side-by-side comparison and decided on no light edge (2026-11-20).
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
  src/scene/           three.js: spheroids, path, photon, spin arrow, camera
  src/pad.js           the momentum pad, a control with the kit's shape ({ el, value, onChange })
  src/momenta.js       the pad's geometry: where the region ends, and where the dot goes when the spins change
  src/spins.js         the spin sliders, held to a + b ≤ 1, with the "Extremal" label
  src/presets.js       the presets and the opening state
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

  The spheroid is drawn as two shells either side of the path: an opaque backdrop just outside, for the far half, and a translucent shaded shell just inside, over the near half. So nothing behind the spheroid shows through it, and it looks the same in both themes.

  The trail is drawn in time order, in two passes: the steps on the spheroid's far face, then those on its near face. Each step draws its dark edge and its colour together, so where strands cross, the newer one passes cleanly over the older, edge and all. Crossing strands lie at the same depth, so a depth test between them gave jagged crossings with no edge between the strands; and an edge drawn as a separate, wider line made the trail look dashed at every joint. The near/far split comes from the spheroid's outward normal at each step, so it costs nothing extra: the same two passes as before, and the same work per frame (measured: 2.7 ms before, 2.65 ms after, on the longest trail).
- [x] Temporary sliders for a, b, Φ and Ψ, and the "Trail length" slider. The opening state is Fig 4.10 until step 3 chooses one.
- [x] Camera controls, both themes, fitting the window, `needsWebGL`. Drag, or the arrow keys, to turn; the buttons, pinching, or + and − to zoom; "Reset view" goes back to the paper's viewpoint, fitted to the current spheroid. The zoom stays put when the orbit changes, so the spheroid visibly grows and shrinks.
- [x] Check in Edge (headless, through Playwright), with no console errors:
  - standalone, at phone width, and embedded as the site does it, in both themes;
  - Figs 4.2, 4.3, 4.4, 4.6, 4.7, 4.8 and 4.10, Tangherlini, near the extremal edge, and momenta with no orbit. Every readout matches the paper review's values;
  - `npm run build` and the validator: valid;
  - the bundle's size: `demo.js` 599 KB, the whole bundle 677 KB of the 2 MB;
  - the longest trail plays at 60 frames a second.
- [x] **Phones, early:** with emulated touch, a swipe on the 3D view turns it and doesn't scroll the page; a swipe anywhere else scrolls. On phones the 3D view now comes first, above the controls, and is square, at most 70% of the window's height, so there is always page to scroll by. Still to try on a real phone.

### 3. Momentum pad and presets (PR, then release `v0.3.0`)
- [x] The pad (`demo/src/pad.js`, its geometry in `demo/src/momenta.js`):
  - the region by ray bisection, tinted by the orbit's radius;
  - direction-and-fraction memory when the spins change;
  - snapping to the Φ = 0 and Ψ = 0 axes, within 6 pixels;
  - linked Φ and Ψ sliders, and arrow keys to move the dot (Shift for steps of 0.1). Each stops at the edge, and the status line says why.

  A point dragged outside the region goes to the edge along its ray from the origin. At an extremal hole the dot stops at 1 + Φ + Ψ = 10⁻⁴, short of the edge, where φ winds about 21 turns per oscillation. In a short window the pad shrinks, so the rail still fits beside the stage: it is 161 px across at 1366×657.
- [x] Presets for every figure, fitted to the paper's images ([research/scripts/presets/](research/scripts/presets/)) and checked side by side with them in Edge. They match strand for strand, colours included. What the fits found:
  - **The twist is real.** For Fig 4.10, the camera that matches the figure's box axes is the fitted one only once the arctan(a/r) twist is allowed for.
  - **Fig 4.8's colours are about π from Fig 4.9's ψ curve,** so that preset follows Fig 4.8.
  - **Fig 4.3 draws its sphere smaller than the orbit's spheroid,** so the red net bulges past the outline at the equator.
  - **Figs 4.4 and 4.5 are separate presets.** They are the same orbit, but different stretches seen from different sides.
- [x] The opening state: Fig 4.5 with three oscillations of trail (Mark chose it from pictures, 2026-11-20).
- [x] The extremal label, with the status line's reason. The readout and the out-of-region message were done in milestone 2.
- [x] The colour legend ring, top left of the 3D view, with a marker at the photon's ψ.
- [x] Checked in Edge, with no console errors: standalone and embedded with `--dk-fit-height`, both themes, phone width; dragging, snapping, the sliders' stops, the keys, the spins' memory, the extremal hold, presets and Reset. On a phone, a drag on the pad moves the dot and doesn't scroll the page. Moving a spin slider redraws the pad in about 16 ms; the longest trail near the extremal edge plays at 60 frames a second.

### 4. Side panel (PR, then release `v0.4.0`)
- [x] The θ, φ, ψ strip against λ (`demo/src/panel/strip.js`):
  - three rows sharing the λ axis, over the visible trail, with the photon at the right-hand end;
  - θ on a fixed scale, 0 to π/2, with dashed lines at the band's edges θ₋ and θ₊;
  - φ and ψ unwrapped, each row fitted to its stretch, with plain numbers on the axes as in the paper;
  - a band under ψ showing the trail's colour at each λ.

  The demo's λ is 0 at θ_min, and its φ is 0 there too. The axes read λ − λ₀, φ + φ₀ and ψ + ψ₀, where λ₀ and φ₀ come from the presets of Figs 4.8 and 4.10, fitted to Figs 4.9 and 4.11, and are 0 for the rest. The 3D view's φ₀ is already in its camera. So the Fig 4.10 preset's panel reproduces Fig 4.11: λ from 1 to 5, the same curves.
- [x] The (φ, ψ) torus tab (`demo/src/panel/torus.js`). The trail wraps at the edges and fades towards its tail, as the 3D trail does, and the colour wheel runs along the ψ axis.
- [x] **Through an axis, the curves stay smooth** (Mark, 2026-11-20: "The curves should be smooth, as they are geodesics. Jumps are purely artifacts of coordinate systems"). The panel uses the 3D view's signed coordinates (equations §7):
  - through the pole, θ swings between −θ₊ and θ₊;
  - through the equator, θ runs on past π/2 to π − θ₋;
  - through both (Φ = Ψ = 0, Figs 4.2 and 4.6), θ = χ/2 and grows steadily.

  φ and ψ have no jumps, and ψ always agrees with the trail's colours. The alternative, θ folded into [0, π/2] with φ or ψ jumping by π at each passage, was shown to Mark as pictures and set aside.
- [x] **Fig 4.8 keeps ψ₀ = 3.683,** so its colours match Fig 4.8, and its panel's ψ is Fig 4.9's curve raised by π. θ and φ match Fig 4.9. (Mark left this to Claude.)
- [x] **A preset's trail is its exact stretch.** The trail slider rounds the length to its step, 0.01 in log λ, which made Fig 4.10's start λ = 1.02. The exact length now holds until the slider moves.
- [x] **Layout.** On wide frames the panel sits beside the 3D view and takes 38% of the stage's width, from 15 to 24rem. On phones it goes below the 3D view, 20rem high. The rail is unchanged, so the frame still fits at 1366×657, with the pad at 161 px.
- [x] **Performance.** The panel walks the refined path through `trace()` in `trail.js`, which shares `fill()`'s window and per-oscillation shifts; nothing is solved again, and nothing is allocated per frame.
  - Each curve is drawn through at most four values per pixel column.
  - A trail of more than 20,000 points is visited at an even stride. Near the extremal edge a trail of λ = 300 has about 240,000 (1,024 a period); Fig 4.10's has 18,000, so it is drawn in full.
  - Measured in Edge, playing at 60 frames a second, the main thread's time per frame, 3D trail included: on Fig 4.10's longest trail, 4.3 ms with the strip and 4.6 ms with the torus; near the extremal edge, 8.8 ms and 10.5 ms.
- [x] Checked in Edge, with no console errors: standalone and embedded, both themes, phone width; every preset; the tabs' arrow keys; Reset. On a phone, a swipe on the panel scrolls the page, and one on the 3D view turns it.

### After step 4: notes on the live demo (PRs, then releases `v0.4.1` to `v0.4.3`)
- [x] **A steady zoom** (Mark, 2026-11-20). A preset used to fit the camera to its spheroid, so the zoom jumped with every preset. Now a preset turns the camera to its viewpoint and keeps the visitor's zoom. The default zoom fits a sphere of radius 2.25, which holds the largest spheroid any spins and momenta give: 2.23, counter-rotating at a or b = 0.99. "Reset view" goes back to it.
- [x] **The trail's length is the visitor's** (Mark, 2026-11-20). Presets no longer set it, and changing it no longer shows "Custom". The opening still has exactly three oscillations of Fig 4.5's trail. A preset's span now gives only the photon's start.
- [x] **Which way the hole turns** (Mark, 2026-11-20: "it's not clear in which direction the black hole is rotating"). An arrow below the spheroid, around a short axle on the z axis, turns clockwise seen from above, towards −φ (`demo/src/scene/spin.js`).
  - It was above the spheroid in v0.4.1, which Mark found too obtrusive: most paths are on the upper half. In v0.4.2 it is below, smaller and a muted grey-blue.
  - It hangs 0.6 below the spheroid, which clears the spheroid's outline from every preset's viewpoint (they look down at up to 36°; the Equal spins preset needs 0.6, the rest less). From higher up, the spheroid hides it.
  - It has a dark outline, so it reads on the dark spheroid and on either theme's page.
  - It turns with the camera, so its head is always at the front left: from the side, the front of the hole moves to the left.
  - It is hidden when a = 0. Checked against the physics: an orbit with Φ = 0 drifts towards −φ, by −0.13 turns per oscillation for a = 0.5, b = 0.
- [x] **Which way the hole turns in the b-plane** (Mark chose it, 2026-11-20). The 3D view shows (x₁, x₂) and only the distance from the axis in the b-plane, so the b-rotation moves nothing in it: it is all in ψ, the colour. So b's arrow goes round the colour wheel, which is the b-plane, towards −ψ: clockwise there. It looks like a's arrow and is hidden when b = 0. Fig 4.6's photon, with Ψ = 0, drifts the way it points.
- [x] **A custom orbit's colours start from red** (Mark, 2026-11-20: the origin of ψ is arbitrary, and red matches pictures of Kerr orbits). A change to the spins or momenta sets ψ₀, φ₀ and λ₀ to 0, so ψ, φ and λ are 0 at θ_min, and an orbit with Ψ = 0 around a hole with b = 0, where ψ stays put, is red all along. Before, a custom orbit kept the last preset's offsets: from the opening, Ψ = 0 gave a green orbit. The colours change once, as the state leaves the preset. v0.4.3.
- [x] **Presets that only differed in length** (Mark, 2026-11-20). Figs 4.4 and 4.5 are the same orbit, so they share one preset, "Figs 4.4, 4.5", with Fig 4.5's viewpoint and colours, which the opening also uses. Fig 4.7 is "two spins" rather than "two spins, a band": the band shows only with a long trail. v0.4.3.
- [x] **The side panel's φ and ψ run from 0 to 2π,** as on the torus (Mark, 2026-11-20: the shifting axes were confusing). This replaces step 4's unwrapped rows, each fitted to its stretch. Where a curve passes a multiple of 2π, it leaves at one edge of its row and comes back at the other. θ does too through both the pole and the equator, where it grows by π each oscillation, so every row's axis is fixed. The ticks are fractions of π on every row.
- [x] Checked in Edge, with no console errors: both themes; Figs 4.2, 4.3, 4.5, 4.6, 4.7, 4.8, 4.10 and 4.12; zooming in, then choosing a preset; "Reset view"; the longest trail, which keeps its length through presets.

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

Settled with Mark on 2026-11-20. The comparison page, <https://claude.ai/artifact/H3VwXgAEs8gpsK1A31fqcv>, which is private to Mark, showed him the colour wheels and the Ψ = 0 pictures. For step 4, <https://claude.ai/artifact/5Nqu8jw1oocXfkeszj2eK8>, also private, showed the side panel's two choices.

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
| Blue strands | No light edge: blue is dim against the spheroid (1.6–1.7:1) but reads by hue (step 2). |
| Opening state | Fig 4.5's orbit and viewpoint, with three oscillations of trail (step 3). |
| The pad's region | Tinted by the orbit's radius (step 3). |
| The side panel through an axis | Smooth, in the 3D view's signed coordinates: "jumps are purely artifacts of coordinate systems" (step 4). |
| Fig 4.8's colour offset | Follows Fig 4.8's colours rather than Fig 4.9's ψ curve; Mark left it to Claude (step 4). |
| Zoom and trail length | The visitor's own: presets leave them alone (after step 4). |
| Which way the hole turns | An arrow below the spheroid for a, where it is out of the paths' way, and one round the colour wheel for b, whose turning shows only in ψ (after step 4). |
| The side panel's φ and ψ | From 0 to 2π, as on the torus, rather than unwrapped (after step 4). |
| A custom orbit's colours | ψ = 0, red, at θ_min, so Ψ = 0 orbits at b = 0 are red, as in pictures of Kerr orbits (after step 4). |

## Risks

- **Phones.** Dragging to rotate competes with scrolling the page. Tested in step 2 with emulated touch: the page scrolls everywhere but on the 3D view, which is at most 70% of the window's height. Still to try on a real phone. The fallback is two-finger rotation on touch screens.
- **The pad on narrow screens.** On phones, the controls rail stacks below the stage (since step 2), and the pad is as wide as the rail, up to 320 px. Beside the stage on a short desktop window it shrinks, down to 144 px.
- **Orbits near an axis.** φ or ψ can swing by about π within one step, which would smear the colour. Subdivision handles it, and Fig 4.10 is the test case.
- **Orbits near the extremal co-rotating edge** wind φ and ψ by tens of turns per oscillation: about 21 at a = b = 0.5, 1 + Φ + Ψ = 10⁻⁴. The path must subdivide by φ as well as ψ. The pad stops the dot at 1 + Φ + Ψ = 10⁻⁴ (step 3).
- **The one-orbit result** is numerical. The solver is written defensively, and logs if it ever finds a second root.
- **Size.** About 0.6 MB of the 2 MB budget, so there is room.
