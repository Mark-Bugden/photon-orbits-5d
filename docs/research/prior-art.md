# Prior art: interactive spherical photon orbits around a 5D Myers–Perry black hole

Research for the "Plotting the Photons" demo (`demo/`, slug `photon-orbits-5d`), based on M. Bugden, *Spherical photon orbits around a 5D Myers–Perry black hole*, Gen. Rel. Grav. 50, 30 (2018), [arXiv:1801.03389](https://arxiv.org/abs/1801.03389). Compiled 2026-11-20.

## How to read this report

Every claim about an existing tool carries one of these labels:

| Label | Meaning |
| --- | --- |
| **[V]** | Verified. We read the page, the source code, the paper or an API response. Some sections were gathered by research sub-agents, who read the sources themselves and are labelled the same way. |
| **[M]** | Measured. I ran it: bundle builds, numerical benchmarks or colour-map metrics. Appendix A gives the method and versions. |
| **[S]** | Seen only in a search-result snippet; the page itself was not opened. |
| **[I]** | Inferred. This is my reading, derivation or opinion. |

No code was copied into the demo. The few short snippets below illustrate a technique, and each one's licence is noted.

---

## 0. Summary

- **Nothing like this demo exists.** No interactive visualiser of Myers–Perry (or any higher-dimensional) geodesics or photon orbits turned up. GitHub has no "myers-perry geodesic" repositories, and none of the 9 papers citing arXiv:1801.03389 re-plots 5D orbits **[V]** (§2, §3). The closest relatives are:
  - Leo Stein's 4D Kerr tool (2016);
  - Ricky Reusser's WebGPU Kerr-geodesics notebook (2026);
  - static Mathematica figures in the 5D papers.
- **Stein's tool** has these parts, all verified from its source (§1):
  - JSXGraph sliders and fixed-step RK4;
  - a 2016 three.js (r74), set up through threestrap;
  - a single 1-pixel red `THREE.Line`, rebuilt every 35 ms by `setTimeout`;
  - a translucent oblate-spheroid horizon.

  Its numerics are clever: a monotonic phase χ in place of θ, plus Δφ from complete elliptic integrals. Its own page names the weak points: fixed-step RK4 "is a problem for orbits close to the spin axis", polar orbits do not work, and the code "is a bit of spaghetti", with memory leaks.
- **There is a better way to compute the orbit than an ODE solver.** At fixed r the orbit has two parts:
  - The latitude motion is periodic. In u = cos²θ and Mino time it is a cubic-potential problem: elliptic in general, and elementary when a = b.
  - φ and ψ are 1D integrals over that motion.

  So: compute one latitude period by quadrature, using a phase variable χ and an exact treatment of the 1/sin²θ and 1/cos²θ spikes near the axes. Then tile copies, rotated by Δφ and Δψ, to build any history length. On the paper's near-axis Fig. 4.10 orbit, this is accurate to better than 1e-7 rad per period with 64 panels and costs about 9 µs **[M]**.

  The same orbit integrated Stein-style (fixed-step RK4 in affine time, Δλ = 0.05) is off by **1.8 rad in ψ per period** **[M]** (§5.5). Every recompute takes under 2 ms, so no Web Worker is needed.
- **Recommended stack:**
  - three.js r187, imported tree-shaken;
  - `OrbitControls`;
  - `Line2`/`LineMaterial`/`LineGeometry` fat lines with per-vertex colour;
  - a preallocated buffer, with `instanceCount` as the history length.

  That adds about **580 KB minified, 141 KB gzipped** **[M]**, under a third of the 2 MB budget. Adding Observable Plot for the side panel brings it to about 855 KB. Draw the animated side panel on a plain 2D canvas instead (§5).
- **Colour for ψ.** Use a perceptually uniform, colour-blind-robust **cyclic** map with a ring-shaped legend. Do not use the HSV hue wheel used in the paper's figures.
  - **Suggested default:** Kovesi *CET-CBTC1*, drawn over a dark orbit spheroid as in the paper's figures, with a thin casing on the line. Alternatives are Crameri *romaO* and *CET-CBC1*.
  - Do not use "two-cycle" maps (CET-C4, C5, CBC2, CBTC2): they give ψ and ψ+π the same colour.
  - No single cyclic map stays above 3:1 contrast against both the light and the dark kit surface and remains distinguishable under colour-vision deficiency (CVD) **[M]**. Plan for a casing or a dark backdrop, and show ψ redundantly in a (φ, ψ) panel (§4).

---

## 1. Leo Stein's Kerr spherical photon orbits tool

**Page:** https://duetosymmetry.com/tool/kerr-circular-photon-orbits/

- **Title and dates [V]:** titled "Kerr Spherical Photon Orbits". It was posted 2016-03-02 and modified 2016-03-05, with a typo correction in 2019 noted on the page.
- **Source [V]:** a Jekyll post, `_posts/2016-02-23-kerr-circular-photon-orbits.md`, in https://github.com/duetosymmetry/web-site on the `customization` branch. The page JavaScript is inline in the post.
- I downloaded the page and every script it loads and read them; everything in this section is **[V]** unless marked otherwise.

### 1.1 Libraries

The page loads these:

| Script | What it is |
| --- | --- |
| `assets/js/vendor/three.min.js` | three.js **revision 74** (Jan 2016). Uses the legacy `THREE.Geometry` API, which was removed from core in r125. |
| `assets/js/vendor/threestrap.min.js` | [threestrap](https://github.com/unconed/threestrap) (Steven Wittens), a three.js bootstrapper. It creates the renderer, scene, camera, render loop and resize handling. |
| `assets/js/vendor/OrbitControls.js` | The old global-script `THREE.OrbitControls`, with touch support. |
| `assets/js/vendor/jsxgraphcore.js` | [JSXGraph](https://jsxgraph.org) 1.12.2, "dual licensed under the GNU LGPL or MIT". The vendored copy has been updated since 2016. It draws the sliders and 2D plots, and also provides the ODE integrator and the root finder. |
| `assets/js/elliptic.js` | Stein's own complete elliptic integrals K(m), E(m) and Π(n, m) by AGM iteration (Garrett 1963). Published separately as [elliptic-integrals-js](https://github.com/duetosymmetry/elliptic-integrals-js) (MIT). |
| MathJax 2.7.1 (CDN) | Equations in the write-up. |

It does **not** use MathBox, d3 or raw WebGL.

### 1.2 How the orbit is computed

**Parameters.** Units are M = 1 and E = 1. Teo's closed forms Φ(a, r) and Q(a, r) give the angular momentum and Carter constant from the spin and radius. The radius slider runs from r₁(a) = 2(1 + cos(⅔ acos(−a))) to r₂(a) = 2(1 + cos(⅔ acos(a))).

**Latitude variable.** The code uses u = cosθ. The latitude potential factorises as Σ²u̇² = a²(u₀² − u²)(u² − u₁²). Following Scott Hughes, Stein writes u = u₀ sin χ with a monotonic phase χ. This removes the sign switch at turning points and the Lipschitz problem.

**The ODE.** It is a 2D system for (χ, φ) in the **affine parameter** ζ (`makechidotphidotRHS`):

```js
// MIT-licensed site repo (see §1.8); abridged
var chidot = a/Sigma * Math.sqrt( U0sq*sinchisq - U1sq );
var phidot = (2.*r*a + (Sigma - 2.*r)*ell/(1. - U0sq*sinchisq)) / (Delta*Sigma);
```

It is integrated with `JXG.Math.Numerics.rungeKutta('rk4', …)`, a **fixed step** of dζ = 0.05:
- The history slider covers 15 to 300 units of ζ, so 300 to 6,000 points.
- Each animation tick adds 0.5 units of ζ (10 steps).
- r and t are not integrated: r is constant, and t "doesn't really add information".

**Frequency ratio.** The ratio Ω_φ/Ω_θ = Δφ/2π is computed in closed form from complete elliptic integrals:

Δφ = 4/√(−u₁²) · [(2r − aΦ)/Δ · K(u₀²/u₁²) + (Φ/a) · Π(u₀², u₀²/u₁²)]

A "Try to match" box root-finds r for a target ratio, using `JXG.Math.Numerics.fzero` on [r₁, r₂].

**Axis problem (question (a) in the brief).** **[V]** Yes, Stein says this, twice:
- Under "Code": "The numerical integrator from JSXGraph that I use is a fixed step-size RK4 integrator. It is sadly not adaptive! This is a problem for orbits close to the spin axis."
- Under "Limitations": "Orbits that pass very close to the spin axis experience rapid variation in φ in a small amount of affine parameter. These numerics will be especially bad. For exactly polar orbits … φ goes through a jump discontinuity of π when going across the pole. This isn't possible with the formulation I've got here."

He also warns that "the simulation is pretty wonky along the Φ = 0 line". The tool does not otherwise deal with the axis: the step size is simply "what seemed like a reasonable step size in most of parameter space".

The same failure shows up, measured, on a 5D orbit in §5.5. Note that the fix is not just "use an adaptive integrator":
- The rapid φ swing near the pole happens where the φ-circle is tiny. The *picture* barely changes, but the *accumulated* φ must be exact, or the rest of the orbit is rotated.
- Integrating the 1/sin²θ spike analytically (or adaptively), and sampling positions uniformly in χ, solves both problems.

### 1.3 Mapping to 3D, horizon, and other scene objects

**Mapping.** Boyer–Lindquist coordinates are treated as oblate spheroidal coordinates, as suggested by Steve Drasco (2016-03-05 update). three.js is Y-up, so the code swaps axes:

```text
(x, y, z) = (√(r²+a²) sinθ sinφ,  r cosθ,  √(r²+a²) sinθ cosφ)
```

**Horizon.** A `SphereGeometry(1, 32, 32)` scaled to (√(r₊²+a²), r₊, √(r₊²+a²)), using `MeshBasicMaterial` in black at opacity 0.5 (transparent, unlit). There is **no ergosphere**.

**Other objects:**
- a blue `ArrowHelper` for the spin vector;
- a translucent `GridHelper` at y = −6;
- an `AmbientLight`;
- a clear colour of `0xf0f0f0`, so a light theme only.

### 1.4 Camera and controls

- threestrap's default perspective camera is placed at (6, 5, 4).
- `THREE.OrbitControls` runs with `rotateSpeed = 0.5`, giving mouse and touch rotate, zoom and pan.
- The parameter sliders are **JSXGraph sliders drawn on a JSXGraph board**, not HTML `<input type="range">` elements.
- **[I]** As a result they are not native form controls, which limits keyboard and screen-reader access.

### 1.5 Responsiveness when sliders move

- **Every JSXGraph `drag` event** synchronously does all of the following:
  - recomputes Δφ;
  - re-integrates the **whole** history from χ = φ = 0;
  - rebuilds the horizon geometry;
  - rebuilds the trajectory geometry.
- There is no throttling, no worker and no precomputation.
- It stays smooth only because the work is tiny: a 2-variable RK4 over at most 6,000 steps, at well under a millisecond per input **[I]**.
- Side effect: the trajectory restarts from its initial point on every slider move.

### 1.6 How the curve is drawn

- A single `THREE.Line` with `LineBasicMaterial({color: 0xff0000, linewidth: 2})`.
- `linewidth` is ignored. The current three.js docs say "WebGL and WebGPU ignore this setting and always render line primitives with a width of one pixel" **[V]** (`src/materials/LineBasicMaterial.js`, r187). So the curve is a 1-pixel, aliased-ish, single-colour polyline with no depth cue.
- Each update creates a new `THREE.Geometry` from an array of `Vector3`s and disposes the old one, so it allocates on every tick.

### 1.7 Animation

- `setTimeout(runTick, 35)` (about 28 Hz, not tied to `requestAnimationFrame`) integrates 0.5 units of ζ forward.
- It `concat`s the new points onto the end of the solution array and `splice`s the same number off the front: a sliding window over an array of arrays.
- φ and χ are reduced mod 2π.
- There is no photon "head" marker and no fading.
- Play/Pause is a button that starts in the playing state.

### 1.8 Licence and reuse

- **The site repo:** GitHub reports MIT for https://github.com/duetosymmetry/web-site **[V]**. However, the `LICENSE` file is the upstream Minimal Mistakes theme's ("Copyright (c) 2016 Michael Rose"), inherited because the repo is a fork. **[I]** Whether that grant clearly covers Stein's own post scripts is ambiguous.
- **elliptic.js:** clearly MIT via elliptic-integrals-js **[V]**.
- **Recommendation [I]:** reimplement rather than copy. The algorithms are short, and equations are not copyrightable. If you adapt anything, credit Stein.

### 1.9 What it does well, and what we can improve

**Does well:**
- One idea, with very few controls (spin, radius, history).
- Immediate feedback.
- The χ phase trick.
- Physically meaningful mapping into oblate spheroidal coordinates.
- A frequency-ratio readout and resonance finder, with a table of presets for clickable resonances.
- An excellent long-form explanation directly under the tool, including an honest limitations section.
- Static 2D companion plots (u₀(r) and Δφ/2π against r).

**Could improve:**

| Area | Problem | Our fix |
| --- | --- | --- |
| Numerics | Fixed-step RK4 in affine time fails near the axis | §5.5 |
| Line | 1 px, single colour, no depth cues, no head marker, no fade | §5.2 to §5.4 |
| Slider moves | The orbit restarts on every move | Keep the phase |
| Animation | `setTimeout` instead of `requestAnimationFrame`; rebuilds geometry every tick and leaks memory | Preallocated buffers |
| Motion | No reduced-motion handling; autoplays | Start paused under `prefers-reduced-motion` |
| Theme and layout | Light theme only; fixed layout | Follow the site theme |
| Controls | Non-native sliders | Native range inputs |
| History units | Arbitrary affine parameter, so the trail's visual length varies with r | Measure history in latitude oscillations |

---

## 2. Other Kerr geodesic and photon-orbit visualisers

Gathered by a research sub-agent; I spot-checked some items.

### 2.1 Other duetosymmetry tools

There is no tools index: https://duetosymmetry.com/tools/ returns 404 and `/tool/` returns 403 **[V]**. All of these are JSXGraph pages **[V]**:

- **Kerr ISCO calculator** (2015): https://duetosymmetry.com/tool/kerr-isco-calculator/
- **Kerr Calculator V2** (2015): https://duetosymmetry.com/tool/kerr-calculator-v2/ , an interactive version of BPT72 Fig. 1. It uses JSXGraph `glider`s on closed-form curves, kept in sync from `board.on('update')`.
- **Poincaré section clicker** (2020, double pendulum): three.js, threestrap, OrbitControls and JSXGraph RK4. Its post repeats "JSXGraph only has fixed-step-size integrators, no adaptive ones".

Rivest & Stein presented "A Web Tool for Visualizing Kerr Spacetime and Bound Geodesics" at APS SES 2022: https://meetings-archive.aps.org/ses/2022/n02/4 **[V title]**. No live URL or repo was found.

### 2.2 Edward Teo

**Teo 2003**, GRG 35, 1909 **[V]**: https://physics.nus.edu.sg/wp-content/uploads/sites/5/2019/10/paper.pdf
- Constants of motion are closed form in r.
- Δφ per latitude oscillation is closed form via complete K and Π.
- Orbits are plotted by **fixed-step RK4 in affine parameter**, with the step chosen so that Δφ and u₀ match the exact values to 5 significant figures.
- Orbits are drawn on an imaginary sphere with numbered direction arrows.

**Teo's web page** **[V]:** https://phyweb.physics.nus.edu.sg/~phyteoe/kerr/ has static GIFs. The "wire-frame model" is a Java applet, which is dead technology.

**Teo 2021**, "Spherical orbits around a Kerr black hole", GRG 53, [arXiv:2007.04022](https://arxiv.org/abs/2007.04022) **[V]**:
- Gives **closed-form Mino-time solutions with Jacobi functions**, e.g. u = √w₁ sn(…λ, k).
- φ(λ) and t(λ) use incomplete Π, with the photon case in Appendix B.
- The figures are plotted from these closed forms.

No animations were found for either paper.

### 2.3 Yukterez (Simon Tyran)

**Mathematica code** **[V]:** https://yukterez.net/org/kerr.txt
- Evolves (t, r, θ, φ, p_r, p_θ) in proper or affine time.
- Because it evolves the momenta, it avoids square-root sign flips at turning points.
- Uses `NDSolve` with `InterpolationOrder→All` and a horizon event.
- Draws with `ParametricPlot3D` and a Hue-faded trail, rasterised frame by frame.

**Licences [V]:**
- The forum page says CC-BY 4: http://notizblock.yukterez.net/viewtopic.php?t=81
- The Wikimedia GIFs are CC BY-SA 4.0, e.g. "Kerr photon orbits with orbital inclination", which shows all spherical photon orbits for a = M: https://commons.wikimedia.org/wiki/File:Kerr_photon_orbits_with_orbital_inclination_thumbnail.gif

### 2.4 Black Hole Perturbation Toolkit and other analytic codes

**KerrGeodesics** (Mathematica, MIT) **[V]:** https://bhptoolkit.org/modules/KerrGeodesics/ , https://github.com/BlackHolePerturbationToolkit/KerrGeodesics
- `KerrGeoOrbit` has two methods: `"FastSpec"` (spectral, Hopper–Forseth–Osburn–Evans 2015) and `"Analytic"` (van de Meent & Warburton, arXiv:1802.05281).
- It covers timelike bound orbits.
- `KerrGeoPhotonSphereRadius` is closed form where possible.

**KerrGeoPy** (Python, MIT) **[V]:** https://github.com/BlackHolePerturbationToolkit/KerrGeoPy , https://kerrgeopy.readthedocs.io , [arXiv:2406.01413](https://arxiv.org/abs/2406.01413)
- Stable orbits use **Fujita & Hikida 2009** ([arXiv:0906.1420](https://arxiv.org/abs/0906.1420)): Jacobi elliptic functions in Mino time. `orbit.trajectory()` returns callables of Mino time.
- Plots are matplotlib 3D, with a linewidth "tail" that decays exponentially (`tau`).
- Timelike only.

**FEW** **[V]:** `few.utils.geodesic` (https://fastemriwaveforms.readthedocs.io/en/latest/tutorial/utility.html) gives only frequencies and constants, no positions.

**Null geodesics in closed form:**
- Gralla & Lupsasca 2020 ([arXiv:1910.12881](https://arxiv.org/abs/1910.12881)): every Kerr null geodesic in Jacobi form **[V]**.
- kgeo (Chael, **GPLv3**) implements it: https://github.com/achael/kgeo **[V]**. Learn from it; do not copy.
- geokerr (Dexter & Agol, [arXiv:0903.0620](https://arxiv.org/abs/0903.0620)) uses Carlson integrals **[V abstract]**.

### 2.5 Browser-based visualisers

**Ricky Reusser, "Geodesics in Kerr Spacetime"** (Mar 2026) **[V source]:** https://rreusser.github.io/notebooks/kerr-geodesics , repo https://github.com/rreusser/notebooks (MIT). This is the closest modern analogue and the most polished.
- **Integration:** Mino time, using the *second-order* form d²r/dλ² = R′(r)/2 and d²θ/dλ² = Θ′(θ)/2 (no sign tracking). Adaptive Cash–Karp RK4(5), tolerance about 1e-8, up to 8,000 steps into a `Float32Array`.
- **Rendering:** raw **WebGPU**, SDF-antialiased GPU lines, and 3-layer depth peeling for the transparent horizon and ergosphere.

**Alexander Nieddu, Kerr-Black-Holes** **[V]:** https://alexnieddu.github.io/Kerr-Black-Holes
- p5.js WEBGL, fixed-step RK4 in proper time, one step per frame.
- The trail is a growing array, redrawn in full every frame, so cost grows over time.
- **No licence.**

**bongo0/Kerr-orbits** **[V README]:** https://github.com/bongo0/Kerr-orbits , bounded photon orbits with p5.js and RK4 in Kerr–Schild Cartesian coordinates. **No licence.**

**Otto Seiskari, black-hole** **[V]:** https://github.com/oseiskar/black-hole , three.js and GLSL. Schwarzschild ray tracing only, using leapfrog on u″(φ). The code is MIT.

**Lensing ray tracers [S]:** Singularity (WebGPU), Gargantua (React Three Fiber) and adriwin06. These make images, not orbit plots.

**analyticphysics.com "Visualizing Aspects of the Kerr Metric"** **[V]:** https://analyticphysics.com/General%20Relativity/Visualizing%20Aspects%20of%20the%20Kerr%20Metric.htm , MathCell. It shows horizons and potentials, but no geodesics.

**Wolfram Demonstrations "3D Kerr Black Hole Orbits"** **[V via cloud notebook]:** https://demonstrations.wolfram.com/3DKerrBlackHoleOrbits/ . Needs Wolfram Player.

### 2.6 Desktop and offline tools

- **EinsteinPy** (MIT) **[V]:** https://docs.einsteinpy.com/en/stable/api/integrators/fantasy.html . Uses the **FANTASY** explicit symplectic integrator (Christian & Chan 2021, [arXiv:2010.02237](https://arxiv.org/abs/2010.02237)), not a Wang integrator. Plotting is matplotlib or Plotly.
- **Imaging codes [S]:** GRay2, Odyssey and Gyoto are GPU/C++ imaging codes, not orbit viewers.
- **Andrew Hamilton's black hole flight simulator** **[V]:** https://jila.colorado.edu/~ajsh/insidebh/bhfs.html . Not public, and "does not yet implement rotating black holes".

### 2.7 Lessons [I]

- **Mino time is the standard tool.** Teo 2021, Fujita–Hikida, KerrGeoPy and Reusser all use it: it decouples the latitude motion. At fixed r, the azimuthal angles are then integrals over θ.
- **Avoid square-root sign flips**, in one of these ways:
  - a phase variable (Stein, Hughes);
  - a second-order form (Reusser);
  - momenta (Yukterez);
  - closed-form `sn` (Teo, Gralla–Lupsasca).
- **Exploit periodicity.** Compute one latitude period and tile it. Validate Δφ against closed forms, as Teo does.
- **Polish benchmark:** Reusser's notebook for line quality and transparency. The only Kerr spherical-photon-orbit tools are Stein's and Teo's dead applet.

---

## 3. Myers–Perry and higher-dimensional black holes

Gathered by a research sub-agent, who read the paper's LaTeX source and PDF; I re-derived the key equations numerically (§5.5).

### 3.1 The paper itself (arXiv:1801.03389) [V]

**Metric (eq. 2.1).** Boyer–Lindquist-type coordinates with x = r², following Frolov–Stojković:

ds² = −dt² + (μ/ρ²)[dt + a sin²θ dφ + b cos²θ dψ]² + ρ²dx²/(4Δ) + ρ²dθ² + (x+a²) sin²θ dφ² + (x+b²) cos²θ dψ²

where:
- ρ² = x + a² cos²θ + b² sin²θ
- Δ = (x+a²)(x+b²) − μx
- x± = ½(μ − a² − b² ± √((μ−a²−b²)² − 4a²b²)), which requires μ ≥ (a+b)²

**[I]** The printed "(x+b)² cos²θ dψ²" is a typo for (x+b²). Papnoi et al. eq. 1 has (x+b²).

**Conserved quantities.** E (set to 1), Φ, Ψ and a Carter-like K. The examples use μ = 1.

**Equations of motion (eqs. 2.3).** With overdot = d/dλ:
- (ρ²θ̇)² = Θ, with Θ = E²(a² cos²θ + b² sin²θ) + K − Φ²/sin²θ − Ψ²/cos²θ
- ρ²φ̇ = Φ/sin²θ − μa(x+b²)ℰ/Δ − (a²−b²)Φ/(x+a²)
- ρ²ψ̇ = Ψ/cos²θ − μb(x+a²)ℰ/Δ + (a²−b²)Ψ/(x+b²)
- ℰ = E + aΦ/(x+a²) + bΨ/(x+b²)

**Spherical orbits.** The conditions 𝒳 = d𝒳/dx = 0 apply to a cubic in x whose coefficients are linear in K. The paper solves for (x, K) given (a, b, μ, Φ, Ψ) "in Mathematica or Maple".

**Integration.** The paper integrates numerically in Hamiltonian form (Liu–Ding–Jing) to avoid sign flips at θ turning points. It does not use Mino time or closed forms. It notes that Diemer et al. solved MP geodesics with Weierstrass functions.

**Visualisation (§4.1).** The paper names three methods:
1. A path in (θ, φ, ψ) space.
2. θ, φ, ψ plotted against λ (Figs 4.9 and 4.11).
3. Pseudo-Cartesian coordinates, collapsed to 3D:
   - x₁ = √(r²+a²) sinθ cos[φ − arctan(a/r)]
   - y₁ = √(r²+a²) sinθ sin[φ − arctan(a/r)]
   - z₁ = √(r²+b²) cosθ

   **ψ is shown as a colour wheel** (Fig. 4.1). The hue order is red at 0, yellow at π/3, green at 2π/3, cyan at π, blue at 4π/3, magenta at 5π/3: an HSV wheel.

Every figure shows a **black spheroid with a white equator ring**, with the orbit drawn on it. The orbit is red when ψ does not move, otherwise hue-coloured.

**Figures and presets.** Parameter sets, all with μ = 1:

| Fig | a | b | Φ | Ψ | Shows |
| --- | --- | --- | --- | --- | --- |
| 4.2 | 0.5 | 0 | 0 | 0 | All latitudes |
| 4.3 | 0.5 | 0 | 2 | 0 | Quasiperiodic equatorial band |
| 4.4, 4.5 | 0.8 | 0 | 0 | 0.5 | Polar orbit; rainbow "petals" |
| 4.6 | 0 | 0.99 | 0 | 0 | "Football" (prolate) hole |
| 4.7 | 0.6 | 0.3 | 1 | −1 | Band in [θ₋, θ₊] |
| 4.8 | 0.9 | 0.1 | 0.8 | 0.8 | Extremal; **cusps** (projection artefacts) and sharp colour changes |
| 4.10, 4.11 | 0.9 | 0.09 | 0.1 | 0.01 | Non-monotonic φ, ψ ("5D Lense–Thirring") |
| 4.12 | 0.5 | 0.25 | ≈ −0.9292 | ≈ 0.9036 | Circular orbit at θ = π/4 |

The sub-agent solved for (x, K) for these sets, and my independent JS solver in §5.5 reproduced the values for Figs 4.7, 4.8 and 4.10 **[M]**. They work as test vectors (Appendix B).

**Two subtleties for the demo [I]:**
- z₁ as written is a *norm*, so it is ≥ 0 for θ ∈ [0, π/2]. Yet in Fig. 4.2 the orbit dips below the equator, so the figure seems to use a *signed* cosθ (θ unfolded to [0, π]). Decide deliberately (§4.4).
- The curve labelled "x" in Figs 4.9 and 4.11 sits at about 1.38 and 1.13, which matches r = √x, not x.

**Earlier talk [V].** Bugden's 2015 Hobart talk (https://maths-people.anu.edu.au/~bugdenm/Hobart2015.pdf) already showed:
- orbits on a shaded sphere, including cusps;
- a 3D path in (θ, φ, ψ) space explaining where the cusps come from.

### 3.2 Papers citing arXiv:1801.03389 [V]

**Sources:**
- INSPIRE (`refersto:recid:1647371`): https://inspirehep.net/api/literature?q=refersto:recid:1647371
- Semantic Scholar: https://api.semanticscholar.org/graph/v1/paper/arXiv:1801.03389/citations
- OpenCitations: https://api.opencitations.net/index/v2/citations/doi:10.1007/s10714-018-2352-x

Together they give **9 genuine citing works**:

1. **Bugden 2020**, "Trapped photons in Schwarzschild–Tangherlini spacetimes", CQG 37 015001 ([arXiv:1909.07298](https://arxiv.org/abs/1909.07298)). Topology of the trapped set.
2. **Bugden**, PhD thesis on T-duality ([arXiv:1904.03583](https://arxiv.org/abs/1904.03583)). A self-citation.
3. **Habibina & Ramadhan 2020** ([arXiv:2007.03211](https://arxiv.org/abs/2007.03211)). 4D nonlinear-electrodynamics black holes; 2D plots only.
4. **Habibina, Jayawiguna & Ramadhan 2021** ([arXiv:2104.12071](https://arxiv.org/abs/2104.12071)). Equatorial 2D orbits.
5. **Cao, Liu & Wu 2022** ([arXiv:2206.09518](https://arxiv.org/abs/2206.09518)). 4D Kerr–Newman with string clouds and quintessence; 3D Cartesian orbit plots and projections. Cites Bugden as "a five-dimensional rotating black hole".
6. **Anjum, Afrin & Ghosh 2023** ([arXiv:2301.06373](https://arxiv.org/abs/2301.06373)). 4D Kerr with perfect-fluid dark matter; Teo-style 3D orbit plots.
7. **Novo, Cunha & Herdeiro 2024**, "Hypershadows…", JHEP 11 (2024) 171 ([arXiv:2410.05390](https://arxiv.org/abs/2410.05390)).
   - Says "Spherical photon orbits in the general MP BH were first discussed in [Bugden]".
   - For a = b, light rings exist at **every radius in an interval**.
   - Figures show the 3D "hypershadow" surface and slices, not trajectories.
8. **J. Yang 2026** ([arXiv:2511.21829](https://arxiv.org/abs/2511.21829)). Numerical hypershadows; point clouds and slices.
9. **E. Plácido Flores**, MSc thesis, UAM-Iztapalapa (2020), DOI 10.24275/uami.2b88qc47g. Metadata only.

**None of them re-plots 5D spherical photon orbits.**

### 3.3 Earlier visualisations of higher-dimensional geodesics [V]

| Work | Method | How orbits, or the 4th dimension, are shown |
| --- | --- | --- |
| Frolov & Stojković 2003, [gr-qc/0301016](https://arxiv.org/abs/gr-qc/0301016) | Separation in x = r²; constants E, Φ, Ψ, K | Effective-potential plot only |
| Gooding & Frolov 2008, [arXiv:0803.1031](https://arxiv.org/abs/0803.1031) | Capture surfaces | 3D surfaces in impact-parameter space, seen from four directions; 2D slices |
| Kagramanova & Reimers 2012, [arXiv:1208.3686](https://arxiv.org/abs/1208.3686) | a = b; Weierstrass functions in **Mino time** | Orbits in (X, Y, Z) with Z = √(u+a²) cosθ cosψ: a *linear projection dropping W*, with no colour. Cones bound the θ motion; spheres mark the horizons. |
| Diemer, Kunz, Lämmerzahl & Reimers 2014, [arXiv:1404.3865](https://arxiv.org/abs/1404.3865) | General a ≠ b; Weierstrass ℘, ζ, σ; Mino time | X, Y, Z, W coordinates with W dropped, in greyscale. Fig. 8 demonstrates **projection artefacts** explicitly. |
| Grunau, Neumann & Reimers 2018, [arXiv:1711.02933](https://arxiv.org/abs/1711.02933) | MP–AdS₅ | Photon region in the X–Z plane; orbit curves with wire-frame horizons |
| Delsate, Rocha & Santarelli 2015, [arXiv:1507.03602](https://arxiv.org/abs/1507.03602) | Equal angular momenta | **Hopf base CP¹ (θ, φ) plus fibre ψ**, with orbits in the (θ, φ) and (ψ, φ) planes. The base motion is a circle. |
| Papnoi et al. 2014, [arXiv:1407.0834](https://arxiv.org/abs/1407.0834); Amir, Singh & Ghosh, [arXiv:1707.09521](https://arxiv.org/abs/1707.09521) | a = b shadows | 2D celestial (α, β) shadows. **Watch the convention:** their 𝒦 equals K − Φ² − Ψ² in Bugden's notation **[I]**. |
| Hertog, Lemmens & Vercnocke 2019, [arXiv:1903.05125](https://arxiv.org/abs/1903.05125) | Ray tracing of 5D black holes and rings | 2D slices of 3D shadows |

**Summary of approaches.** Ways of showing the 4th dimension in the literature:
- colour (Bugden only);
- a linear projection that drops one coordinate;
- coordinate-space paths;
- coordinates plotted against λ;
- Hopf base and fibre;
- slices.

Nobody uses a (φ, ψ) torus, an interactive 4D rotation, or animation **[V for the papers above; I for "nobody"]**.

---

## 4. Showing a 4D curve in 3D, and colour as a dimension

Sources were read by a research sub-agent. I added my own measurements, in CAM02-UCS with colorspacious (Appendix A). The sub-agent's independent CIELAB computation gives the same ranking.

### 4.1 Cyclic map for ψ, or diverging map for w?

- **ψ is an angle, so it needs a cyclic map.** matplotlib: cyclic maps are for "phase angle, wind direction, or time of day"; "L\* should change monotonically from start to middle, and inversely from middle to end"; HSV is "a poor choice" (https://matplotlib.org/stable/users/explain/colors/colormaps.html) **[V]**. Crameri: "Use with circular data only!" (https://www.fabiocrameri.ch/cycliccolourmaps/) **[V]**.
- **w is signed and does not wrap, so it would need a diverging map** centred on 0, with symmetric limits ±r: "divergent about a central value (e.g., centred about zero)" (Crameri, Shephard & Heron 2020, *Nat. Commun.* 11, 5444, https://pmc.ncbi.nlm.nih.gov/articles/PMC7595127) **[V]**.
- **Kovesi's points** ([arXiv:1509.03700](https://arxiv.org/abs/1509.03700); https://colorcet.com/) **[V]**:
  - "the standard colour circle is not a good map".
  - Good cyclic maps pass through four distinct high-chroma colours "with some symmetry of lightness values".
  - Two-phase maps accept "a colour ambiguity corresponding to phase angles of 0 and 180 degrees", which is fatal for ψ.
  - Colour-blind-safe maps must lie in the 2D protan/deutan (lightness × blue–yellow) or tritan colour space.
  - **[I]** So a CVD-safe *cyclic* map must vary in lightness. Every constant-lightness hue wheel folds over for deuteranopes.
- **Lighting.** Keep the curve **unlit** (`LineMaterial` is unlit), so that lightness carries only data and never shading. Kovesi's argument for orthogonality between colour and shading **[V]** applies here.

### 4.2 Measured comparison of cyclic maps [M]

Each map was sampled at 256 points. Columns:
- **J′:** CAM02-UCS lightness range.
- **unif:** coefficient of variation of the perceptual step size; 0 is perfectly uniform.
- **<3:1:** share of the map below WCAG 3:1 non-text contrast against the kit's `--color-surface`: light #fbfaf6, dark #1a1b20. White and a navy #0f1b33 gave almost identical figures.
- **opp ΔE:** the minimum colour difference between ψ and ψ+π, for normal vision / deuteranopia / protanopia / tritanopia (Machado 2009, 100% severity). Below about 10, opposite phases start to look alike.

| Map (licence) | J′ | unif | <3:1 light | <3:1 dark | opp ΔE N / D / P / T |
| --- | --- | --- | --- | --- | --- |
| HSV hue wheel (as in the paper's figures) | 31–97 | 0.49 | 52% | 15% | 70 / 25 / 51 / 53 |
| matplotlib twilight | 16–89 | 0.13 | 37% | 43% | 42 / 33 / 32 / 41 |
| Crameri **romaO** (MIT) | 36–88 | 0.17 | 45% | 33% | 46 / **42 / 42 / 42** |
| Crameri **vikO** (MIT) | 22–82 | 0.19 | 34% | 48% | 46 / 40 / 36 / 50 |
| Crameri brocO / corkO / bamO | 23–86 | 0.15–0.20 | 37–41% | 36–43% | 32–33 / 23–32 / 28–33 / **7–22** |
| cmocean **phase** (MIT; isoluminant) | 57 | 0.00 | **0%** | **0%** | 55 / **2** / 16 / 17 |
| CET-C10 (isoluminant) | 68–71 | 0.18 | 100% | 0% | 35 / 3 / 6 / 9 |
| CET-C1, C2, C6, C7, C8, C9, C11 (CC-BY 4.0) | various | 0.20–0.41 | 26–100% | 0–25% | 39–50 / **4–25** / 6–19 / 11–44 |
| CET-C4, C5, CBC2, CBTC2 (two-cycle) | — | — | — | — | **≈1 for everyone: unusable for ψ** |
| **CET-CBC1** (blue–white–yellow–black; protan/deutan-safe) | 21–94 | 0.13 | 45% | 34% | 52 / **51 / 50** / 37 |
| **CET-CBTC1** (cyan–white–red–grey; built for tritanopes) | 45–98 | 0.16 | 66% | **4% (min 2.9:1)** | 40 / 29 / 21 / **45** |

The sub-agent's measure (the minimum CIE76 ΔE between colours ≥ 60° apart) gives the same ordering. For CBTC1 it found 20/19/12/19, for romaO 38/23/25/15, and 1 for deuteranopes on phase or an OKLCH circle **[V/C from the sub-agent]**.

### 4.3 No single map does everything; strategies [M + I]

**The constraint.** To clear 3:1 against *both* #fbfaf6 and #1a1b20, a colour's luminance must lie between about 0.13 and 0.29 (L\* about 43–60).

**The test.** I squeezed the CVD-safe maps into that band, rescaling CAM02-UCS lightness to 48–60:
- They pass 3:1 on both themes.
- But opposite-phase ΔE under CVD falls to about 11–19, from 40–50.
- Uniformity worsens (unif 0.29–0.42).
- The sub-agent got the same result from a hand-built blue/yellow loop at L\* 41–62 (sep60 ≈ 11).

So no single cyclic map is all of these: CVD-safe, ≥ 3:1 on both themes, and finely resolved. Pick one of these strategies:

**A. One CVD-safe map, plus a casing.**
- Draw the same `Line2` geometry about 2 px wider, in the theme's ink or contrast colour, *behind* the coloured line. This is standard cartographic line casing.
- The legend is identical in both themes.
- The near-white or near-black parts of the cycle read as a "hollow" line, which is acceptable.

**B. A dark orbit surface in both themes.** This is the paper's look: a black spheroid with a white equator ring.
- Inside the silhouette, the backdrop is then dark whatever the page theme, so a light-ish cyclic map such as **CET-CBTC1** (4% below 3:1 on dark) works there.
- Only the silhouette edge and the page around it need the casing.
- **[I] My preferred default.** It is faithful to the paper and needs only one map.

**C. A different map per theme.**
- The sub-agent proposes CBTC1 for dark plus a lightness-mirrored prototype for light.
- Stops for dark: `#26bbd7 #a3ddec #f9f8f8 #fac7bf #e99488 #a57972 #636465 #538e9c`.
- Stops for light: `#007995 #004958 #262626 #62332d #a55249 #aa746c #8b8b8b #3f8c9c`.
- The light prototype is **not yet perceptually equalised**, and the legend changes with the theme.

**D. A constant-lightness map** (cmocean phase, or an OKLCH circle).
- It passes 3:1 everywhere and is immune to shading.
- But it fails deuteranopes, so ψ must also be encoded another way (WCAG 1.4.1).

**Whichever strategy, also do these:**
- Show ψ redundantly in a linked **(φ, ψ) periodic square** (§4.4), which satisfies WCAG 1.4.1 Use of Color (https://www.w3.org/WAI/WCAG22/Understanding/use-of-color.html) **[V]**.
- Test with Chrome DevTools → Rendering → "Emulate vision deficiencies" (https://developer.chrome.com/docs/devtools/rendering/apply-effects) **[V]**.

**Credit and storage.** The CET maps are CC-BY 4.0 (https://colorcet.com/download/index.html) **[V]**, so credit Kovesi in "How this works". The Crameri maps are MIT (https://zenodo.org/records/1243862) **[V]**. A 256-entry table costs about 1–3 KB in the bundle.

**Aligning colour with the hidden coordinate [I, from the sub-agent].** CBTC1's lightness extremes sit near ψ = π/2 (white) and 3π/2 (dark grey), with cyan at 0 and red at π. So *lightness* tracks sign(w) = sign(sin ψ), and hue tracks cos ψ. This matters most in a linear projection that drops w (§4.4). There, ψ and −ψ coincide in 3D, and lightness, the channel everyone reads reliably, tells them apart. For other maps, rotate the table by a quarter turn. colorcet ships `_s25` shifted variants **[V]**.

### 4.4 Showing the fourth dimension: projections and companion views

**Precedents [V]:**
- **Niles Johnson's Hopf fibration video** (https://nilesjohnson.net/hopf.html; CC BY-NC 4.0):
  - computed in Sage and rendered with Tachyon;
  - "each point on the two-sphere has been given a distinct color and the circle which projects to it is drawn in the same color";
  - S³ drawn via a modified stereographic projection.
- **Dimensions** (Leys, Ghys, Alvarez), chapters 7–8: Hopf fibration, stereographic projection, Villarceau circles (https://www.dimensions-math.org/Dim_CH7_E.htm).
- **Hanson's Meshview:** "any arbitrary palette can be used to color code the 4D depth of each point", plus a "4D rolling ball" for full 4D rotation from mouse input (https://homes.luddy.indiana.edu/hansona/papers/meshview.pdf).
- **GL4D:** a survey of techniques: "wireframe representations, hyperplane slicing, color coding, view transformations, projection, and animation" (https://www.cse.cuhk.edu.hk/~cwfu/papers/GL4D/paper/GL4D.pdf).
- **Marc ten Bosch, 4D Toys:** 3D slices with a slider for the 4th axis (https://marctenbosch.com/ndphysics , https://4dtoys.com/).
- **MathBox** (Wittens, MIT; dormant since 2023): https://github.com/unconed/mathbox
- **A three.js fat-line Hopf fibration demo:** https://discourse.threejs.org/t/the-hopf-fibration-from-topology/19145

**Geometry of our curve [I].** With the flat-space embedding X = R_a sinθ cosφ, Y = R_a sinθ sinφ, Z = R_b cosθ cosψ, W = R_b cosθ sinψ (R_a = √(x+a²), R_b = √(x+b²)):
- A spherical orbit lies on a 3D ellipsoid, a squashed S³.
- Each θ = const level is a (Clifford-like) torus whose two angles are **exactly (φ, ψ)**.
- The orbit fills the shell between the θ₋ and θ₊ tori.

**Main-view options:**

| Option | What it shows | Pros | Cons |
| --- | --- | --- | --- |
| **Paper's view:** (X, Y, \|Z+iW\|) = (R_a sinθ cosφ, R_a sinθ sinφ, R_b cosθ), colour = ψ | A dome, or a spheroid if θ is unfolded | Faithful to arXiv:1801.03389; ψ fully encoded by colour; the constant-r spheroid is a natural backdrop | Uses the hemisphere θ ∈ [0, π/2] only, unless cosθ is signed or unfolded as Fig. 4.2 appears to do. Cusps are projection artefacts (Fig. 4.8). |
| Linear projection dropping W (Kagramanova–Reimers; DKLR) | (X, Y, R_b cosθ cosψ), colour = sign(W) or ψ | Familiar 3D look; a 4D rotation slider in the (Z, W) plane gives motion parallax | ψ and −ψ overlap, and the picture changes with the rotation |
| Stereographic projection of the S³ from a point the orbit never reaches (e.g. θ = 0 when Φ ≠ 0) | The whole orbit, with no information lost | Matches the Hopf and Dimensions pictures; the θ-shell becomes nested tori of revolution; frees colour for time | Distorts sizes; abstract for non-specialists |

**Companion views, in priority order [I]:**
1. **(φ, ψ) periodic square (flat torus)** with a cursor tracking the photon:
   - Shows winding numbers and quasi-periodicity at a glance.
   - Doubles as the non-colour encoding of ψ.
   - Nobody in the literature uses it.
2. **θ, φ, ψ against time** (the paper's Figs 4.9 and 4.11), over the visible window.
3. Optionally, a (Z, W) rotation slider or a stereographic toggle, for the curious.

### 4.5 Legend [I, with V for the CSS details]

- **A ring, not a bar.** Draw it with CSS `conic-gradient`, which starts at 12 o'clock and runs clockwise (https://developer.mozilla.org/en-US/docs/Web/CSS/gradient/conic-gradient) **[V]**.
  - For the maths convention (0 at 3 o'clock, anticlockwise), use `from 90deg` with the stops reversed.
  - Use 32–64 explicit stops sampled from the same table, rather than relying on the browser's interpolation.
- **Labels:** 0, π/2, π and 3π/2 outside the ring.
- **A live needle or dot** at the photon's current ψ.
- **A 1 px border**, so the white and black sectors stay visible in both themes.
- **Accessibility:** `role="img"` and an `aria-label`.
- **Placement:** next to the stage, not inside the WebGL canvas, so it stays crisp and themeable.

---

## 5. Implementation tooling and performance

### 5.1 Budget and measured bundle sizes [M]

**The budget [V].** The demo-kit contract allows a single `demo.js` and `demo.css`, no CDN, and a bundle under 2 MB (1 KB = 1000 B). The validator sums raw on-disk file sizes (`bin/validate.js`), so **minified, not gzipped, size is what counts**. The current stub is 7.2 KB JS + 10.2 KB CSS + about 68 KB of images.

**The measurements.** I built each scenario with the repo's own Vite (8.3.8, Rolldown with the Oxc minifier), with code splitting off as in `demo/vite.config.js`. Each entry used only the named imports, so tree-shaking applied.

| Scenario | Minified | gzip −9 |
| --- | ---: | ---: |
| three r187: WebGLRenderer, Scene, PerspectiveCamera, BufferGeometry, Line, Mesh, SphereGeometry, MeshBasicMaterial | 529 KB | 129 KB |
| … + `OrbitControls` | 552 KB | 134 KB |
| … + `Line2`, `LineMaterial`, `LineGeometry` | **580 KB** | **141 KB** |
| … + MeshStandard and Phong materials, lights, fog, Points | 589 KB | 143 KB |
| three + OrbitControls + `TubeGeometry`/`CatmullRomCurve3` (instead of Line2) | 570 KB | 139 KB |
| three + OrbitControls + MeshLine 3.3.1 (instead of Line2) | 564 KB | 137 KB |
| `import * as THREE` + the same addons (no tree-shaking) | 781 KB | 194 KB |
| Observable Plot 0.6 (`plot`, `lineY`, `ruleY`) alone | 275 KB | 90 KB |
| three + OrbitControls + Line2 + Plot | 855 KB | 231 KB |
| OGL 1.0.11 (Renderer, Camera, Orbit, Polyline, Sphere, Program, Mesh) | 60 KB | 17 KB |
| regl 2.1.1 (npm `main`, debug checks included) | 121 KB | 39 KB |
| @stdlib `ellipj` + `ellipk` + `ellipe` | 43 KB | 16 KB |
| `dopri` 0.0.13 (Dormand–Prince) | 14 KB | 4.5 KB |
| `odex` 3.0.0-rc.4 (Gragg–Bulirsch–Stoer) | 14 KB | 4.5 KB |

**What the numbers say:**
- `WebGLRenderer` pulls in nearly all of three's shader library, so **tree-shaking buys only about 25%**. Lit materials cost almost nothing extra once the renderer is in.
- three + controls + fat lines fit comfortably, with about 1.3 MB of headroom.
- Plot is affordable but not free. It is also SVG-based, so a poor fit for a panel redrawn every frame (§5.9).
- **OGL** is the lightweight three.js-like alternative, 10× smaller. It has its own `Polyline` (screen-space thick lines) and `Orbit` control. But it has a much smaller community and fewer guarantees **[I]**.

### 5.2 Ways to draw a smooth, thick, per-vertex-coloured 3D curve

| Option | Thick? | Per-vertex colour | Anti-aliasing | Updating | Notes |
| --- | --- | --- | --- | --- | --- |
| three `Line` + `LineBasicMaterial` | **No**, always 1 px **[V]** | Yes (`vertexColors`) | MSAA only | `setDrawRange(start, count)` on a preallocated `BufferAttribute`: true sliding window, zero copies | Fine as a fallback or for a faint "full orbit" ghost |
| three **`Line2` / `LineMaterial` / `LineGeometry`** (addons) | Yes: screen-space, width in CSS px, or `worldUnits` | **Yes, RGB only**: `instanceColorStart/End` are `vec3`; opacity is a single uniform **[V]** | MSAA + `alphaToCoverage` **[V]** | Instanced segments; `geometry.instanceCount` limits the count drawn **[V]**, but there is no "start" offset | Official example (`webgl_lines_fat`) uses per-vertex HSL colours. Fog supported **[V]**. Recommended. |
| three `TubeGeometry` | Yes, real 3D tube; shades with lights | Yes, via a colour attribute | MSAA | Rebuilds geometry (radial × tubular vertices) on every change; costly for 20k+ points | Looks great as a static "hero" path with lighting. Width is in world units, so it thins with distance (a natural depth cue). |
| MeshLine 3.3.1 (pmndrs; MIT) | Yes | **No per-vertex colour**: only a 2–3-stop `gradient`, a texture `map`, or `alphaMap` **[V]** | Own shader | `setPoints()` rebuilds JS arrays; `advance()` shifts | Has a per-point width callback (tapering tails). Weaker fit for colour-by-ψ. |
| OGL `Polyline` | Yes | Via a custom shader | Own | Rebuild or update the attribute | Smallest bundle |
| regl or raw WebGL2 with a custom instanced line shader | Yes | Anything, e.g. the colour map looked up from ψ in the shader | SDF in the fragment shader | Full control | Most work; this is what Reusser does in WebGPU |
| 2D canvas + hand-rolled projection | Yes (`lineWidth`) | Per segment (one stroke per segment, or gradients) | Canvas AA | Re-stroke every frame; 10k+ segments per frame is slow | No depth buffer, so occlusion by the horizon means painter's sorting by hand. Not recommended for 3D. |

**[I] Recommendation:** three.js `Line2`/`LineMaterial` for the trail, a plain `Line` for an optional faint full-orbit ghost, and `Mesh`es for the horizon and orbit surfaces.

If you later want colour-by-ψ evaluated per *fragment* (exact hue sweeps within a segment), patch `LineMaterial` with `onBeforeCompile`, or replace it with a small custom instanced line shader. Feed it unwrapped ψ as a per-vertex scalar and a 256-entry colour-map texture.

### 5.3 `Line2`/`LineMaterial` details and pitfalls (three r187) [V unless marked]

**Resolution.** `LineSegments2.onBeforeRender` sets `material.resolution` from `renderer.getViewport()` (CSS pixels) on every render. So `linewidth` is in **CSS pixels** and devicePixelRatio is handled automatically.
- Older tutorials say to update `resolution` on resize by hand; that is no longer needed.
- If you render into a sub-viewport, such as an inset, the resolution follows the viewport.

**`setPositions` and `setColors` allocate.**
- Each call creates a new `InstancedInterleavedBuffer` and recomputes the bounding box and sphere, so **don't call them per frame**.
- Allocate once at the maximum size with `setUsage(DynamicDrawUsage)`. Write into `geometry.attributes.instanceStart.data.array` (interleaved start xyz and end xyz, 6 floats per segment), set `needsUpdate` (or `addUpdateRange`) and adjust `geometry.instanceCount`.
- Then either recompute the bounding sphere or set `line.frustumCulled = false`. **[I]** Otherwise a stale bounding sphere can cull the line wrongly.

**The format doubles the data.** `LineGeometry` turns N points into N − 1 segments with duplicated endpoints (8 vertices and 18 indices per segment instance). 20k segments is still trivial for any GPU **[I]**.

**Colour has no per-vertex alpha.**
- A fading tail must blend the colour toward the background colour, recomputed on theme change; or use a shader patch.
- With `transparent: true` and opacity < 1, the overlapping round end caps at each joint blend twice and show "beads" **[I, from the shader structure]**. Keep the trail opaque and fade by colour.

**Anti-aliasing.** Create the renderer with `antialias: true` and set `material.alphaToCoverage = true`. The official example does this **[V]**.

**Colour space.**
- three treats vertex colours as linear-sRGB.
- Convert colour-map entries from sRGB, as the official example does with `color.setHSL(t, 1, 0.5, THREE.SRGBColorSpace)`. Otherwise the map's perceptual uniformity is lost: it shows up washed out or too dark **[I]**.
- Precompute a 256-entry linear-RGB table once.

**Sliding windows.** `instanceCount` always draws from instance 0, and WebGL2 has no portable base-instance offset **[I]**. For a moving window, either:
- rewrite the visible window into the buffer each frame (cheap: 20k segments × 24 B is about 0.5 MB); or
- draw everything and hide out-of-window segments in the shader, using a per-segment parameter attribute and `uHead`/`uHistory` uniforms.

### 5.4 Updating the curve and drawing the history trail [I, with measured costs]

**On a parameter change** (a, b, Φ, Ψ):
1. Solve for (x, K).
2. Compute one latitude-period table at N_draw nodes (256–512) in χ, plus Δφ, Δψ and the period in Mino time (or affine or coordinate time).
3. Keep the *phase* (the current χ, or the number of periods elapsed) so that the photon does not jump back to the start.

**Each frame:**
- Advance the head parameter.
- Tile periods, rotating φ by nΔφ and ψ by nΔψ, to fill the visible window into the preallocated buffer.
- Set `instanceCount`.

Tiling 20,480 points into a `Float32Array` takes **1.2 ms in Node 24 (V8)** **[M]**.

**History slider.** Measure history in **latitude oscillations** (θ-periods) rather than raw affine parameter. It reads the same at every r and maps directly onto the tiling. Show the equivalent Δλ or Δt in the status line if you like.

**Point budget.** For under 0.25 px chord error on a sphere about 250 px across, you need about 70+ segments per great-circle turn **[I]** (sagitta s²/8R). Cyclic colour sweeps want ≥ 64 vertices per ψ turn. So 256–512 per θ-period is plenty, and 100 periods is still only 25k–50k segments.

**Adaptive refinement.** Near the axes, φ (or ψ) swings by ~π over a tiny χ interval. The *position* barely moves, but the *colour* (ψ) can jump across the whole wheel within one segment. Either:
- subdivide where consecutive |Δψ| exceeds ~10°, which is cheap because ψ(χ) is available in closed form plus quadrature; or
- interpolate ψ (not RGB) in the shader.

**The photon head:**
- a small sphere or sprite at the interpolated head position;
- a short bright section at the front, fading toward the background colour at the back;
- under reduced motion, no auto-play, but scrubbing the head still works.

**Coalescing input.** Coalesce slider `input` events into one recompute per animation frame, using demo-kit's single rAF loop, which calls `render` after control changes. The measured compute cost (§5.5) is far below a frame, so no worker, throttling or debouncing is needed.

### 5.5 Computing the orbit: methods compared [M]

**Equations used.** These are Bugden 2018 eqs. 2.3, as transcribed by the sub-agent, with E = 1, Mino time dγ = dλ/ρ², and u = cos²θ. The φ and ψ lines are verified; the substitution is my derivation, checked numerically:

```text
(du/dγ)² = 4P(u),   P(u) = u(1−u)(b² + (a²−b²)u + K) − Φ²u − Ψ²(1−u)
dφ/dγ = Φ/(1−u) + C_φ,   C_φ = −μa(x+b²)ℰ/Δ − (a²−b²)Φ/(x+a²)
dψ/dγ = Ψ/u     + C_ψ,   C_ψ = −μb(x+a²)ℰ/Δ + (a²−b²)Ψ/(x+b²)
dλ/dγ = ρ² = x + b² + (a²−b²)u
```

**Why this works:**
- P is a **cubic** in u when a ≠ b, which is Jacobi/Weierstrass territory (Kagramanova–Reimers; DKLR).
- When a = b it is quadratic, so u(γ) is a pure sinusoid and everything is elementary.

**Substitution [I, checked numerically].** Write u = u₁ + (u₂ − u₁) sin²χ, where u₁ and u₂ are the roots bounding the allowed band. Then P = (u − u₁)(u₂ − u)L(u) with L **linear**, and

dχ/dγ = √L(u(χ)) > 0

This is a 5D analogue of Stein's χ: no turning points, no sign flips. Switching the independent variable to χ turns φ, ψ, λ and t into **plain integrals of explicit functions of χ**. One θ-period is χ ∈ [0, π], and φ(χ + nπ) = φ(χ) + nΔφ.

**Handling the axis spikes.** As u₂ → 1 (Φ → 0) or u₁ → 0 (Ψ → 0), the terms Φ/(1−u) and Ψ/u spike. Subtract the singular part, which integrates exactly:

∫ dχ / (α − β sin²χ) = atan2(√(α−β) sin χ, √α cos χ) / √(α(α−β))

The remainder Φ(1/√L(u) − 1/√L(u₂))/(1−u) is bounded, so plain Gauss quadrature on a uniform χ grid converges fast.

**Results** (Node 24; errors in Δφ and Δψ per θ-period against a 2×10⁶-node reference). Fig. 4.10 has θ ∈ [0.049, 1.565], so it passes near both axes:

| Method | Fig 4.7 (generic) | Fig 4.8 (extremal) | **Fig 4.10 (near both axes)** |
| --- | --- | --- | --- |
| Fixed RK4, affine λ, Δλ = 0.05 (Stein's setting; about 55–66 steps/period) | ~1e-9 | ~3e-8 | **Δφ 3.6e-3, Δψ 1.8 rad** |
| Fixed RK4, Mino time, 64 steps/period | 3e-10 | 7e-9 | Δψ 0.56 rad |
| Fixed RK4, Mino time, 256 steps/period | 1e-12 | 3e-11 | Δψ 0.09 rad |
| Fixed RK4, Mino time, 1024 steps/period | 1e-13 | 6e-13 | Δψ 1.6e-5 |
| **Singularity-subtracted 3-pt Gauss, 64 panels** | 3e-14 | 5e-14 | **Δφ 6e-8, Δψ 4e-9** |
| **Singularity-subtracted 3-pt Gauss, 256 panels** | 3e-14 | 5e-14 | **1e-13 / 1e-8** |

**Timings:**
- One θ-period by subtracted quadrature (128 panels): **about 9 µs**.
- An unoptimised (x, K) solver (20,000-point scan plus bisection): about 0.3 ms. Use Newton continuation from the previous root while dragging.
- Hand-written DOPRI5 on a structurally identical test system with Φ = 0.002: **168 steps and 0.07 ms per period** at tolerance 1e-9, and it holds the tolerance. So adaptive RK45 *also* works, provided it integrates one period and tiles, and uses dense output for drawing (DOPRI5's step points are too sparse to draw directly).
- Fixed-step quadrature or RK4 without the subtraction needs more than 4,000 points per period at Φ = 0.002.

**Recommendation [I]:** at fixed r, **don't integrate an ODE at all**:
1. Use closed-form u(χ).
2. Compute φ, ψ, λ (and t) by singularity-subtracted Gauss quadrature over one period.
3. Tile.

This is simpler than RK45 in a worker. It is exact over long histories (no drift) and handles the axes. The frequency ratios Δφ/2π and Δψ/2π come for free, giving 5D analogues of Stein's Ω_φ/Ω_θ readout and of the resonance finder.

**Other notes:**
- **Exact closed forms.** If you want closed forms instead of quadrature, θ(γ) is Jacobi `sn²` and φ, ψ are incomplete Π. You would need Carlson R_F/R_J (DLMF §19.36, https://dlmf.nist.gov/19.36 ; Carlson 1995, [arXiv:math/9409227](https://arxiv.org/abs/math/9409227)), about 80 lines. The quadrature route avoids this.
- **demo-kit's `rk4`** allocates a new `Float64Array` per step and copies each slope **[V]**. It is fine for small jobs, but avoid it in a per-frame path.
- **Symplectic integrators** (FANTASY, Wang) matter for long-time energy conservation of *generic* orbits. They are irrelevant here, because the orbit is reduced to quadrature **[I]**.
- **Checks.** Test the extremal case, where the horizon's square root is exactly zero, and clamp `Math.sqrt(Math.max(0, …))`. My first run of Fig 4.8 got NaN without the clamp **[M]**.

### 5.6 Libraries for elliptic functions and ODEs [V for licence and version; M for size]

| Library | Provides | Licence | Bundled |
| --- | --- | --- | --- |
| [elliptic-integrals-js](https://github.com/duetosymmetry/elliptic-integrals-js) (Stein) | Complete K, E, Π (AGM) | MIT | ~1 KB |
| [@stdlib/math-base-special-ellipj](https://github.com/stdlib-js/math-base-special-ellipj), `-ellipk`, `-ellipe` | Jacobi sn, cn, dn, am; complete K, E | Apache-2.0 | 43 KB min together |
| [dopri](https://github.com/mrc-ide/dopri-js) 0.0.13 | Dormand–Prince 5(4) with dense output | MIT | 14 KB |
| [odex](https://github.com/littleredcomputer/odex-js) 3.0.0-rc.4 | Gragg–Bulirsch–Stoer with dense output | BSD-2-Clause | 14 KB |
| [ode45-cash-karp](https://github.com/scijs/ode45-cash-karp) (scijs) | Adaptive Cash–Karp | MIT | small; last published 2022 |

**More findings from the sub-agent [V]:**
- **stdlib** has only `ellipj`, `ellipk` and `ellipe`: no incomplete F, E or Π, and no Carlson functions. Its npm packages bundle fine with Vite **[M]**. Do not use its jsDelivr ESM build, which imports from a CDN.
- **[paulmasson/math](https://github.com/paulmasson/math)** (MIT) has Carlson R_C/R_D/R_F/R_G/R_J, incomplete F, E and Π, sn, cn, dn, am and Weierstrass functions, with complex arguments. But it is a single global `build/math.js` of 145 KB, not on npm and not tree-shakeable.
- **math.js** has no elliptic functions. **numeric.js** has none.
- **npm `elliptic` is elliptic-curve cryptography**, a trap.

**[I] Recommendation:** none are needed for the quadrature route. If you go closed-form, hand-write Jacobi `sn`/`cn`/`dn` by descending Landen/AGM (DLMF §22.20, https://dlmf.nist.gov/22.20 ; about 40 lines) and Carlson R_F/R_J. A tiny tested module beats a dependency here.

### 5.7 Scene: horizon, orbit surface, depth cues and controls

**Surfaces [I].** Following the paper's look:
- an **orbit spheroid** (the constant-r surface; oblate if a² > b², prolate otherwise) in a dark, matte, slightly translucent shade;
- a white equator ring;
- the **horizon** spheroid inside it.

Draw the opaque trail first. Then draw the surfaces with `transparent: true` and `depthWrite: false`, front faces only, with an explicit `renderOrder`: the inner horizon, then the outer orbit surface. Both are centred at the origin, so three's per-object distance sort cannot order them.

**Avoid z-fighting [I].** The trail lies *on* the orbit surface. Shrink that surface by about 0.5%, or use `polygonOffset`.

**Depth cues:**
- the segments behind the translucent surface are dimmed by it;
- `Scene.fog` (supported by `LineMaterial` **[V]**);
- optionally a second pass of the trail with `depthFunc: GreaterDepth` at low opacity, so the hidden half reads as "behind" rather than vanishing;
- a cheap Fresnel rim on the spheroid makes its 3D shape readable.

**Controls [V]:**
- `OrbitControls` supports damping, touch (one finger rotates, two fingers dolly or pan), and keyboard via `listenToKeyEvents`. The arrow keys **pan**; **Ctrl, Meta or Shift + arrow** rotates.
- It sets `touch-action: none` on its element. **[I]** On a phone the canvas then swallows page scrolling, and the wheel zooms instead of scrolling the article. Options: require a click or focus before wheel-zoom, `zoomToCursor` off, a modest stage height (the manifest says 560 px), and a "Reset view" button.
- **[I]** With `enableDamping`, the controls need frames to settle while the demo is paused. Re-render on the controls' `change` event until it stops, or disable damping when paused.

**Pixel ratio [I].** `renderer.setPixelRatio(Math.min(devicePixelRatio, 2))` caps fill cost on phones. Line widths stay in CSS px.

### 5.8 Accessibility

See §4 for colour-blind safety. Further points:

- **Reduced motion [V].** demo-kit's frame already reads `prefers-reduced-motion` and starts paused, and its Play/Pause satisfies WCAG 2.2.2 (Pause, Stop, Hide).
  - **[I]** Also disable auto-rotate and the camera damping "flick", and let the head be scrubbed by a slider.
- **WCAG 2.2.2 (Pause, Stop, Hide)** requires a pause control for automatic motion lasting more than 5 s **[V]**: https://www.w3.org/WAI/WCAG22/Understanding/pause-stop-hide.html
- **WCAG 2.3.3 (Animation from Interactions):** respecting `prefers-reduced-motion` is a sufficient technique **[V]**: https://www.w3.org/WAI/WCAG22/Understanding/animation-from-interactions.html
- **WCAG 2.5.7 (Dragging Movements, AA):** drag-to-rotate needs a non-drag alternative **[V]**: https://www.w3.org/WAI/WCAG22/Understanding/dragging-movements.html
  - **[I]** Add rotate and zoom buttons or a "Reset view" button.
  - OrbitControls' keyboard rotation (Ctrl/Shift + arrows) is undiscoverable, so document it in "How this works".
- **WCAG 1.4.11 (Non-text Contrast)** asks 3:1 for the line against its background **[V]**: https://www.w3.org/WAI/WCAG22/Understanding/non-text-contrast.html . See §4.3.
- **Native controls [I].** Use kit sliders (native inputs): keyboard operable, labelled, with formatted values. Avoid canvas-drawn sliders like Stein's JSXGraph ones.
  - Native `<input type=range>` has the implicit slider role **[V]**.
  - Add `aria-valuetext` for symbols, e.g. "Φ = 0.80".
  - The pattern is described at https://www.w3.org/WAI/ARIA/apg/patterns/slider/ **[V]**.
- **Text alternative [I].** Give the stage `role="img"` and a live `aria-label`, or use the status line, to describe the current orbit. For example: "Orbit at r = 1.33, θ between 37.6° and 51.5°, Δφ/2π = 0.336, Δψ/2π = −0.608."
- **Status line [I].** Use it for the "no spherical orbit for these momenta" state rather than silently clamping.

### 5.9 Fit with demo-kit 0.1.0 [V]

- **Stage.** Use `stage: 'div'` and append the three.js canvas. `render({ width, height, pixelRatio, reason, … })` runs on control changes, resizes, theme changes and every animation frame while playing; there is one rAF loop and it stops off-screen.
  - Call `renderer.setSize(width, height, false)` on `reason === 'resize'`.
  - Re-read tokens on `reason === 'theme'`, for the background, ink and casing colours.
- **Manifest.** Set `"needsWebGL": true`; it is currently `false`.
- **Workers.** `runInWorker` spawns a fresh Worker from a Blob URL **per call** and serialises the function **[V]**. That suits one-off long jobs, not per-input recomputes. It isn't needed here (§5.5).
- **Side panel.** `plot()` (Observable Plot, +275 KB) re-renders SVG. For an animated θ, φ, ψ against λ panel, either draw the curves with Plot once per parameter change and move a cursor overlay, or hand-draw on a small `<canvas>` (cheapest, and theme tokens are easy).
- **Fonts.** The kit's tokens name Fraunces and IBM Plex but don't ship them **[V]**, so a standalone demo uses the fallback stacks. Third-party fonts are not allowed in the bundle anyway.

---

## 6. Recommended architecture [I]

1. **Physics module** (pure functions, unit-tested against Appendix B):
   - `solve(a, b, Φ, Ψ, μ = 1) → {x, K} | null`, using a closed-form K(x) and a root of 𝒳 with Newton continuation;
   - `band(x, K) → {u₁, u₂, L}`;
   - `period(N) → tables of (χ, u, φ, ψ, λ, t)` plus Δφ, Δψ, Λ, T, by singularity-subtracted Gauss quadrature;
   - `fill(window, out: Float32Array)` by tiling.
2. **Rendering:**
   - three.js r187 (tree-shaken), `OrbitControls`;
   - one preallocated `Line2` (max ~50k segments, `DynamicDrawUsage`) with `alphaToCoverage`;
   - a head marker;
   - translucent orbit and horizon spheroids, plus the equator ring;
   - optionally a faint 1-px `Line` "ghost" of the long-term orbit;
   - colours from a 256-entry linear-RGB LUT of a cyclic CVD-friendly map, plus a ring legend;
   - a casing pass: the same geometry, 2 px wider, in the theme ink colour, drawn first.
3. **Update policy:**
   - Recompute synchronously in `render` when `reason === 'controls'` (under 2 ms).
   - Advance and refill the window when `reason === 'tick'`.
   - Never allocate per frame.
   - Keep the phase across parameter changes.
4. **Side panel:** a small 2D canvas with two parts:
   - a (φ, ψ) periodic square with the trail and a head dot, which encodes ψ without colour and shows winding;
   - θ, φ, ψ over the visible window, with a moving cursor.

   Plus Δφ/2π and Δψ/2π readouts.
5. **Presets:** the paper's Figs 4.2–4.12, including the cusp, Lense–Thirring and constant-θ cases. Add an a = b "every radius is a light ring" mode (Novo et al. 2024).

## 7. Decisions to make

- **What is z?** One of:
  - the paper's norm z₁ = √(r²+b²) cosθ ≥ 0, giving a dome with colour = ψ (faithful to the formula; the orbit "reflects" at the equator when Ψ = 0);
  - a signed, unfolded cosθ (matches Fig 4.2);
  - a linear projection z = √(r²+b²) cosθ cosψ with colour = sinψ or w (DKLR style; cusps appear and vanish as you rotate in 4D).

  A 4D-rotation slider in the (z, w) plane would turn this choice into a teaching moment.
- **Which parameter animates the head and measures history?** Options are θ-periods, Mino time, affine λ, or BL coordinate time t (what a distant observer would time).
- **Inputs.** (Φ, Ψ) can give zero or several orbits. Consider a 2D (Φ, Ψ) pad with the valid region shaded, instead of two independent sliders.

---

## Appendix A: methods and versions [M]

- **Bundle sizes:** a scratch project outside the repo, with three 0.187.0, vite 8.3.8, ogl 1.0.11, regl 2.1.1, meshline 3.3.1, @observablehq/plot 0.6.x, @stdlib/math-base-special-ellipj 0.4.1, -ellipk 0.2.4, dopri 0.0.13 and odex 3.0.0-rc.4.
  - Each scenario was built as a single-entry build with `codeSplitting: false`.
  - Sizes are `wc -c` of `demo.js` and `gzip -9c | wc -c`.
- **Numerics:** Node 24.20.0 (V8), plain JS, averaged over 100–1000 repetitions after a warm-up.
  - The reference used 10-point Gauss–Legendre on 200,000 panels.
  - DOPRI5 was hand-written with the standard Dormand–Prince tableau.
- **Colour metrics:** Python 3.13 with colorspacious (CAM02-UCS; CVD simulation at 100% severity, Machado et al. 2009), colorcet 3.2.1, cmocean 3.0.3, cmcrameri and matplotlib 3.11.2.
  - Each map was sampled at 256 points.
  - Contrast is the WCAG 2 ratio.
- **Page and source reads:** the Stein page HTML and every script it loads (three r74, threestrap, OrbitControls, JSXGraph 1.12.2, elliptic.js) were downloaded and read on 2026-11-20.
- **Fat-line behaviour:** read from three r187 `examples/jsm/lines/*.js`, `src/materials/LineBasicMaterial.js` and `examples/jsm/controls/OrbitControls.js`.

## Appendix B: test vectors [M; (x, K) also V from the sub-agent's independent solve]

All with μ = 1 and E = 1. Δφ and Δψ are per θ-period. They are my computation and have not been checked against the paper, which does not list them.

| Case | a | b | Φ | Ψ | x (r) | K | θ band | Δφ/2π | Δψ/2π |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Fig 4.7 | 0.6 | 0.3 | 1 | −1 | 1.7546 (1.3246) | 4.0204 | [0.656, 0.899] | 0.336346 | −0.607973 |
| Fig 4.8 | 0.9 | 0.1 | 0.8 | 0.8 | 1.8990 (1.3781) | 4.5222 | [0.388, 1.154] | 0.306357 | 0.535533 |
| Fig 4.10 | 0.9 | 0.09 | 0.1 | 0.01 | 1.2674 (1.1258) | 3.3606 | [0.049, 1.565] | 0.268751 | 0.464868 |

The sub-agent also reports these (x, K) pairs: Fig 4.2 (1.75, 3.75), Fig 4.3 (2.75, 5.75), Fig 4.4 (1.3168, 3.2404) and Fig 4.12 (1.4553, 3.2040, θ ≈ π/4).
