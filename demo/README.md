# Demo

The interactive demo for this project. It is built with [demo-kit](https://github.com/Mark-Bugden/demo-kit), which is pinned in `package.json`.

## Run it locally

```sh
npm install
npm run dev
npm test          # node:test: the physics, against the reference model's fixtures, the drawn path, the pad's geometry and the panel's ticks
```

## Build and check the bundle

```sh
npm run build
npx demo-kit validate dist/
npx vite preview            # serves dist/ at http://localhost:4173
```

To release a change, bump `version` in `package.json` in its PR (CI fails if the demo changed and that version is already released). When the PR merges, the Demo workflow builds the bundle, sets the manifest's `built` to that day, releases it as `v<version>` with `dist/` zipped as `demo-<slug>.zip`, and starts the site's build, which puts it live. For the bundle contract and the kit's API, see the [demo-kit README](https://github.com/Mark-Bugden/demo-kit#readme).

The files to edit:

| File | What to change |
| --- | --- |
| `public/manifest.json` | The demo's metadata. |
| `src/main.js` | The controls, the animation clock, the status line and the `mount()` call. |
| `src/pad.js`, `src/momenta.js` | The momentum pad, and its geometry: where the allowed region ends, and where the dot goes when the spins change. `test/momenta.test.js` checks the geometry. |
| `src/spins.js` | The spin sliders, held to a + b ≤ 1, with the "Extremal" label. |
| `src/presets.js` | The presets and the opening state. The figures' values were fitted to the paper's images by `../docs/research/scripts/presets/`. |
| `src/scene/` | The 3D view, with three.js: the spheroids, the trail, the photon, the spin arrow (`spin.js`) and the camera. |
| `src/panel/` | The side panel beside the 3D view (`panel.js`): θ, φ and ψ against λ (`strip.js`) and the (φ, ψ) torus (`torus.js`), on a 2D canvas, with their axes' ticks (`ticks.js`, checked by `test/ticks.test.js`). Also the colour legend drawn over the 3D view (`legend.js`). |
| `src/trail.js`, `src/colour.js` | The path ready to draw: `fill()` writes the 3D trail and `trace()` walks it for the panel. And the paper's colour wheel. `test/trail.test.js` checks them. |
| `src/physics/` | The orbits: a port of the Python reference model in `../src/photon_orbits/`, explained in `../docs/equations.md`. Change both together and regenerate `../tests/fixtures/orbits.json`, which `test/physics.test.js` replays. |
| `public/data/`, `public/model/` | Runtime files. |
| `public/static.webp`, `thumb.webp`, `og.jpg`, and the light versions `static-light.webp` and `thumb-light.webp` | The still images: the opening state, shown while the demo loads, and the project card and link preview, Fig 4.10's orbit. `../scripts/stills.py` renders them from the built demo: run `npm run build`, then `uv run --with playwright --with pillow python scripts/stills.py` from the repo root, then build again. |
