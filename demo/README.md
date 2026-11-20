# Demo

The interactive demo for this project. It is built with [demo-kit](https://github.com/Mark-Bugden/demo-kit), which is pinned in `package.json`.

## Run it locally

```sh
npm install
npm run dev
npm test          # node:test: the physics, against the reference model's fixtures, and the drawn path
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
| `src/scene/` | The 3D view, with three.js: the spheroids, the trail, the photon and the camera. |
| `src/trail.js`, `src/colour.js` | The path ready to draw, and the paper's colour wheel. `test/trail.test.js` checks them. |
| `src/physics/` | The orbits: a port of the Python reference model in `../src/photon_orbits/`, explained in `../docs/equations.md`. Change both together and regenerate `../tests/fixtures/orbits.json`, which `test/physics.test.js` replays. |
| `public/data/`, `public/model/` | Runtime files. |
| `public/static.webp`, `thumb.webp`, `og.jpg` | Images. |
