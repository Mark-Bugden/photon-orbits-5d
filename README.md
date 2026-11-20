# Plotting the Photons

An interactive demo of spherical photon orbits around a five-dimensional rotating (Myers–Perry) black hole, for [mark-bugden.github.io](https://mark-bugden.github.io/). It is based on the paper:

> M. Bugden, "Spherical photon orbits around a 5D Myers-Perry black hole", *General Relativity and Gravitation* **50** (2018), [doi:10.1007/s10714-018-2352-x](https://doi.org/10.1007/s10714-018-2352-x), [arXiv:1801.03389](https://arxiv.org/abs/1801.03389).

The visitor sets the hole's two spins, $a$ and $b$, and picks the photon's two angular momenta, $\Phi$ and $\Psi$, on a pad that shows which pairs have a spherical orbit. The demo then draws the orbit on its spheroid, with the fourth dimension, $\psi$, shown as colour.

**Status:** milestones 1 to 4 of [the plan](docs/plan.md) are done: the physics, whose equations are set out in [docs/equations.md](docs/equations.md), the 3D view, the momentum pad with presets for the paper's figures, and the side panel, which plots θ, φ and ψ against λ and the path on the (φ, ψ) torus. Polish and the 1.0 release come next, in milestone 5.

## Repository layout

| Path | What it is |
| --- | --- |
| [`demo/`](demo/) | The interactive demo, built with [demo-kit](https://github.com/Mark-Bugden/demo-kit). [`demo/src/physics/`](demo/src/physics/) is the JavaScript port of the reference model, and [`demo/src/scene/`](demo/src/scene/) the 3D view, drawn with [three.js](https://threejs.org/). Its README explains how to run and build the demo. |
| [`src/photon_orbits/`](src/photon_orbits/) | Python: the reference model of the orbits, which the demo's physics is a port of. |
| [`tests/`](tests/) | Python tests, and `fixtures/orbits.json`: reference results that the demo's physics must reproduce. |
| [`docs/`](docs/) | [The plan](docs/plan.md), [the equations as implemented](docs/equations.md), and [the research](docs/research/) behind them: a line-by-line review of the paper with its errata, and a survey of existing tools. |

## The reference model

With [uv](https://docs.astral.sh/uv/):

```sh
uv sync
uv run pytest                             # about 20 seconds
uv run python -m photon_orbits.fixtures   # rewrites tests/fixtures/orbits.json, after any change to the model
```

The model has no dependencies beyond Python 3.12. The tests also use numpy, scipy and sympy, as independent checks:

- the roots of the radial quartic;
- elliptic-integral periods;
- a full integration of the geodesic equations from the metric.

```python
from photon_orbits import solve_orbit, tabulate, frequencies

orbit = solve_orbit(a=0.6, b=0.3, Phi=1.0, Psi=-1.0)   # Fig 4.7 of the paper
orbit.x, orbit.K                                        # x = r² = 1.754557, K = 4.020430
table = tabulate(orbit)                                 # one latitude oscillation, in Mino time σ and λ
frequencies(orbit, table)                               # turns of φ and ψ per oscillation: (0.336346, -0.607973)
```

## The demo

```sh
cd demo
npm install
npm run dev       # live preview at http://localhost:5173
npm test          # the physics, against the reference model's fixtures, the drawn path, the pad and the panel
```

`src/photon_orbits/` and `demo/src/physics/` do the same operations in the same order. If you change one, change the other, regenerate the fixtures, and run both test suites. CI runs both on every pull request.

## Licence

[MIT](LICENSE).
