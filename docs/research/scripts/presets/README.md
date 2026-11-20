# Fitting the presets to the paper's figures

The demo's presets for Figs 4.2–4.12 (`demo/src/presets.js`) were fitted to the paper's own images, so each one reproduces its figure: the same stretch of path, the same colours and the same viewpoint. These scripts did the fitting.

They read the figures as PNGs in `../figs/`, as the other research scripts do. Those images come from the paper's arXiv source and aren't committed.

## What is fitted

The demo draws the path without the paper's constant twist arctan(a/r), with ψ = 0 and φ = 0 at θ_min (χ = 0). For each figure the fit finds:

- **the stretch of path:** the angle variable χ at its start, and the number of oscillations shown;
- **the camera:** its elevation, and its azimuth about the z axis. Turning the path about z is the same as turning the camera, so the azimuth also takes up the start's φ and the paper's twist;
- **the colour offset ψ₀:** ψ at χ = 0.

The camera is Mathematica's: perspective, 3.38 box lengths from the centre, z up.

## The scripts

Run them from the repository's root, in this order.

| Script | What it does |
| --- | --- |
| `measure.py` | Finds each figure's sphere outline, the equator ring's front arc (its depth gives the elevation), and the coloured path's pixels and hues. Writes `measured.json`. |
| `fit.py` | For figures showing a few oscillations: a grid search over start, length, elevation and azimuth, then Nelder–Mead. The score is the mean distance from our path to the paper's pixels and back, plus how much the colour offset varies along the path. |
| `fast.py` | For figures showing many oscillations. Oscillation k is oscillation 0 turned by k times the turn of φ per oscillation, so one oscillation's projection is computed once for every azimuth, and any stretch scores as a sum of cumulative sums. The best few of each length are then rescored exactly and refined. |
| `derive.py` | Prints the presets in the demo's terms, and cross-checks Figs 4.8 and 4.10 against the starts fitted from Figs 4.9 and 4.11 (`../fit_figs.out`). |

```sh
uv run --with pillow python docs/research/scripts/presets/measure.py
C='[0,0.393,0.785,1.178,1.571,1.963,2.356,2.749,3.142,3.534,3.927,4.320,4.712,5.105,5.498,5.890]'
C2='[0,0.785,1.571,2.356,3.142,3.927,4.712,5.498,6.283,7.069,7.854,8.639,9.425,10.21,10.996,11.781]'
uv run --with pillow --with matplotlib python docs/research/scripts/presets/fit.py \
  "{\"4.2\": {\"chis\":$C2, \"oscs\":[1.75,2,2.25], \"els\":[14,18,22], \"free\":[0,1,2,3], \"max_osc\": 5},
    \"4.4\": {\"chis\":$C, \"oscs\":[1.5,1.75,2,2.25,2.5], \"els\":[30,34,38], \"free\":[0,1,2,3], \"max_osc\": 5},
    \"4.6\": {\"chis\":$C2, \"oscs\":[1.5,1.75,2], \"els\":[-3,0,3], \"free\":[0,1,2,3], \"max_osc\": 5},
    \"4.10\": {\"chis\":$C, \"oscs\":[1.25,1.5,1.75], \"els\":[26,30,34], \"free\":[0,1,2,3], \"max_osc\": 5},
    \"4.12\": {\"chis\":[0], \"oscs\":[0.75,1,1.25,1.5,1.75], \"els\":[28,32,36], \"free\":[1,2,3], \"max_osc\": 4}}"
uv run --with pillow --with matplotlib python docs/research/scripts/presets/fast.py \
  '{"4.3": {"els":[10,12,14,16], "min_osc":45, "max_osc":50, "chi_steps":96},
    "4.5": {"els":[33,37,41,45], "min_osc":8, "max_osc":35},
    "4.7": {"els":[22,26,30,34], "min_osc":8, "max_osc":35},
    "4.8": {"els":[26,30,34,38], "min_osc":1, "max_osc":9}}'
uv run --with pillow --with matplotlib python docs/research/scripts/presets/derive.py
```

Each fit writes `fit-<fig>.json` and `fit-<fig>.png` here: the paper, our path, and the two overlaid. [`fits.json`](fits.json) collects the results that `presets.js` uses.

## From the fits to the presets

The presets take the fitted values, with these exceptions:

- **Fig 4.2** starts at exactly χ = π, at the equator, and shows exactly two oscillations (fitted: 3.19 and 2.006).
- **Fig 4.6** shows two oscillations from the pole, the caption's "one latitudinal oscillation". The fit can't fix the length, since the rest lies behind the spheroid.
- **Fig 4.10** takes its start and colours from the Fig 4.11 fit, which plots the same stretch against λ: λ from 1 to 5. That start is χ = 4.242 (the 3D fit's is 4.212), and its ψ there is 2.957 (the 3D fit's is 2.978).
- **Fig 4.12** shows one full turn of φ. The fit only sees the front of the ring.
- **The side panel** (milestone 4) reads like Figs 4.9 and 4.11 through two more values, from `../fit_figs.out`: λ₀, the demo's λ at the paper's λ = 0, and φ₀, φ at θ_min. Fig 4.10 has λ₀ = 0.7509 and φ₀ = −1.5887, and its ψ₀ already matches Fig 4.11. Fig 4.8 has λ₀ = 0.076 and φ₀ = 1.6122, and keeps the ψ₀ of its own colours, 3.683, so its panel's ψ is Fig 4.9's raised by about π.

## What the fits show about the figures

- **The twist is real.** Fig 4.10's box shows its camera azimuth to be about −12°. The Fig 4.11 fit's φ, with our fitted azimuth, puts it at 347.5° when the arctan(a/r) twist is included, and at 26.1° without it.
- **Fig 4.8 doesn't share Fig 4.9's start in φ and ψ.** Its stretch starts at θ_min, as Fig 4.9's fit says, but its colours are about π from Fig 4.9's ψ curve (3.72 against 0.45), and its box doesn't match the azimuth that Fig 4.9's φ implies. So the preset follows Fig 4.8.
- **In Fig 4.3 the drawn sphere is smaller than the orbit's spheroid.** The red net bulges past the black outline at the equator. Its length was fixed by counting strands: about 12.7 cross each row across the middle of the band, which takes about 47–50 oscillations.
