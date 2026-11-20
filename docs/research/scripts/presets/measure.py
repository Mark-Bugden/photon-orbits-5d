"""Measure each paper figure: the sphere's outline, the equator ring (so the elevation), and the coloured path pixels."""

import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageOps
from scipy import ndimage

FIGS = Path(__file__).resolve().parent.parent / "figs"
OUT = Path(__file__).resolve().parent
NAMES = {"4.2": "nomomenta", "4.3": "phiequals2", "4.4": "firstcolour", "4.5": "manycolour", "4.6": "football",
         "4.7": "band", "4.8": "cusp", "4.10": "confused", "4.12": "boring"}
SCALE = 0.5


def load(name):
    image = ImageOps.flip(Image.open(FIGS / f"{name}.png").convert("RGB"))
    image = image.resize((round(image.width * SCALE), round(image.height * SCALE)), Image.LANCZOS)
    return np.asarray(image).astype(float) / 255


def hsv(rgb):
    mx = rgb.max(axis=2)
    mn = rgb.min(axis=2)
    s = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1e-9), 0)
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    d = np.maximum(mx - mn, 1e-9)
    h = np.where(mx == r, ((g - b) / d) % 6, np.where(mx == g, (b - r) / d + 2, (r - g) / d + 4)) / 6
    return h, s, mx


def measure(fig):
    rgb = load(NAMES[fig])
    h, s, v = hsv(rgb)
    black = v < 0.16
    path = (s > 0.5) & (v > 0.6)
    light = (rgb.min(axis=2) > 0.6) & (s < 0.2)
    blob = black | path
    labels, count = ndimage.label(ndimage.binary_closing(blob, iterations=2))
    sizes = ndimage.sum(blob, labels, range(1, count + 1))
    sphere = labels == (1 + int(np.argmax(sizes)))
    sphere = ndimage.binary_fill_holes(sphere)
    rows, cols = np.nonzero(sphere)
    cx = 0.5 * (cols.min() + cols.max())
    cy = 0.5 * (rows.min() + rows.max())
    rx = 0.5 * (cols.max() - cols.min())
    ry = 0.5 * (rows.max() - rows.min())
    # The equator ring's front arc: light pixels inside the sphere's box, below the centre, near the middle column.
    inside = ndimage.binary_dilation(sphere, iterations=3)
    ring = light & inside
    band = ring[:, int(cx) - 3:int(cx) + 4]
    ring_rows = np.nonzero(band.any(axis=1))[0]
    below = ring_rows[ring_rows > cy]
    above = ring_rows[ring_rows < cy]
    ring_low = float(below.mean()) if len(below) else None
    ring_high = float(above.mean()) if len(above) else None
    py, px = np.nonzero(path & inside)
    return {
        "fig": fig, "size": [rgb.shape[1], rgb.shape[0]], "cx": cx, "cy": cy, "rx": rx, "ry": ry,
        "ring_low": ring_low, "ring_high": ring_high,
        "path": {"x": px.tolist(), "y": py.tolist(), "hue": h[py, px].tolist()},
    }


if __name__ == "__main__":
    out = {}
    for fig in (sys.argv[1:] or NAMES):
        m = measure(fig)
        el = math.degrees(math.asin(min(1, (m["ring_low"] - m["cy"]) / m["rx"]))) if m["ring_low"] else None
        print(f"Fig {fig}: centre ({m['cx']:.1f}, {m['cy']:.1f}), semi-axes {m['rx']:.1f} x {m['ry']:.1f}, "
              f"ring front at {m['ring_low']}, back at {m['ring_high']}, elevation ~ {el}, path px {len(m['path']['x'])}")
        out[fig] = m
    (OUT / "measured.json").write_text(json.dumps(out))
