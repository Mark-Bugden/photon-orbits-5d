"""Renders the demo's still images from the built demo, in a real browser.

- static.webp: the opening state, the zero momenta orbit with the side panel beside it, shown
  while the demo loads and without JavaScript.
- thumb.webp and og.jpg: the 3D view alone, on CARD's orbit and trail, at the default zoom, for the
  project card and link previews. The page's address ends in ?card, so the demo leaves out the spin
  arrow.
- static-light.webp and thumb-light.webp: the same in the light theme, which the site shows under
  its light theme. og.jpg stays dark only: link previews have no theme.

It serves demo/dist/ locally, drives it with Playwright, and writes the images to demo/public/,
fitted onto their theme's page colour and kept under the bundle contract's size limits. The demo
starts paused under reduced motion, so every still is of the same moment. Build first:

    cd demo && npm run build && cd ..
    uv run --with playwright --with pillow python scripts/stills.py            # Microsoft Edge
    uv run --with playwright --with pillow python scripts/stills.py chrome     # or another channel

Then rebuild, so dist/ picks up the new images.
"""

import functools
import io
import math
import sys
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from PIL import Image
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "demo" / "dist"
OUT = ROOT / "demo" / "public"
THEMES = ("dark", "light")
BACKGROUND = {"dark": (10, 14, 20), "light": (245, 246, 248)}  # --color-bg in each theme
# The card's orbit: a preset, by its name in the drop-down, and the trail's length in λ.
CARD = {"preset": "Equal spins (a = b)", "trail": 100}
# Wide enough for the side panel to sit beside the 3D view; the stage is then about 860 px wide.
VIEWPORT = {"width": 1200, "height": 1000}
SCALE = 2
READY = "/r = /.test(document.querySelector('.dk-status')?.textContent)"
# Every still: the stage at the image's shape, without the zoom buttons.
STILL_STYLE = (
    "[data-demo] .dk-stage { height: auto !important; aspect-ratio: var(--still-shape, 16 / 9) !important; }"
    " [data-demo] .po-view-buttons { display: none !important; }"
)
# The card: the 3D view alone, without the side panel or the colour key. ?card leaves out the spin arrow.
CARD_STYLE = "[data-demo] .po-panel, [data-demo] .po-legend { display: none !important; }"


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def serve():
    handler = functools.partial(QuietHandler, directory=str(DIST))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server


def fit(png, width, height, theme):
    image = Image.open(io.BytesIO(png)).convert("RGB")
    scale = min(width / image.width, height / image.height)
    image = image.resize((round(image.width * scale), round(image.height * scale)), Image.LANCZOS)
    canvas = Image.new("RGB", (width, height), BACKGROUND[theme])
    canvas.paste(image, ((width - image.width) // 2, (height - image.height) // 2))
    return canvas


def save(image, name, fmt, limit, **options):
    for quality in range(90, 40, -5):
        buffer = io.BytesIO()
        image.save(buffer, fmt, quality=quality, **options)
        if limit is None or buffer.tell() < limit:
            break
    else:
        raise SystemExit(f"{name}: could not get under {limit / 1000:.0f} KB")
    (OUT / name).write_bytes(buffer.getvalue())
    print(f"{name}: {image.width}×{image.height}, quality {quality}, {buffer.tell() / 1000:.1f} KB")


def open_page(browser, url, style):
    page = browser.new_page(viewport=VIEWPORT, device_scale_factor=SCALE, reduced_motion="reduce")
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(url)
    page.add_style_tag(content=STILL_STYLE + style)
    page.wait_for_function(READY, timeout=30000)
    return page, errors


def shoot(page, themes=THEMES):
    """Screenshots the stage in each theme: {theme: png}."""
    shots = {}
    for theme in themes:
        page.evaluate(f"document.documentElement.dataset.theme = '{theme}'")
        page.wait_for_timeout(500)  # the view repaints in the theme's colours
        shots[theme] = page.locator(".dk-stage").screenshot()
    return shots


def main():
    channel = sys.argv[1] if len(sys.argv) > 1 else "msedge"
    if not (DIST / "index.html").exists():
        raise SystemExit("No demo/dist/: run `npm run build` in demo/ first")
    server = serve()
    url = f"http://127.0.0.1:{server.server_port}/"
    with sync_playwright() as p:
        browser = p.chromium.launch(channel=channel)

        # The opening state, as a visitor first sees it.
        page, errors = open_page(browser, url, "")
        opening = shoot(page)
        page.close()

        # The card: CARD's orbit and trail.
        page, more = open_page(browser, url + "?card", CARD_STYLE)
        errors += more
        page.get_by_label("Preset").select_option(label=CARD["preset"])
        page.get_by_label("Trail length").evaluate(
            "(el, value) => { el.value = value; el.dispatchEvent(new Event('input', { bubbles: true })); }",
            math.log10(CARD["trail"]),
        )
        card = shoot(page)
        page.evaluate("document.documentElement.style.setProperty('--still-shape', '1200 / 630')")
        og = shoot(page, ("dark",))["dark"]
        browser.close()
    server.shutdown()
    if errors:
        raise SystemExit(f"The demo reported errors: {errors}")

    for theme in THEMES:
        suffix = "" if theme == "dark" else "-light"
        save(fit(opening[theme], 1200, 675, theme), f"static{suffix}.webp", "WEBP", None, method=6)
        save(fit(card[theme], 1200, 675, theme), f"thumb{suffix}.webp", "WEBP", 80_000, method=6)
    save(fit(og, 1200, 630, "dark"), "og.jpg", "JPEG", 150_000, optimize=True, progressive=True)


if __name__ == "__main__":
    main()
