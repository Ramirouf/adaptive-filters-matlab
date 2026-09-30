"""Read the curve markers off Figs. 3a and 3b of the paper.

The paper's figures put markers on the fKF, sKF and KF curves (SG has none). This script
renders page 6 of the publisher PDF at 400 dpi, finds the markers by colour, and converts
their pixel centres to (t [s], misalignment [dB]) with the figure gridlines. The output,
data/paper_fig3_markers.csv, is what verify/compare_paper.m checks our curves against.

The PDF is not in this repository (it is copyrighted); pass its path:
    python data/digitize_fig3.py /path/to/paper.pdf
Needs pdftoppm (poppler), numpy, scipy and pillow.

Reading resolution: one pixel is 0.065 dB vertically, and 0.013 s (Fig. 3a) or 0.0013 s
(Fig. 3b) horizontally.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

DPI = 400
# crop: panel box as fractions of the page (left, top, right, bottom).
# x0, xs: pixel column of t = 0 and pixels per second; y0, ys: pixel row of 0 dB and pixels
# per dB. Measured on the gridlines of the 400 dpi crop (t = 0, 2, ..., 12 s or 0, 0.2, ...;
# 0, -10, -20 dB). legend: box to ignore, in crop pixels. steep_until: the KF stars are
# found by how wide the cyan trace is, which only works once the curve has flattened.
PANELS = {
    "fig3a": dict(crop=(0.52, 0.035, 0.975, 0.215), x0=221.5, xs=157.9 / 2.0, y0=192.0,
                  ys=15.4, legend=(840, 170, 1160, 410), steep_until=1.0),
    "fig3b": dict(crop=(0.52, 0.215, 0.975, 0.39), x0=221.5, xs=157.9 / 0.2, y0=91.0,
                  ys=15.45, legend=(840, 70, 1160, 310), steep_until=0.17),
}
COLOR = {"fKF": (251, 0, 0), "sKF": (191, 0, 190), "KF": (10, 192, 191)}


def in_box(x, y, box):
    return box[0] < x < box[2] and box[1] < y < box[3]


def filled_markers(im, rgb, p):
    """fKF circles and sKF crosses: blobs that survive an erosion thicker than the line."""
    mask = np.abs(im - np.array(rgb)).sum(2) < 60
    lab, n = ndimage.label(ndimage.binary_erosion(mask, np.ones((5, 5))))
    out = []
    for i in range(1, n + 1):
        ys, xs = np.where(lab == i)
        if ys.size >= 15 and not in_box(xs.mean(), ys.mean(), p["legend"]):
            out.append((xs.mean(), ys.mean()))
    return out


def star_markers(im, rgb, p):
    """KF stars: columns where the cyan trace is much wider than the 3-pixel line."""
    mask = np.abs(im - np.array(rgb)).sum(2) < 90
    mask[:, :int(p["x0"] + p["steep_until"] * p["xs"])] = False
    cols = np.where(mask.sum(0) >= 10)[0]
    out = []
    for g in np.split(cols, np.where(np.diff(cols) > 3)[0] + 1):
        if g.size == 0:
            continue
        ys = np.where(mask[:, g].any(1))[0]
        x, y = g.mean(), 0.5 * (ys.min() + ys.max())
        if ys.max() - ys.min() <= 20 and not in_box(x, y, p["legend"]):
            out.append((x, y))
    return out


def main(pdf):
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdftoppm", "-f", "6", "-l", "6", "-r", str(DPI), "-png", str(pdf),
                        f"{tmp}/page"], check=True)
        page = Image.open(next(Path(tmp).glob("page*.png"))).convert("RGB")
    W, H = page.size
    rows = []
    for panel, p in PANELS.items():
        c = p["crop"]
        im = np.asarray(page.crop((int(c[0] * W), int(c[1] * H), int(c[2] * W), int(c[3] * H))))
        im = im.astype(int)
        for variant, rgb in COLOR.items():
            pts = star_markers(im, rgb, p) if variant == "KF" else filled_markers(im, rgb, p)
            for x, y in sorted(pts):
                rows.append((panel, variant, (x - p["x0"]) / p["xs"], -(y - p["y0"]) / p["ys"]))

    out = Path(__file__).resolve().parent / "paper_fig3_markers.csv"
    with open(out, "w") as f:
        f.write("% Markers read off Figs. 3a and 3b of the paper (p. 6) by data/digitize_fig3.py.\n")
        f.write("% panel, variant, t_sec, m_db   (m_db = mean misalignment in dB; SG has no markers)\n")
        for panel, variant, t, m in rows:
            f.write(f"{panel}, {variant}, {t:.4f}, {m:.2f}\n")
    print(f"wrote {out}: {len(rows)} markers")


if __name__ == "__main__":
    main(sys.argv[1])
