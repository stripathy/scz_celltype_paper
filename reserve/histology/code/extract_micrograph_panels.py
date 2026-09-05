#!/usr/bin/env python3
"""Extract the four lipofuscin-suppressed representative composites embedded in
results/microscopy/RNAscope_fig_representative_inkscape_SJT.svg as PNGs, in
display orientation, so the R figure script can draw them with house-style
labels and the manual cell markers (results/microscopy/marker_coordinates.csv).

The SVG stores each composite with a vertical-flip transform
(matrix(1,0,0,-1,...)), so the raw base64 PNGs are upside down relative to
what Inkscape displays; we flip them back here. Panel order follows the SVG
layout and marker_coordinates.csv: 0 = Control L2/3, 1 = SCZ L2/3,
2 = Control L5/6, 3 = SCZ L5/6 (donors 1047 and 1211; sites 5/6 and 12/17).

Usage (from histology/):
    python3 code/extract_micrograph_panels.py
Outputs:
    results/microscopy/panels/panel{0..3}_{group}_{layer}.png
"""

import base64
import io
import re
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
MICRO = HERE.parent / "results" / "microscopy"
SVG = MICRO / "RNAscope_fig_representative_inkscape_SJT.svg"
OUT = MICRO / "panels"
NAMES = ["Control_L23", "SCZ_L23", "Control_L56", "SCZ_L56"]


def main() -> None:
    OUT.mkdir(exist_ok=True)
    svg = SVG.read_text()
    tags = re.findall(r"<image[^>]*?>", svg, flags=re.S)
    imgs = []
    for tag in tags:
        m = re.search(r'xlink:href="data:image/png;base64,\s*([^"]+)"', tag)
        if not m:
            continue
        x = float(re.search(r' x="([^"]+)"', tag).group(1))
        y = float(re.search(r' y="([^"]+)"', tag).group(1))
        imgs.append((y, x, base64.b64decode(m.group(1))))
    # SVG y is negative and grows more negative for lower rows; sort rows top
    # to bottom (less negative first), then columns left to right.
    imgs.sort(key=lambda t: (-t[0], t[1]))
    assert len(imgs) == 4, f"expected 4 embedded images, found {len(imgs)}"
    for i, (_, _, data) in enumerate(imgs):
        im = Image.open(io.BytesIO(data)).convert("RGB")
        im = im.transpose(Image.FLIP_TOP_BOTTOM)   # undo the SVG flip
        p = OUT / f"panel{i}_{NAMES[i]}.png"
        im.save(p)
        print(p.name, im.size)


if __name__ == "__main__":
    main()
