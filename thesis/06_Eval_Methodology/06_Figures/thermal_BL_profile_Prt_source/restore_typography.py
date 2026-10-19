#!/usr/bin/env python3
"""Restore legible typography while preserving the archived plot's curve pixels.

This is a graphic reconstruction, not a replot from numerical data. The source
is the original high-resolution PNG embedded in slide 21 of Weekly 260803.pptx.
The Pr_t sensitivity study's plotting script and raw profile data were not
located. All curve pixels are retained; only the old legend regions (which do
not intersect the curves) are replaced by the matching background grid.

Run: python3 restore_typography.py
Dependencies: numpy, Pillow, matplotlib. Outputs are next to this directory.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from PIL import Image
import numpy as np

HERE = Path(__file__).resolve().parent
OUTPUT = HERE.parent / "thermal_BL_profile_Prt"
SOURCE = HERE / "original_export.png"
EXPECTED_SHA256 = "9708da54cd3e973f5ca2d962801fb3143dac98e87d9ea43852aa47e65a63a297"
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED_SHA256
source = np.asarray(Image.open(SOURCE).convert("RGB"))

# Original panel frames and tick-centre geometry, measured in export pixels.
# The image artist maps original pixel centres back into the original axes.
# At y=1062: T/T_inf=1.0; at y=253: T/T_inf=2.75.
y0, y1 = 160, 1104
pixel_y = lambda y: 1.0 + (1062.0 - y) * 1.75 / 809.0
panels = [
    dict(model="Spalart–Allmaras", x0=170, x1=1618, tick0=189, tick1=1584,
         mask=(180, 690, 865, 1085), blank_column=710,
         retau=[467.0, 456.3, 447.9], pos=(0.125, 0.17, 0.39, 0.73)),
    dict(model=r"$k$–$\omega$ SST", x0=1717, x1=3165, tick0=1736, tick1=3131,
         mask=(1728, 2240, 865, 1085), blank_column=2260,
         retau=[465.5, 459.6, 455.1], pos=(0.59, 0.17, 0.39, 0.73)),
]
colors = ["black", "#0072b2", "#e69f00", "#cc79a7"]
# Requested reduction from the previous 11/12 pt typography.
FONT_SCALE = 0.70
TEXT_FONT_PT = round(11 * FONT_SCALE, 2)
TITLE_FONT_PT = round(12 * FONT_SCALE, 2)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": TEXT_FONT_PT,
    "axes.titlesize": TITLE_FONT_PT, "axes.labelsize": TEXT_FONT_PT,
    "xtick.labelsize": TEXT_FONT_PT, "ytick.labelsize": TEXT_FONT_PT,
    "legend.fontsize": TEXT_FONT_PT, "pdf.fonttype": 42,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
})
fig = plt.figure(figsize=(170 / 25.4, 75 / 25.4), facecolor="white")
audit = []
for index, p in enumerate(panels):
    restored = source.copy()
    mx0, mx1, my0, my1 = p["mask"]
    # These strips are outside all curves and labels. Combining their gray
    # levels reconstructs the existing white background and grid intersections.
    # No data curve lies inside these original legend rectangles.
    background_row = source[1090, mx0:mx1]
    background_col = source[my0:my1, p["blank_column"]]
    replacement = np.minimum(background_col[:, None, :], background_row[None, :, :])
    restored[my0:my1, mx0:mx1] = replacement
    changed = np.any(restored != source, axis=2)
    outside = changed.copy()
    outside[my0:my1, mx0:mx1] = False
    assert not outside.any(), "Pixels outside the old legend must remain exact."
    crop = restored[y0:y1, p["x0"]:p["x1"]]
    xcoord = lambda x: (x - p["tick0"]) / (p["tick1"] - p["tick0"])
    extent = (xcoord(p["x0"] - .5), xcoord(p["x1"] - .5),
              pixel_y(y1 - .5), pixel_y(y0 - .5))
    ax = fig.add_axes(p["pos"])
    ax.imshow(crop, extent=extent, origin="upper", interpolation="none", aspect="auto", zorder=0)
    ax.set_xlim(extent[:2]); ax.set_ylim(extent[2:])
    ax.set_xticks(np.arange(0, 1.01, .2)); ax.set_yticks(np.arange(1.0, 2.76, .25))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.tick_params(axis="both", pad=3, length=3)
    ax.set_xlabel(r"$U/U_\infty$", labelpad=4)
    if index == 0:
        ax.set_ylabel(r"$T/T_\infty$", labelpad=4)
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.set_title(p["model"], pad=6)
    labels = [r"DNS ($Re_\tau=450$)"] + [
        rf"$Pr_t={prt:g}$ ($Re_\tau={retau:.1f}$)"
        for prt, retau in zip((.6, .8, 1.), p["retau"])
    ]
    handles = [Line2D([], [], color=c, lw=1.5) for c in colors]
    legend = ax.legend(handles, labels, loc="lower left", bbox_to_anchor=(.02, .01),
                       frameon=False, borderpad=0, borderaxespad=0,
                       handlelength=1.35, handletextpad=.5, labelspacing=.24)
    audit.append({"model": p["model"], "original_crop_xyxy": [p["x0"], y0, p["x1"], y1],
                  "old_legend_replaced_xyxy": [mx0, my0, mx1, my1],
                  "pixels_changed_outside_legend": int(outside.sum()),
                  "retau": p["retau"], "legend_base_font_pt": TEXT_FONT_PT})
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for text in fig.findobj(matplotlib.text.Text):
    if not text.get_visible() or not text.get_text():
        continue
    bounds = text.get_window_extent(renderer)
    assert bounds.x0 >= -0.1 and bounds.y0 >= -0.1, (text.get_text(), bounds)
    assert bounds.x1 <= fig.bbox.x1 + 0.1 and bounds.y1 <= fig.bbox.y1 + 0.1, (text.get_text(), bounds)
for ext in ("pdf", "png"):
    options = {"dpi": 300} if ext == "png" else {}
    fig.savefig(OUTPUT.with_suffix('.' + ext), **options)
plt.close(fig)
manifest = {
    "method": "Typography restoration using unchanged raster curves; not a numerical-data replot",
    "source_presentation": "/mnt/c/Users/MichaelGnadlinger/OneDrive - Hypersonica/Hypersonica Team - Documents/02 Technical/03 Flight physics/16 Theses/2026/01 Gnadl/Docs_and_Presentations/Weekly 260803.pptx",
    "source_slide": 21, "source_zip_member": "ppt/media/image13.png",
    "source_png_sha256": EXPECTED_SHA256, "source_size_pixels": [3189, 1232],
    "preserved_previous_manuscript_png_sha256": hashlib.sha256((HERE / "previous_manuscript.png").read_bytes()).hexdigest(),
    "raw_data_status": "Original Pr_t=0.6/0.8/1.0 plotting script and profile data not located in accessible WSL/Windows sources; no curves digitized, refitted or substituted.",
    "output_width_mm": 170, "output_height_mm": 75,
    "font_points": {"axis_labels": TEXT_FONT_PT, "tick_labels": TEXT_FONT_PT, "legend": TEXT_FONT_PT, "panel_titles": TITLE_FONT_PT},
    "font_revision": "User-requested 30% reduction from the previous 11 pt text / 12 pt titles; figure dimensions and curves unchanged.",
    "math_note": "Math subscripts retain conventional smaller size.",
    "panels": audit,
}
(HERE / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(f"Saved {OUTPUT.with_suffix('.pdf')} and PNG; ordinary labels {TEXT_FONT_PT:g} pt and titles {TITLE_FONT_PT:g} pt at 170 mm.")
