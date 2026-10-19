#!/usr/bin/env python3
"""Typeset Figure 5.3 at 170 mm width with 10 pt text and unchanged curves.

The original high-resolution raster retains every curve and DNS marker exactly.
Only the old legend, in a curve-free area, is replaced by the matching grid.
The panels are embedded unchanged in a PDF with vector ticks, labels and legend.
No curves are digitized, refitted or regenerated from a different data revision.
Run with Python, matplotlib, NumPy and Pillow; outputs go to the parent folder.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from PIL import Image
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "original_export.png"
OUTPUT = HERE.parent / "Griffin_Wallfunction_Edge_Quantities"
SOURCE_HASH = "9bae26b99dc882f786a63d7d9cd3bc953415249ac23af19dfd3f177ee0fe1533"
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_HASH
source = np.asarray(Image.open(SOURCE).convert("RGB"))
restored = source.copy()

# The old legend covers only the upper-left, curve-free region of panel 1.
# A clear row and column retain the existing grid and its antialiasing.
mx0, mx1, my0, my1 = 118, 809, 60, 216
restored[my0:my1, mx0:mx1] = np.minimum(
    source[my0:my1, 900][:, None, :],
    source[230, mx0:mx1][None, :, :])
outside = np.any(restored != source, axis=2)
outside[my0:my1, mx0:mx1] = False
assert not outside.any()

font_dir = Path('/mnt/c/texlive/2026/texmf-dist/fonts/opentype/public/tex-gyre')
family = 'DejaVu Sans'
if font_dir.is_dir():
    for face in ('regular', 'italic', 'bold', 'bolditalic'):
        font_manager.fontManager.addfont(font_dir / f'texgyreheros-{face}.otf')
    family = 'TeX Gyre Heros'
plt.rcParams.update({
    'font.family': family, 'font.size': 10, 'axes.labelsize': 10,
    'xtick.labelsize': 10, 'ytick.labelsize': 10, 'legend.fontsize': 10,
    'legend.title_fontsize': 10, 'mathtext.fontset': 'cm',
    'axes.linewidth': .65, 'savefig.facecolor': 'white',
})
fig, axes = plt.subplots(3, 1, figsize=(170 / 25.4, 200 / 25.4),
                         sharex=True, constrained_layout=True)
fig.get_layout_engine().set(rect=(0, 0, 1, .84), h_pad=.04, hspace=.04)

# Coordinates are source-image pixel centres, not inferred numerical data.
# Ticks coincide with the original grid lines; both axis ranges are retained.
x0, x1 = 112, 1971
xticks = [194.5, 456.5, 717.5, 978.5, 1240.5, 1501.5, 1763.5]
panels = [
    (53, 710, [70.5, 192.5, 314.5, 437.5, 559.5, 681.5],
     ['10', '8', '6', '4', '2', '0'], r'Reported $\delta$ [mm]'),
    (736, 1393, [790.5, 875.5, 960.5, 1045.5, 1130.5, 1215.5, 1300.5, 1384.5],
     ['860', '855', '850', '845', '840', '835', '830', '825'], r'Reported $U_e$ [m/s]'),
    (1419, 2076, [1500.5, 1645.5, 1790.5, 1935.5],
     ['80', '75', '70', '65'], r'Reported $T_e$ [K]'),
]
for ax, (y0, y1, yticks, ylabels, ylabel) in zip(axes, panels):
    crop = restored[y0:y1, x0:x1]
    extent = (x0-.5, x1-.5, y1-.5, y0-.5)
    ax.imshow(crop, extent=extent, interpolation='none', aspect='auto')
    ax.set_xlim(extent[:2]); ax.set_ylim(extent[2:])
    ax.set_xticks(xticks, [f'{v / 10:.1f}' for v in range(7)])
    ax.set_yticks(yticks, ylabels)
    ax.set_ylabel(ylabel, labelpad=5)
    ax.tick_params(length=3, width=.65, pad=3)
axes[-1].set_xlabel('Wall-face x [m]')

nominal = np.array([1, 2, 5, 10, 20, 30, 40, 60, 100, 150, 200])
colors = plt.cm.viridis(LogNorm(vmin=1, vmax=200)(nominal))
handles = [Line2D([], [], color=c, lw=1.15, label=str(n))
           for n, c in zip(nominal, colors)]
handles.extend([
    Line2D([], [], color='black', ls='--', lw=1.6, label='Wall-resolved SA 5e-4_dx5'),
    Line2D([], [], color='tab:red', ls='none', marker='*', ms=8,
           label='DNS, Griffin edge'),
    Line2D([], [], color='tab:red', ls='none', marker='D', ms=5,
           markerfacecolor='none', label=r'DNS at reported $\delta=3.6$ mm'),
])
legend = fig.legend(handles=handles, loc='upper center', ncol=4,
                    bbox_to_anchor=(.5, 1), frameon=False,
                    title=r'Nominal $y^+$ / reference', handlelength=1.6,
                    handletextpad=.5, columnspacing=1.1, labelspacing=.25,
                    borderaxespad=.2)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
visible = [legend.get_title(), *legend.get_texts()]
for ax in axes:
    visible.extend([ax.xaxis.label, ax.yaxis.label,
                    *ax.get_xticklabels(), *ax.get_yticklabels()])
for text in visible:
    if not text.get_visible() or not text.get_text():
        continue
    box = text.get_window_extent(renderer)
    assert box.x0 >= 0 and box.y0 >= 0, (text.get_text(), box)
    assert box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height, (text.get_text(), box)
    assert text.get_fontsize() == 10
assert not legend.get_window_extent(renderer).overlaps(axes[0].get_tightbbox(renderer))
for extension in ('pdf', 'png'):
    fig.savefig(OUTPUT.with_suffix('.' + extension), dpi=300)
plt.close(fig)
(HERE / 'font_layout_audit.json').write_text(json.dumps({
    'method': 'Retained raster curves with new vector typography; no numerical replot.',
    'original_export': '/home/mgnadlinger/Xue_et_al/flat_plate/M6Tw025/validation/outputs_Griffin_SA/reported_edge_quantities.png',
    'original_generator': '/home/mgnadlinger/WF_validation/LowH_M3/Griffin_SA/analyze_griffin_edge_diagnostics.py',
    'source_sha256': SOURCE_HASH,
    'source_size_pixels': list(Image.open(SOURCE).size),
    'width_mm': 170, 'height_mm': 200, 'base_font_pt': 10,
    'font_family': family,
    'pixels_changed_outside_old_legend': int(outside.sum()),
    'old_legend_rectangle_xyxy': [mx0, my0, mx1, my1],
    'original_panel_crops_xyxy': [[x0, p[0], x1, p[1]] for p in panels],
    'all_labels_inside_canvas': True, 'legend_clear_of_panels': True,
    'source_reconstruction_check': 'Original generator and relocated data were found. Excluding flagged diagnostics nearly reproduced the export, but exact curve equivalence was not established; the retained raster is therefore used unchanged.',
}, indent=2) + '\n')
print('Saved Figure 5.3 at 170 mm with 10 pt text and preserved raster curves.')
