"""Replot Figure 4.32 from retained original CSV with 10% larger text.

Original source: /home/mgnadlinger/Cao_lamZPGFP/validation/plot_first_cell_acoustic_re.py
CSV: validation/mesh_sensitivity_output_1550/1550_dx_1p5e-1_wall_normal_first_cell_acoustic_re.csv
The snapshot was verified byte-identical by tmp/notation_figures_20261019/cao_audit.json.
Original plotted-coordinate hashes are checked before exporting. Keep the
2274 x 1065 canvas and manuscript inclusion width, so font scaling is exact.
Requires NumPy, Matplotlib and LaTeX with lmodern and dvipng on PATH.
"""
from pathlib import Path
import csv
import hashlib
import json
import os
if os.name == 'nt':
    os.environ['PATH'] = r'C:\texlive\2026\bin\windows' + os.pathsep + os.environ['PATH']
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCES = HERE / 'acoustic_re_sources'
SCALE = 1.10
DPI = 220
original = json.loads((SOURCES / 'original_audit.json').read_text())
with (SOURCES / 'data.csv').open(newline='') as stream:
    rows = list(csv.DictReader(stream))
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Computer Modern Roman'],
    'text.usetex': True, 'text.latex.preamble': r'\usepackage{lmodern}',
    'font.size': 12*SCALE, 'axes.labelsize': 13*SCALE,
    'xtick.labelsize': 11*SCALE, 'ytick.labelsize': 11*SCALE,
    'legend.fontsize': 11*SCALE, 'lines.linewidth': 1.8,
    'axes.linewidth': .8, 'axes.grid': True, 'grid.alpha': .24,
    'grid.linewidth': .7, 'figure.facecolor': 'white', 'savefig.facecolor': 'white',
})
fig, ax = plt.subplots(figsize=(2274/DPI, 1065/DPI), dpi=DPI)
colors = [plt.get_cmap('viridis')(x) for x in np.linspace(.15, .82, 4)]
for case, color, record in zip(dict.fromkeys(row['case'] for row in rows), colors, original['curves']):
    selected = [row for row in rows if row['case'] == case]
    ax.semilogy([float(row['x_over_L']) for row in selected],
                [float(row['Re_delta_a_mean']) for row in selected],
                color=color, label=record['label'])
ax.axvline(1., color='0.4', linestyle=':', linewidth=1.)
ax.set_xlim(original['axes'][0]['xlim'])
ax.set_ylim(original['axes'][0]['ylim'])
ax.set_xlabel(r'$x/L$ (corner at $x/L=1$)')
ax.set_ylabel(r'$Re_{\mathrm{cell1}}$')
ax.minorticks_on()
ax.grid(True, which='minor', alpha=.10, linewidth=.5)
ax.legend(loc='upper left', framealpha=.95, labelspacing=.2, borderpad=.25)
fig.subplots_adjust(left=.078, right=.986, bottom=.132, top=.992)
hashes = [hashlib.sha256(np.column_stack((line.get_xdata(), line.get_ydata())).astype('<f8').tobytes()).hexdigest() for line in ax.lines]
assert hashes == [entry['sha256'] for entry in original['curves']], 'Plotted data changed'
for extension in ('png', 'pdf'):
    fig.savefig(HERE / ('05_first_cell_acoustic_re.' + extension), dpi=DPI)
(SOURCES/'font_update.json').write_text(json.dumps({
    'font_scale': SCALE, 'axis_label_pt': 13*SCALE, 'tick_and_legend_pt': 11*SCALE,
    'canvas_pixels': [2274, 1065], 'plotted_data_hashes_match_original': True,
    'data_sha256': hashlib.sha256((SOURCES/'data.csv').read_bytes()).hexdigest(),
    'source_script': original['original_script'], 'matplotlib_version': matplotlib.__version__,
}, indent=2)+'\n')
print('Figure 4.32 regenerated; all original plotted-coordinate hashes match.')
