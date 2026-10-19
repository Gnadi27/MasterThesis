#!/usr/bin/env python3
r"""Replot the two detailed wall-resolved profile figures at thesis print size.

Run with Python, matplotlib and numpy. The adjacent JSON cache contains exact
plotted float64 coordinates captured from the original WSL plotting scripts;
source hashes, station metadata and pixel-exact source/crop checks are retained
alongside it. No simulation, normalization or profile interpolation is changed.

The PDF width is exactly 170 mm. With width=\linewidth in the thesis, all base
text, ticks and legends remain 10.5 pt. PNG previews are also written. Only the
upper two panels of the original LowH_M3 figure are included, as in the thesis.
"""
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker
import numpy as np

HERE = Path(__file__).resolve().parent
FONT_SIZE_PT = 10.5
WIDTH_MM = 170.0


def configure_fonts(font_dir):
    if font_dir.is_dir():
        for face in ('regular', 'italic', 'bold', 'bolditalic'):
            font_manager.fontManager.addfont(font_dir / f'texgyreheros-{face}.otf')
        family = 'TeX Gyre Heros'
    else:
        family = 'DejaVu Sans'
    plt.rcParams.update({
        'font.family': family, 'font.size': FONT_SIZE_PT,
        'axes.labelsize': FONT_SIZE_PT, 'axes.titlesize': FONT_SIZE_PT,
        'xtick.labelsize': FONT_SIZE_PT, 'ytick.labelsize': FONT_SIZE_PT,
        'legend.fontsize': FONT_SIZE_PT, 'mathtext.fontset': 'cm',
        'axes.linewidth': 0.65, 'pdf.fonttype': 3,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    })
    return family


def load_data():
    path = HERE / 'wall_resolved_profile_data.json'
    provenance = json.loads((HERE / 'wall_resolved_profile_provenance.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == provenance['data_sha256']
    data = json.loads(path.read_text())
    assert len(data['M6Tw025_BL_profiles']['axes']) == 4
    assert len(data['Part4_LowH_M3_BL_profiles']['axes']) == 2
    for record in data.values():
        assert record['unmodified_source_rerender_pixel_exact']
        for axis in record['axes']:
            assert len(axis['lines']) == 3
            for line in axis['lines']:
                assert len(line['x']) == len(line['y'])
                assert np.all(np.isfinite(line['x']))
                assert np.all(np.isfinite(line['y']))
    return data


def make_figure(record, output_dir):
    four_panels = len(record['axes']) == 4
    height_mm = 144.0 if four_panels else 90.0
    fig, axes = plt.subplots(2 if four_panels else 1, 2,
                             figsize=(WIDTH_MM / 25.4, height_mm / 25.4),
                             constrained_layout=True)
    axes = axes.flat
    # An external shared legend avoids covering the enlarged plots. The lower
    # note retains the source's provisional DNS reconstruction qualification.
    layout_rect = (0, 0, 1, 0.92) if four_panels else (0, 0.125, 1, 0.76)
    fig.get_layout_engine().set(rect=layout_rect, w_pad=0.03, h_pad=0.035,
                               hspace=0.05, wspace=0.05)
    point_count = 0
    panel_audit = []
    for ax, original in zip(axes, record['axes']):
        for line in original['lines']:
            displayed, = ax.plot(line['x'], line['y'], label=line['label'],
                                  color=line['color'], linestyle=line['linestyle'],
                                  marker=line['marker'], markevery=line['markevery'],
                                  markerfacecolor=line['color'], markersize=3.0,
                                  linewidth=1.5, zorder=line['zorder'])
            # Cache values are round-tripped exactly and are used without any
            # resampling, smoothing, normalization or changed endpoint rule.
            np.testing.assert_array_equal(displayed.get_xdata(), line['x'])
            np.testing.assert_array_equal(displayed.get_ydata(), line['y'])
            point_count += len(line['x'])
        title = original['title']
        if title == 'Dimensional temperature-velocity relation':
            title = 'Temperature–velocity relation'
        ax.set(title=title, xlabel=original['xlabel'], ylabel=original['ylabel'],
               xscale=original['xscale'], yscale=original['yscale'],
               xlim=original['xlim'], ylim=original['ylim'])
        if original['xscale'] == 'log':
            ax.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=12))
            ax.xaxis.set_minor_locator(ticker.LogLocator(
                base=10, subs=np.arange(2, 10), numticks=100))
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax.grid(True, which='major', color='#D0D0D0', linewidth=0.55)
        ax.grid(True, which='minor', color='#ECECEC', linewidth=0.35)
        panel_audit.append({'original_title': original['title'],
                            'curves': len(original['lines']),
                            'points': sum(len(line['x']) for line in original['lines']),
                            'exact_cache_coordinates': True,
                            'original_axis_limits_retained': True})
    handles, labels = fig.axes[0].get_legend_handles_labels()
    legend = fig.legend(handles, labels, loc='upper center', ncol=3,
                         bbox_to_anchor=(0.5, 1.0), frameon=False,
                         handlelength=1.6, handletextpad=0.5, columnspacing=1.3,
                         borderaxespad=0.3)
    note = None
    if not four_panels:
        closure = abs(record['DNS_reconstruction_provisional_endpoint_closure_pct'])
        note = fig.text(0.5, 0.025,
                        r'DNS reconstruction assumed $\delta_{99}$ endpoint closure'
                        + f'\n(provisional; {closure:.2f}% mismatch).',
                        ha='center', va='bottom', fontsize=FONT_SIZE_PT)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_records = []
    visible_text = list(fig.texts) + list(legend.get_texts())
    for ax in fig.axes:
        visible_text.append(ax.title)
        for axis in (ax.xaxis, ax.yaxis):
            visible_text.extend((axis.label, axis.get_offset_text()))
            # Matplotlib creates additional ticks beyond the view limits but
            # does not draw them. Audit only ticks within the plotted range.
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if lo <= tick.get_loc() <= hi:
                    visible_text.extend((tick.label1, tick.label2))
    for artist in visible_text:
        if not artist.get_visible() or not artist.get_text():
            continue
        bbox = artist.get_window_extent(renderer)
        # Empty minor-tick labels have already been skipped. All visible
        # labels, panel titles, legend and notes must fit the fixed canvas.
        assert bbox.x0 >= -0.01 and bbox.y0 >= -0.01, (artist.get_text(), bbox)
        assert bbox.x1 <= fig.bbox.x1 + 0.01, (artist.get_text(), bbox)
        assert bbox.y1 <= fig.bbox.y1 + 0.01, (artist.get_text(), bbox)
        assert artist.get_fontsize() >= FONT_SIZE_PT
        text_records.append({'text': artist.get_text(), 'base_font_pt': artist.get_fontsize()})
    # Prevent the shared legend or qualification from touching any panel text.
    legend_box = legend.get_window_extent(renderer)
    note_box = note.get_window_extent(renderer) if note is not None else None
    for ax in fig.axes:
        panel_box = ax.get_tightbbox(renderer)
        assert not panel_box.overlaps(legend_box)
        if note_box is not None:
            assert not panel_box.overlaps(note_box)
    stem = output_dir / record['stem']
    fig.savefig(stem.with_suffix('.pdf'), metadata={
        'Creator': Path(__file__).name,
        'Title': f"{record['case']}: wall-resolved profiles compared with DNS"})
    fig.savefig(stem.with_suffix('.png'), dpi=300)
    plt.close(fig)
    return {'width_mm': WIDTH_MM, 'height_mm': height_mm,
            'base_font_pt': FONT_SIZE_PT, 'panels': panel_audit,
            'verified_exact_points': point_count,
            'all_visible_text_inside_canvas': True,
            'legend_and_note_clear_of_panels': True,
            'text': text_records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=HERE)
    parser.add_argument('--font-dir', type=Path, default=Path(
        '/mnt/c/texlive/2026/texmf-dist/fonts/opentype/public/tex-gyre'))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    family = configure_fonts(args.font_dir)
    data = load_data()
    audit = {'font_family': family, 'figures': {
        stem: make_figure(record, args.output_dir) for stem, record in data.items()}}
    (args.output_dir / 'wall_resolved_profile_layout_audit.json').write_text(
        json.dumps(audit, indent=2) + '\n')
    print('Replotted both profile figures; verified all 2592 retained plotted points exactly.')


if __name__ == '__main__':
    main()
