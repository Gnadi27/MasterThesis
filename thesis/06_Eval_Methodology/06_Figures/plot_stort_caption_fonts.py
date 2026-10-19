"""Regenerate Figure 6.2 at 170 mm width with caption-sized (10 pt) text.

Run with Python, NumPy, Matplotlib and openpyxl. The adjacent stort_sources
directory retains the original calculation script and source workbook.
Only display notation, typography and layout change; all plotted samples,
axis limits, and the source calculation CSV are checked against the source.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import runpy
import sys
import tempfile
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.text import Text
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCES = HERE / "stort_sources"
WIDTH_IN = 170 / 25.4
HEIGHT_IN = 4.05
FONT_PT = 10.0
ORIGINAL_CSV_SHA256 = "4d256a621c2d454720001061383d7d6cf9645dd73d4dc85afb5ba73e66362b92"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plotted_state(fig):
    return [{"xlim": list(ax.get_xlim()), "ylim": list(ax.get_ylim()),
             "xscale": ax.get_xscale(), "yscale": ax.get_yscale(),
             "curves": [{"label": line.get_label(),
                         "points": len(line.get_xdata()),
                         "sha256": hashlib.sha256(np.column_stack(
                             (line.get_xdata(), line.get_ydata())).astype("<f8").tobytes()).hexdigest()}
                        for line in ax.lines]} for ax in fig.axes]


def check_layout(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    visible = []
    for ax in fig.axes:
        visible.extend(ax.texts)
        for axis, bounds in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            if not axis.get_visible():
                continue
            visible.append(axis.label)
            for tick in axis.get_major_ticks():
                if bounds[0] <= tick.get_loc() <= bounds[1]:
                    visible.extend((tick.label1, tick.label2))
    for legend in fig.legends:
        visible.extend(legend.get_texts())
        box = legend.get_window_extent(renderer)
        assert box.y0 > fig.axes[0].get_window_extent(renderer).y1, "Legend overlaps plotting area"
    checked = []
    for text in visible:
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height, text.get_text()
        assert text.get_fontsize() == FONT_PT, text.get_text()
        checked.append(text.get_text())
    return checked


def main():
    script = SOURCES / "plot_stort_recovery.py"
    workbook = SOURCES / "STORT_experiment_data.xlsx"
    with tempfile.TemporaryDirectory(prefix="stort-fonts-") as temporary:
        argv = [str(script), "--workbook", str(workbook), "--output-dir", temporary]
        with patch.object(sys, "argv", argv), patch.object(Figure, "savefig"):
            source = runpy.run_path(str(script), run_name="__main__")
        csv_hash = digest(Path(temporary) / "STORT_wall_recovery_ratio.csv")
        assert csv_hash == ORIGINAL_CSV_SHA256, "Source calculation CSV changed"
    fig, ax = source["fig"], source["ax"]
    original = plotted_state(fig)

    # The manuscript's cropped source export displays the legend and plot,
    # without the source figure's title, equation or two file-detail footnotes.
    for text in list(fig.texts):
        text.remove()
    for legend in list(fig.legends):
        legend.remove()
    fig.set_size_inches(WIDTH_IN, HEIGHT_IN)
    fig.subplots_adjust(left=0.085, right=0.793, bottom=0.13, top=0.78)
    for axes in fig.axes:
        axes.tick_params(axis="both", labelsize=FONT_PT, pad=3)
        axes.yaxis.labelpad = 5
        for line in axes.lines:
            line.set_linewidth(1.25 if axes is ax else 1.0)
    ax.lines[-1].set_linewidth(0.7)
    fig.axes[2].spines["right"].set_position(("outward", 45))
    for text in fig.findobj(Text):
        text.set_text(text.get_text().replace(r"T_\mathrm{ad}", "T_r"))
        text.set_fontsize(FONT_PT)
    ax.set_ylabel(r"$T_w/T_r$", fontsize=FONT_PT)
    ax.texts[0].set_text(r"$T_w=T_r$")
    ax.texts[0].set_position((43, 1.025))

    handles, labels = ax.get_legend_handles_labels()
    trajectory_lines = source["trajectory_lines"]
    fig.legend(handles + trajectory_lines,
               labels + [line.get_label() for line in trajectory_lines],
               loc="upper left", bbox_to_anchor=(0.085, 0.995),
               borderaxespad=0, ncol=2, frameon=False, fontsize=FONT_PT,
               handlelength=2.2, handletextpad=0.6, columnspacing=1.2,
               labelspacing=0.45)
    final = plotted_state(fig)
    assert final == original, "Numerical plot content changed"
    visible_text = check_layout(fig)
    for extension in ("png", "pdf"):
        fig.savefig(HERE / f"STORT_Twall_Tad.{extension}", dpi=300)
    audit = {
        "original_script": "/home/mgnadlinger/ZPG_FP_overview/STORT/recovery_temperature/plot_stort_recovery.py",
        "original_workbook": "/home/mgnadlinger/ZPG_FP_overview/STORT/STORT_experiment_data.xlsx",
        "script_snapshot": str(script.relative_to(HERE)),
        "script_sha256": digest(script),
        "workbook_snapshot": str(workbook.relative_to(HERE)),
        "workbook_sha256": digest(workbook),
        "calculation_csv_sha256": csv_hash,
        "calculation_csv_identical_to_original": csv_hash == ORIGINAL_CSV_SHA256,
        "matplotlib_version": matplotlib.__version__,
        "figure_width_mm": 170,
        "figure_height_mm": HEIGHT_IN * 25.4,
        "font_size_pt_at_manuscript_width": FONT_PT,
        "curves_and_axes_unchanged": final == original,
        "visible_text_inside_canvas": visible_text,
        "plotted_state": final,
        "display_changes": "10 pt labels, ticks, annotation and legend at 170 mm manuscript width; T_r notation retained; layout and line widths adjusted for the smaller physical canvas.",
    }
    (SOURCES / "font_layout_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    plt.close(fig)
    print("Saved STORT_Twall_Tad.png and .pdf; plotted samples and axes unchanged.")


if __name__ == "__main__":
    main()
