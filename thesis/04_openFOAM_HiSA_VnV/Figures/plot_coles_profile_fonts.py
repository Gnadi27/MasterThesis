"""Regenerate the two Coles velocity figures with caption-sized plot text.

Uses the accepted iteration-20000 CSVs and the original NASA plotting
functions in the adjacent Files/ValData archive. The original WSL simulations
are /home/mgnadlinger/Coles_TurbZPGFP/{kwSST,SA}; no simulation is rerun.
Run with Python, NumPy, pandas, Matplotlib and TeX Live (Latin Modern).
Only typography, tick density and layout change; plotted coordinates and
axis limits are checked against the original plotting functions.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
from unittest.mock import patch

if os.name == "nt":
    os.environ["PATH"] = r"C:\texlive\2026\bin\windows" + os.pathsep + os.environ["PATH"]
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "coles-mpl-cache"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plotted_state(fig):
    return [{"xlim": list(ax.get_xlim()), "ylim": list(ax.get_ylim()),
             "xscale": ax.get_xscale(), "yscale": ax.get_yscale(),
             "curves": [hashlib.sha256(np.column_stack(
                 (line.get_xdata(), line.get_ydata())).astype("<f8").tobytes()).hexdigest()
                        for line in ax.lines]} for ax in fig.axes]


def check_layout(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text = list(fig.texts)
    for ax in fig.axes:
        text.extend((ax.xaxis.label, ax.yaxis.label))
        for axis, bounds in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            for tick in axis.get_major_ticks():
                if bounds[0] <= tick.get_loc() <= bounds[1]:
                    text.extend((tick.label1, tick.label2))
        legend = ax.get_legend()
        text.extend(legend.get_texts())
        box = legend.get_window_extent(renderer).padded(2)
        for line in ax.lines:
            if line.get_path().transformed(line.get_transform()).intersects_bbox(box, filled=False):
                raise ValueError(f"Legend overlaps curve: {line.get_label()}")
    for item in text:
        if item.get_visible() and item.get_text():
            box = item.get_window_extent(renderer)
            if (box.x0 < 0 or box.y0 < 0 or box.x1 > fig.bbox.width
                    or box.y1 > fig.bbox.height):
                raise ValueError(f"Text outside figure: {item.get_text()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=WORKSPACE.parents[1] /
                        "Files/ValData/Tunel/03 TurbZPGFP/01 Coles")
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    validation = args.data_root / "validation"
    script = validation / "compare_nasa_figures.py"
    spec = importlib.util.spec_from_file_location("coles_original", script)
    cmp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cmp)
    accepted = validation / "outputs_nasa_matched_20000"
    summary_path = accepted / "validation_summary.json"
    summary = json.loads(summary_path.read_text())
    assert summary["iteration"] == 20000 and summary["convergence_accepted"]
    cmp.U_INF = summary["freestream"]["U_m_s"]
    cmp.CASES = {
        "kwSST": {"label": "OpenFOAM/HiSA SST", "color": "#1f77b4", "marker": "o"},
        "SA": {"label": "OpenFOAM/HiSA SA", "color": "#2ca02c", "marker": "s"},
    }
    cmp.REFERENCES = {
        "Str-Uns_SST": {"label": "NASA Str-Uns SST", "case_dir": "Str-Uns_SST",
                        "color": "#9ecae1", "linestyle": "-"},
        "Str-Uns_SA": {"label": "NASA Str-Uns S-A", "case_dir": "Str-Uns_SA",
                       "color": "#98df8a", "linestyle": "-"},
    }
    refs = {key: cmp.load_reference(args.data_root, key) for key in cmp.REFERENCES}
    experiment, _, _, _ = cmp.load_experiment(args.data_root)
    cases = {}
    sources = [script, validation / "replot_nasa_matched.py", summary_path]
    for key in cmp.CASES:
        source = accepted / "data" / key / "profile_x0p5461.csv"
        profile = pd.read_csv(source)
        cases[key] = {"profile": profile[profile.y_m <= .03].copy()}
        sources.append(source)
    sources.extend(args.data_root / "NASA Downloads" / "Cases" / key / filename
                   for key in cmp.REFERENCES for filename in ("FP4p5.BL.dat", "FP4p5.Cf.dat"))
    sources.extend(args.data_root / "NASA Downloads" / "Results" / filename
                   for filename in ("FP4p5-BL.dat", "FP4p5-Cf.dat"))
    audit = {"original_simulation_root": "/home/mgnadlinger/Coles_TurbZPGFP",
             "source_sha256": {str(p): digest(p) for p in sources},
             "U_inf_m_s": cmp.U_INF, "figures": {}}
    args.output.mkdir(parents=True, exist_ok=True)
    savefig = Figure.savefig

    def save(fig, path, **kwargs):
        before = plotted_state(fig)
        profiles = Path(path).stem == "fig6_velocity_profiles"
        name = "06_Coles_velocity_profiles" if profiles else "06_Coles_velocity_experiment"
        # Preserve the old PDF page dimensions so the inclusion scales do not
        # change. At 481.89 pt and 0.85*481.89 pt manuscript widths, both
        # figures have approximately 10 PDF-point labels, ticks and legends.
        width, height = (751.68, 327.745) if profiles else (578.88, 387.361)
        font = 15.55 if profiles else 14.1
        fig.set_size_inches(width / 72, height / 72)
        if fig._suptitle:
            fig._suptitle.set_text("")
        for ax in fig.axes:
            ax.set_title("")
            if ax.get_xlabel() == "y [m]":
                ax.set_xlabel(r"$y$ [m]")
            if ax.get_ylabel() == "y [m]":
                ax.set_ylabel(r"$y$ [m]")
            ax.xaxis.label.set_fontsize(font)
            ax.yaxis.label.set_fontsize(font)
            ax.tick_params(axis="both", which="both", labelsize=font)
            ax.xaxis.get_offset_text().set_fontsize(font)
            ax.yaxis.get_offset_text().set_fontsize(font)
            ax.legend(loc="upper left" if profiles else "lower right", fontsize=font, borderpad=.4,
                      labelspacing=.35, handlelength=1.6, handletextpad=.5)
        if profiles:
            fig.axes[1].xaxis.set_major_locator(MultipleLocator(.002))
            fig.subplots_adjust(left=.105, right=.985, bottom=.19, top=.84, wspace=.13)
        else:
            fig.subplots_adjust(left=.12, right=.98, bottom=.155, top=.88)
        fig.text(.5, .965, r"Boundary-layer profiles at $x=0.5461$ m",
                 ha="center", va="top", fontsize=font)
        assert before == plotted_state(fig), "Plotted data or axis limits changed"
        check_layout(fig)
        for extension in ("pdf", "png"):
            savefig(fig, args.output / f"{name}.{extension}", dpi=220)
        audit["figures"][name] = {
            "font_size_pt": font, "page_size_pdf_pt": [width, height],
            "effective_font_size_pdf_pt": font * 481.89 * (1 if profiles else .85) / width,
            "plotted_state_unchanged": True, "plotted_state": before,
        }

    with patch.object(Figure, "savefig", save), patch.object(Figure, "tight_layout"):
        cmp.plot_fig6(args.output, refs, cases)
        cmp.plot_fig7(args.output, refs, cases, experiment)
    (args.output / "coles_profile_font_provenance.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value["effective_font_size_pdf_pt"]
                      for key, value in audit["figures"].items()}, indent=2))


if __name__ == "__main__":
    main()
